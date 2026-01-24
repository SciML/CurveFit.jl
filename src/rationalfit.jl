@concrete struct RationalPolynomial
    numerator <: AbstractVector
    denominator <: AbstractVector
end

function (rpoly::RationalPolynomial)(x::Number)
    return evalpoly(x, rpoly.numerator) / evalpoly(x, rpoly.denominator)
end

function __linear_rational_matrix!(A, x, y, p, q)
    @inbounds for i in axes(x, 1)
        A[i, 1] = true
        @simd ivdep for k in 1:p
            A[i, k + 1] = x[i]^k
        end
        @simd ivdep for k in 1:q
            A[i, p + 1 + k] = -y[i] * x[i]^k
        end
    end
    return
end

function __rational_fit_residual!(p::Integer, q::Integer)
    return let p = p, q = q
        (resid, coeffs, x) -> __rational_fit_residual!(resid, coeffs, x, p, q)
    end
end

function __rational_fit_residual!(resid, coeffs, x, p::Integer, q::Integer)
    num = view(coeffs, 1:(p + 1))
    den = vcat(one(eltype(x)), view(coeffs, (p + 2):(p + q + 1)))

    @inbounds @simd ivdep for i in eachindex(resid)
        resid[i] = evalpoly(x[i], num) / evalpoly(x[i], den)
    end

    return resid
end

# Common Solve Interface
@concrete struct LinearRationalFitCache <: AbstractCurveFitCache
    mat <: AbstractMatrix
    linsolve_cache
    prob <: CurveFitProblem
    alg <: RationalPolynomialFitAlgorithm
    kwargs
end

@concrete struct NonlinearRationalFitCache <: AbstractCurveFitCache
    initial_guess_cache <: Union{Nothing, LinearRationalFitCache}
    nonlinear_cache
    prob <: CurveFitProblem
    alg <: RationalPolynomialFitAlgorithm
    kwargs
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::RationalPolynomialFitAlgorithm; kwargs...
    )
    @assert !is_nonlinear_problem(prob) "Rational polynomial fitting doesn't work with \
                                         nlfunc specification."
    sigma_not_supported(prob)

    coeffs_length = alg.num_degree + alg.den_degree + 1

    if alg.alg isa AbstractLinearAlgorithm
        @assert prob.u0 === nothing "Rational polynomial fit doesn't support initial \
                                   guess (u0) specification"

        A = similar(prob.x, length(prob.x), coeffs_length)
        return LinearRationalFitCache(
            A, init(LinearProblem(A, prob.y), alg.alg; kwargs...), prob, alg, kwargs
        )
    end

    initial_guess_cache = if prob.u0 !== nothing
        nothing
    else
        A = similar(prob.x, length(prob.x), coeffs_length)
        LinearRationalFitCache(
            A, init(LinearProblem(A, prob.y), alg.alg; kwargs...), prob, alg, kwargs
        )
    end
    nonlinear_cache = init(
        NonlinearCurveFitProblem(
            NonlinearFunction{true}(
                __rational_fit_residual!(alg.num_degree, alg.den_degree);
                resid_prototype = similar(prob.x)
            ),
            similar(prob.x, coeffs_length),
            prob.x,
            prob.y
        ),
        __FallbackNonlinearFitAlgorithm(alg.alg);
        kwargs...
    )
    return NonlinearRationalFitCache(
        initial_guess_cache, nonlinear_cache, prob, alg, kwargs
    )
end

function CommonSolve.solve!(cache::LinearRationalFitCache)
    __linear_rational_matrix!(
        cache.mat, cache.prob.x, cache.prob.y, cache.alg.num_degree, cache.alg.den_degree
    )
    cache.linsolve_cache.A = cache.mat
    sol = solve!(cache.linsolve_cache)
    resid = if sol.resid === nothing
        # Linear problem: y ≈ p/q => y*q - p ≈ 0 (linearized residual)
        # But StatsAPI expects y - p/q (nonlinear residual) or linearized?
        # Standard definition is y - y_pred.
        # So we should compute y - p(x)/q(x) using the fitted params.

        # We need to construct the RationalPolynomial to eval it.
        # Helper function from rationalfit.jl isn't easily accessible inside solve! without alloc.
        # But we can reuse the logic from call:
        rpoly = RationalPolynomial(
            view(sol.u, 1:(cache.alg.num_degree + 1)),
            vcat(one(eltype(sol.u)), view(sol.u, (cache.alg.num_degree + 2):length(sol.u)))
        )
        cache.prob.y .- rpoly.(cache.prob.x)
    else
        sol.resid
    end
    return CurveFitSolution(cache.alg, sol.u, resid, cache.prob, sol.retcode)
end

function CommonSolve.solve!(cache::NonlinearRationalFitCache)
    if cache.initial_guess_cache !== nothing
        sol = solve!(cache.initial_guess_cache)
        u0 = sol.u # Linear fit returns all coefficients [p0...pn, q1...qn]
    else
        u0 = cache.prob.u0
    end

    # Re-create problem to avoid reinit! crash
    nonlinear_prob = NonlinearCurveFitProblem(
        NonlinearFunction{true}(
            __rational_fit_residual!(cache.alg.num_degree, cache.alg.den_degree);
            resid_prototype = similar(cache.prob.x)
        ),
        u0,
        cache.prob.x,
        cache.prob.y
    )

    sol = solve(nonlinear_prob, __FallbackNonlinearFitAlgorithm(cache.alg.alg); cache.kwargs...)

    return CurveFitSolution(
        cache.alg, sol.u, sol.resid, cache.prob, sol.retcode, sol.original
    )
end

function (sol::CurveFitSolution{<:RationalPolynomialFitAlgorithm})(x::Number)
    return RationalPolynomial(
        view(sol.u, 1:(sol.alg.num_degree + 1)),
        vcat(
            one(eltype(sol.u)),
            view(sol.u, (sol.alg.num_degree + 2):(length(sol.u)))
        )
    )(x)
end
