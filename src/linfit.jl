function __linear_fit_internal(
        fnx::F1, x::AbstractArray{T1}, fny::F2, y::AbstractArray{T2}, sigma::Union{AbstractArray, Nothing}
    ) where {F1, F2, T1, T2}
    T = promote_type(T1, T2)

    sw, sx2, sxy, sx, sy = zero(T), zero(T), zero(T), zero(T), zero(T)
    @simd ivdep for i in eachindex(x, y)
        fn_xi = fnx(x[i])
        fn_yi = fny(y[i])

        w = isnothing(sigma) ? one(T) : inv(sigma[i]^2)
        sw += w

        sx += w * fn_xi
        sy += w * fn_yi
        sx2 = muladd(w * fn_xi, fn_xi, sx2)
        sxy = muladd(w * fn_xi, fn_yi, sxy)
    end

    det = (sw * sx2 - sx * sx)
    a0 = (sx2 * sy - sxy * sx) / det
    a1 = (sw * sxy - sx * sy) / det

    return (a0, a1)
end

function __vandermondepoly!(A, x, n)
    A[:, 1] .= 1
    @inbounds for i in 1:n
        @simd ivdep for k in axes(A, 1)
            A[k, i + 1] = A[k, i] * x[k]
        end
    end
    return
end

# Default Solver
@concrete struct GenericLinearFitCache <: AbstractCurveFitCache
    prob <: CurveFitProblem
    kwargs
    alg <: LinearCurveFitAlgorithm
end

function CommonSolve.init(prob::CurveFitProblem, alg::LinearCurveFitAlgorithm; kwargs...)
    @assert !is_nonlinear_problem(prob) "Linear curve fitting only works with linear \
                                         problems"
    @assert prob.u0 === nothing "Linear fit doesn't support initial guess (u0) \
                               specification"

    return GenericLinearFitCache(prob, kwargs, alg)
end

function CommonSolve.solve!(cache::GenericLinearFitCache)
    b, a = __linear_fit_internal(
        cache.alg.xfun, cache.prob.x, cache.alg.yfun, cache.prob.y, cache.prob.sigma
    )
    y_pred = cache.alg.yfun_inverse.(b .+ a .* cache.alg.xfun.(cache.prob.x))
    resid = cache.prob.y .- y_pred
    return CurveFitSolution(
        cache.alg, (a, b), resid, cache.prob, ReturnCode.Success
    )
end

function (sol::CurveFitSolution{<:LinearCurveFitAlgorithm})(x)
    a, b = sol.u
    return sol.alg.yfun_inverse.(b .+ a .* sol.alg.xfun.(x))
end

# Polynomial Fit
@concrete struct PolynomialFitCache <: AbstractCurveFitCache
    vandermondepoly_cache <: AbstractMatrix
    linsolve_cache
    prob <: CurveFitProblem
    alg <: PolynomialFitAlgorithm
    kwargs
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::PolynomialFitAlgorithm; kwargs...
    )
    @assert !is_nonlinear_problem(prob) "Linear curve fitting only works with linear \
                                         problems"
    @assert prob.u0 === nothing "Polynomial fit doesn't support initial guess \
                               (u0) specification"
    sigma_not_supported(prob)

    vandermondepoly_cache = similar(prob.x, length(prob.x), alg.degree + 1)
    linsolve_cache = init(
        LinearProblem(vandermondepoly_cache, prob.y), alg.linsolve_algorithm; kwargs...
    )
    return PolynomialFitCache(vandermondepoly_cache, linsolve_cache, prob, alg, kwargs)
end

function CommonSolve.solve!(cache::PolynomialFitCache)
    __vandermondepoly!(cache.vandermondepoly_cache, cache.prob.x, cache.alg.degree)
    cache.linsolve_cache.A = cache.vandermondepoly_cache
    sol = solve!(cache.linsolve_cache)
    # Always compute residuals manually as LinearSolve may return incorrect residuals
    resid = cache.prob.y .- cache.vandermondepoly_cache * sol.u
    return CurveFitSolution(cache.alg, sol.u, resid, cache.prob, sol.retcode)
end

function (sol::CurveFitSolution{<:PolynomialFitAlgorithm})(x)
    return evalpoly.(x, Ref(sol.u))
end
