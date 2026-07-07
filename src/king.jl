function __king_fun!(resid, p, x)
    @inbounds @simd ivdep for i in eachindex(resid)
        resid[i] = p[1] + p[2] * x[2, i]^(p[3]) - x[1, i]^2
    end
    return nothing
end

# Common Solve Interface for KingCurveFitAlgorithm
@concrete struct KingFitCache <: AbstractCurveFitCache
    prob <: CurveFitProblem
    alg <: KingCurveFitAlgorithm
    kwargs
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::KingCurveFitAlgorithm;
        alias::CurveFitAliasSpecifier = CurveFitAliasSpecifier(), kwargs...
    )
    @assert !is_nonlinear_problem(prob) "King's law fitting doesn't work with nlfunc specification."
    @assert prob.u0 === nothing "King's law fit doesn't support initial guess (u0) specification"
    sigma_not_supported(prob)
    bounds_not_supported(prob)
    prob = _alias_inputs(prob, alias)
    return KingFitCache(prob, alg, kwargs)
end

function CommonSolve.solve!(cache::KingFitCache)
    # Fit sqrt(U) = b + a*E^2, then recover King's constants A and B:
    # sqrt(U) = -A/B + (1/B)*E^2  =>  B = 1/a, A = -b/a
    b, a = __linear_fit_internal(abs2, cache.prob.x, sqrt, cache.prob.y, nothing)
    B = 1 / a
    A = -b * B
    y_pred = ((cache.prob.x .^ 2 .- A) ./ B) .^ 2
    resid = cache.prob.y .- y_pred
    return CurveFitSolution(cache.alg, (A, B), resid, cache.prob, ReturnCode.Success)
end

function (sol::CurveFitSolution{<:KingCurveFitAlgorithm})(x)
    A, B = sol.u
    return @. ((x^2 - A) / B)^2
end

# Common Solve Interface for ModifiedKingCurveFitAlgorithm
@concrete struct ModifiedKingFitCache <: AbstractCurveFitCache
    initial_guess_cache <: Union{Nothing, KingFitCache}
    prob <: CurveFitProblem
    alg <: ModifiedKingCurveFitAlgorithm
    kwargs
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::ModifiedKingCurveFitAlgorithm;
        alias::CurveFitAliasSpecifier = CurveFitAliasSpecifier(), kwargs...
    )
    @assert !is_nonlinear_problem(prob) "Modified King's law fitting doesn't work with \
                                         nlfunc specification."
    sigma_not_supported(prob)

    prob = _alias_inputs(prob, alias)
    initial_guess_cache = if prob.u0 !== nothing
        nothing
    else
        # KingCurveFitAlgorithm doesn't support bounds so we need to remove them first
        nobounds_prob = if !isnothing(prob.lb) || !isnothing(prob.ub)
            x = @set prob.lb = nothing
            @set! x.ub = nothing
            x
        else
            prob
        end
        init(nobounds_prob, KingCurveFitAlgorithm(); kwargs...)
    end

    return ModifiedKingFitCache(initial_guess_cache, prob, alg, kwargs)
end

function CommonSolve.solve!(cache::ModifiedKingFitCache)
    if cache.initial_guess_cache !== nothing
        sol = solve!(cache.initial_guess_cache)
        u0 = [sol.u[1], sol.u[2], 0.5]
    else
        u0 = cache.prob.u0
    end

    # Re-create the nonlinear problem with the computed u0 to avoid reinit! issues
    # with the composite cache in NonlinearSolve.
    nonlinear_prob = NonlinearCurveFitProblem(
        NonlinearFunction{true}(
            __king_fun!;
            resid_prototype = similar(cache.prob.x)
        ),
        u0,
        stack((cache.prob.x, cache.prob.y); dims = 1),
        nothing;
        lb = cache.prob.lb, ub = cache.prob.ub
    )

    sol = solve(nonlinear_prob, __FallbackNonlinearFitAlgorithm(cache.alg.alg); cache.kwargs...)

    # The inner solve fits E² = A + B·Uⁿ, so its residual `A + B·Uⁿ − E²` is a
    # difference of E² values. We instead report `residuals` as a difference of
    # velocities, `U − Û`, so they match `fitted`/`sol(x)` (i.e.
    # `residuals == prob.y .- fitted`). `vcov` rebuilds the E² residual it needs
    # from `sol.u`.
    A, B, n = sol.u
    resid = @. cache.prob.y - ((cache.prob.x^2 - A) / B)^(1 / n)
    return CurveFitSolution(cache.alg, sol.u, resid, cache.prob, sol.retcode, sol.original)
end

function (sol::CurveFitSolution{<:ModifiedKingCurveFitAlgorithm})(x)
    return @. ((x^2 - sol.u[1]) / sol.u[2])^(1 / sol.u[3])
end
