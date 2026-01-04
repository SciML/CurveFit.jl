function __king_fun!(resid, p, x)
    @inbounds @simd ivdep for i in eachindex(resid)
        resid[i] = p[1] + p[2] * x[2, i]^(p[3]) - x[1, i]^2
    end
    return nothing
end

# Common Solve Interface
@concrete struct ModifiedKingFitCache <: AbstractCurveFitCache
    initial_guess_cache <: Union{Nothing, GenericLinearFitCache}
    nonlinear_cache
    prob <: CurveFitProblem
    alg <: ModifiedKingCurveFitAlgorithm
    kwargs
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::ModifiedKingCurveFitAlgorithm; kwargs...
    )
    @assert !is_nonlinear_problem(prob) "Modified King's law fitting doesn't work with \
                                         nlfunc specification."

    initial_guess_cache = if prob.u0 !== nothing
        nothing
    else
        init(prob, KingCurveFitAlgorithm(); kwargs...)
    end

    nonlinear_cache = init(
        NonlinearCurveFitProblem(
            NonlinearFunction{true}(
                __king_fun!;
                resid_prototype = similar(prob.x)
            ),
            similar(prob.x, 3),
            stack((prob.x, prob.y); dims = 1),
            nothing
        ),
        __FallbackNonlinearFitAlgorithm(alg.alg);
        kwargs...
    )
    return ModifiedKingFitCache(initial_guess_cache, nonlinear_cache, prob, alg, kwargs)
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
        nothing
    )

    sol = solve(nonlinear_prob, __FallbackNonlinearFitAlgorithm(cache.alg.alg); cache.kwargs...)
    return CurveFitSolution(cache.alg, sol.u, sol.resid, cache.prob, sol.retcode, sol.original)
end

function (sol::CurveFitSolution{<:ModifiedKingCurveFitAlgorithm})(x::Number)
    return ((x .^ 2 .- sol.u[1]) ./ sol.u[2]) .^ (1 ./ sol.u[3])
end
