@concrete struct NonlinearFunctionWrapper{iip}
    target
    f
end

SciMLBase.isinplace(::NonlinearFunctionWrapper{iip}) where {iip} = iip

__wrap_nonlinear_function(f::NonlinearFunction, ::Nothing) = f
function __wrap_nonlinear_function(f::NonlinearFunction, target)
    internal_f = NonlinearFunctionWrapper{SciMLBase.isinplace(f)}(target, f.f)
    @set! f.f = internal_f
    @set! f.resid_prototype = similar(target)
    return f
end

(nlf::NonlinearFunctionWrapper{false})(p, X) = nlf.f(p, X) .- nlf.target

function (nlf::NonlinearFunctionWrapper{true})(resid, p, X)
    nlf.f(resid, p, X)
    resid .-= nlf.target
    return resid
end

# NLLS Solvers
@concrete struct GenericNonlinearCurveFitCache <: AbstractCurveFitCache
    prob <: CurveFitProblem
    cache
    alg
    kwargs
end

function SciMLBase.reinit!(cache::GenericNonlinearCurveFitCache; u0 = nothing, x = nothing, y = nothing, kwargs...)
    if !isnothing(u0)
        kwargs = (; kwargs..., u0)
    end

    # x becomes `p` (parameter) in the NonlinearLeastSquaresProblem
    if !isnothing(x)
        kwargs = (; kwargs..., p = x)
    end

    # Update `y` inplace, which is stored in NonlinearFunctionWrapper.target
    if !isnothing(y)
        wrapper = cache.cache.prob.f.f
        copyto!(wrapper.target, y)
    end

    if !isempty(kwargs)
        reinit!(cache.cache; kwargs...)
    end

    return cache
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::__FallbackNonlinearFitAlgorithm; kwargs...
    )
    @assert is_nonlinear_problem(prob) "Nonlinear curve fitting only works with nonlinear \
                                        problems"
    @assert prob.u0 !== nothing "Nonlinear curve fitting requires an initial guess (u0)"

    return GenericNonlinearCurveFitCache(
        prob,
        init(
            NonlinearLeastSquaresProblem(
                __wrap_nonlinear_function(prob.nlfunc, prob.y), prob.u0, prob.x
            ),
            alg.alg;
            kwargs...
        ),
        alg,
        kwargs
    )
end

function CommonSolve.solve!(cache::GenericNonlinearCurveFitCache)
    sol = solve!(cache.cache)
    return CurveFitSolution(cache.alg, sol.u, sol.resid, cache.prob, sol.retcode, sol)
end

function (sol::CurveFitSolution{<:__FallbackNonlinearFitAlgorithm})(x)
    return sol.prob.nlfunc(sol.u, x)
end
