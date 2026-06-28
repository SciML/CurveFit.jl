@concrete struct NonlinearFunctionWrapper{iip}
    target
    sigma
    f
end

SciMLBase.isinplace(::NonlinearFunctionWrapper{iip}) where {iip} = iip

_unwrap_nonlinear_function(f::NonlinearFunctionWrapper) = f
_unwrap_nonlinear_function(f::NonlinearSolveBase.AutoSpecializeCallable) = _unwrap_nonlinear_function(f.orig)
_unwrap_nonlinear_function(f::NonlinearSolveBase.BoundedWrapper) = _unwrap_nonlinear_function(f.f.f)
_unwrap_nonlinear_function(f) = f

# If `target` is nothing then we can completely ignore sigma
__wrap_nonlinear_function(f::NonlinearFunction, ::Nothing, _) = f
function __wrap_nonlinear_function(f::NonlinearFunction, target, sigma)
    internal_f = NonlinearFunctionWrapper{SciMLBase.isinplace(f)}(target, sigma, f.f)
    @set! f.f = internal_f
    @set! f.resid_prototype = similar(target)
    return f
end

# Out-of-place
function (nlf::NonlinearFunctionWrapper{false})(p, X)
    resid = nlf.target .- nlf.f(p, X)

    if !isnothing(nlf.sigma)
        resid ./= nlf.sigma
    end

    return resid
end

# In-place
function (nlf::NonlinearFunctionWrapper{true})(resid, p, X)
    nlf.f(resid, p, X)
    resid .= nlf.target .- resid

    if !isnothing(nlf.sigma)
        resid ./= nlf.sigma
    end

    return resid
end

# NLLS Solvers
@concrete struct GenericNonlinearCurveFitCache <: AbstractCurveFitCache
    prob <: CurveFitProblem
    cache
    u0
    alg
    kwargs
end

function SciMLBase.reinit!(cache::GenericNonlinearCurveFitCache; u0 = nothing, x = nothing, y = nothing, sigma = nothing, kwargs...)
    if !isnothing(u0)
        kwargs = (; kwargs..., u0)
        copyto!(cache.u0, u0)
    end

    # x becomes `p` (parameter) in the NonlinearLeastSquaresProblem
    if !isnothing(x)
        kwargs = (; kwargs..., p = x)
    end

    # Update `y` inplace
    wrapper = _unwrap_nonlinear_function(cache.cache.prob.f.f)
    if !isnothing(y)
        copyto!(wrapper.target, y)
    end

    # Update `sigma` inplace
    if !isnothing(sigma)
        copyto!(wrapper.sigma, sigma)
    end

    reinit!(cache.cache; kwargs...)

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
                __wrap_nonlinear_function(prob.nlfunc, prob.y, prob.sigma), prob.u0, prob.x;
                lb = prob.lb, ub = prob.ub
            ),
            alg.alg;
            kwargs...
        ),
        copy(prob.u0),
        alg,
        kwargs
    )
end

function CommonSolve.solve!(cache::GenericNonlinearCurveFitCache)
    inner = _get_cache(cache)
    x = inner.p
    sol = solve!(cache.cache)

    y = cache.prob.y
    sigma = cache.prob.sigma

    wrapped_f = _unwrap_nonlinear_function(inner.prob.f.f)
    if wrapped_f isa NonlinearFunctionWrapper
        y = wrapped_f.target
        sigma = wrapped_f.sigma
    end

    # Reconstruct the problem with the current settings. We can't copy
    # cache.prob because the cache may have been reinit()'d in which case
    # cache.prob will be out of date and will give wrong results for the stats
    # functions that use it.
    prob = CurveFitProblem(
        x,
        y,
        sigma,
        cache.prob.nlfunc,
        cache.u0,
        cache.prob.lb,
        cache.prob.ub
    )
    return CurveFitSolution(cache.alg, sol.u, sol.resid, prob, sol.retcode, sol)
end

function (sol::CurveFitSolution{<:__FallbackNonlinearFitAlgorithm})(x)
    return sol.prob.nlfunc(sol.u, x)
end

function _get_cache(cache::GenericNonlinearCurveFitCache)
    inner = cache.cache
    return if inner isa NonlinearSolveBase.NonlinearSolvePolyAlgorithmCache
        inner.caches[inner.current]
    else
        inner
    end
end

function Base.show(io::IO, ::MIME"text/plain", cache::GenericNonlinearCurveFitCache)
    inner = cache.cache
    is_polyalg = inner isa NonlinearSolveBase.NonlinearSolvePolyAlgorithmCache
    current_cache = _get_cache(cache)

    context = (:compact => true, :limit => true)

    println(io, "GenericNonlinearCurveFitCache(")

    algstr = if !isnothing(current_cache.alg)
        NonlinearSolveBase.Utils.clean_sprint_struct(current_cache.alg, 4)
    else
        "nothing"
    end
    print(io, "    alg = ")
    if is_polyalg
        print(io, "[NonlinearSolvePolyAlgorithm] ")
    end
    println(io, algstr, ",")

    # Current parameter values
    ustr = sprint(show, current_cache.u; context)
    println(io, "    u = ", ustr, ",")

    # Residual
    resids = NonlinearSolveBase.get_fu(current_cache)
    residstr = sprint(show, resids; context)
    println(io, "    residual = ", residstr, ",")

    # Inf-norm of residual
    normval = LinearAlgebra.norm(resids, Inf)
    normstr = sprint(show, normval; context)
    println(io, "    inf-norm(residual) = ", normstr, ",")

    # Number of steps
    println(io, "    nsteps = ", inner.stats.nsteps, ",")

    # Return code
    println(io, "    retcode = ", current_cache.retcode)
    print(io, ")")
    return nothing
end
