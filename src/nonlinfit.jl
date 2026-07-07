@concrete struct NonlinearFunctionWrapper{iip}
    target
    sigma
    f
end

SciMLBase.isinplace(::NonlinearFunctionWrapper{iip}) where {iip} = iip

__wrap_nonlinear_function(f::NonlinearFunction, ::Nothing, ::Nothing) = f
function __wrap_nonlinear_function(f::NonlinearFunction, target, sigma)
    internal_f = NonlinearFunctionWrapper{SciMLBase.isinplace(f)}(target, sigma, f.f)
    @set! f.f = internal_f
    @set! f.resid_prototype = similar(isnothing(target) ? sigma : target)
    return f
end

# Out-of-place
function (nlf::NonlinearFunctionWrapper{false})(p, X)
    # `target` is only nothing for residual-style problems (y = nothing), in
    # which case the wrapper exists purely to apply `sigma`.
    if isnothing(nlf.target)
        return nlf.f(p, X) ./ nlf.sigma
    end

    resid = nlf.target .- nlf.f(p, X)

    if !isnothing(nlf.sigma)
        resid ./= nlf.sigma
    end

    return resid
end

# In-place
function (nlf::NonlinearFunctionWrapper{true})(resid, p, X)
    nlf.f(resid, p, X)

    if !isnothing(nlf.target)
        resid .= nlf.target .- resid
    end
    if !isnothing(nlf.sigma)
        resid ./= nlf.sigma
    end

    return resid
end

# NLLS Solvers
@concrete struct GenericNonlinearCurveFitCache <: AbstractCurveFitCache
    prob <: CurveFitProblem
    cache
    x <: AbstractArray
    y <: Union{AbstractArray, Nothing}
    sigma <: Union{AbstractArray, Nothing}
    u0
    alg
    kwargs
end

function SciMLBase.reinit!(cache::GenericNonlinearCurveFitCache; u0 = nothing, x = nothing, y = nothing, sigma = nothing, kwargs...)
    if !isnothing(u0)
        @assert size(u0) == size(cache.u0) "reiniting `u0` must keep the same size"
        kwargs = (; kwargs..., u0)
        copyto!(cache.u0, u0)
    end

    # x becomes `p` (parameter) in the NonlinearLeastSquaresProblem
    if !isnothing(x)
        @assert size(x) == size(cache.x) "reiniting `x` must keep the same size"
        copyto!(cache.x, x)
        kwargs = (; kwargs..., p = cache.x)
    end

    # Update `y` inplace
    if !isnothing(y)
        @assert !isnothing(cache.y) "cannot reinit `y` for a problem created without a `y`"
        @assert size(y) == size(cache.y) "reiniting `y` must keep the same size"
        copyto!(cache.y, y)
    end

    # Update `sigma` inplace
    if !isnothing(sigma)
        @assert !isnothing(cache.sigma) "cannot reinit `sigma` for a problem created without a `sigma`"
        @assert size(sigma) == size(cache.sigma) "reiniting `sigma` must keep the same size"
        copyto!(cache.sigma, sigma)
    end

    reinit!(cache.cache; kwargs...)

    return cache
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::__FallbackNonlinearFitAlgorithm;
        alias::CurveFitAliasSpecifier = CurveFitAliasSpecifier(), kwargs...
    )
    @assert is_nonlinear_problem(prob) "Nonlinear curve fitting only works with nonlinear \
                                        problems"
    @assert prob.u0 !== nothing "Nonlinear curve fitting requires an initial guess (u0)"

    # By default the input arrays are copied into the cache so that reinit!'ing
    # the cache (which mutates them in place) doesn't overwrite the user's data.
    prob = _alias_inputs(
        prob, alias; default_x = false, default_y = false, default_sigma = false
    )

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
        prob.x,
        prob.y,
        prob.sigma,
        prob.u0,
        alg,
        kwargs
    )
end

function CommonSolve.solve!(cache::GenericNonlinearCurveFitCache)
    sol = solve!(cache.cache)

    # Reconstruct the problem with the current settings. We can't use
    # cache.prob because the cache may have been reinit()'d in which case
    # cache.prob will be out of date and will give wrong results for the stats
    # functions that use it.
    prob = CurveFitProblem(
        cache.x,
        cache.y,
        cache.sigma,
        cache.prob.nlfunc,
        copy(cache.u0),
        cache.prob.lb,
        cache.prob.ub
    )
    return CurveFitSolution(cache.alg, sol.u, sol.resid, prob, sol.retcode, sol)
end

function (sol::CurveFitSolution{<:__FallbackNonlinearFitAlgorithm})(x)
    if SciMLBase.isinplace(sol.prob) && x isa AbstractArray
        out = similar(sol.resid, length(x))
        sol.prob.nlfunc(out, sol.u, x)
        return out
    end
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
