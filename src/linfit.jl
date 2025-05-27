function __linear_fit_internal(
        fnx::F1, x::AbstractArray{T1}, fny::F2, y::AbstractArray{T2}
) where {F1, F2, T1, T2}
    T = promote_type(T1, T2)
    m = length(x)

    sx2, sy2, sxy, sx, sy = zero(T), zero(T), zero(T), zero(T), zero(T)
    @simd ivdep for i in eachindex(x, y)
        fn_xi = fnx(x[i])
        fn_yi = fny(y[i])
        sx += fn_xi
        sy += fn_yi
        sx2 = muladd(fn_xi, fn_xi, sx2)
        sy2 = muladd(fn_yi, fn_yi, sy2)
        sxy = muladd(fn_xi, fn_yi, sxy)
    end

    a0 = (sx2 * sy - sxy * sx) / (m * sx2 - sx * sx)
    a1 = (m * sxy - sx * sy) / (m * sx2 - sx * sx)

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
@concrete struct GenericLinearFitCache
    prob <: LinearCurveFitProblem
    kwargs
end

function CommonSolve.init(prob::LinearCurveFitProblem, ::Nothing; kwargs...)
    return GenericLinearFitCache(prob, kwargs)
end

function CommonSolve.solve!(cache::GenericLinearFitCache)
    b, a = __linear_fit_internal(
        cache.prob.xfun, cache.prob.x, cache.prob.yfun, cache.prob.y
    )
    return LinearCurveFitSolution(nothing, (a, b), cache.prob)
end

function (sol::LinearCurveFitSolution{<:Nothing})(x::Number)
    a, b = sol.coeffs
    return sol.prob.yfun_inverse(b + a * sol.prob.xfun(x))
end

# Polynomial Fit
@concrete struct PolynomialFitCache
    vandermondepoly_cache <: AbstractMatrix
    linsolve_cache
    prob <: LinearCurveFitProblem
    alg <: PolynomialFitAlgorithm
    kwargs
end

function CommonSolve.init(
        prob::LinearCurveFitProblem, alg::PolynomialFitAlgorithm; kwargs...
)
    @assert prob.xfun===identity "Polynomial fit only works with LinearCurveFitProblem \
                                  with xfun = identity"
    @assert prob.yfun===identity "Polynomial fit only works with LinearCurveFitProblem \
                                  with yfun = identity"
    vandermondepoly_cache = similar(prob.x, length(prob.x), alg.degree + 1)
    linear_problem = LinearProblem(vandermondepoly_cache, prob.y)
    linsolve_cache = init(linear_problem, alg.linsolve_algorithm; kwargs...)
    return PolynomialFitCache(vandermondepoly_cache, linsolve_cache, prob, alg, kwargs)
end

function CommonSolve.solve!(cache::PolynomialFitCache)
    __vandermondepoly!(cache.vandermondepoly_cache, cache.prob.x, cache.alg.degree)
    cache.linsolve_cache.A = cache.vandermondepoly_cache
    sol = solve!(cache.linsolve_cache)
    return LinearCurveFitSolution(cache.alg, sol.u, cache.prob)
end

function (sol::LinearCurveFitSolution{<:PolynomialFitAlgorithm})(x::Number)
    return evalpoly(x, sol.coeffs)
end
