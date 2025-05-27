@concrete struct RationalPolynomial
    numerator <: AbstractVector
    denominator <: AbstractVector
end

function (rpoly::RationalPolynomial)(x::Number)
    return evalpoly(x, rpoly.numerator) / evalpoly(x, rpoly.denominator)
end

# Common Solve Interface
@concrete struct LinearRationalFitCache
    mat <: AbstractMatrix
    linsolve_cache
    prob <: CurveFitProblem
    alg <: RationalPolynomialFitAlgorithm
    kwargs
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::RationalPolynomialFitAlgorithm; kwargs...
)
    @assert !is_nonlinear_problem(prob) "Rational polynomial fitting doesn't work with \
                                         nlfunc specification."
    @assert prob.xfun===identity "Rational polynomial fit only works with \
                                  xfun = identity"
    @assert prob.yfun===identity "Rational polynomial fit only works with \
                                  yfun = identity"

    if alg.alg isa AbstractLinearAlgorithm
        A = similar(prob.x, length(prob.x), alg.num_degree + alg.den_degree + 1)
        return LinearRationalFitCache(
            A, init(LinearProblem(A, prob.y), alg.alg; kwargs...), prob, alg, kwargs
        )
    end

    error("TODO: Nonlinear rational polynomial fitting is not implemented yet")
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

function CommonSolve.solve!(cache::LinearRationalFitCache)
    __linear_rational_matrix!(
        cache.mat, cache.prob.x, cache.prob.y, cache.alg.num_degree, cache.alg.den_degree
    )
    cache.linsolve_cache.A = cache.mat
    sol = solve!(cache.linsolve_cache)
    return CurveFitSolution(cache.alg, sol.u, cache.prob, sol.retcode)
end

function (sol::CurveFitSolution{<:RationalPolynomialFitAlgorithm})(x::Number)
    return RationalPolynomial(
        view(sol.coeffs, 1:(sol.alg.num_degree + 1)),
        vcat(
            one(eltype(sol.coeffs)),
            view(sol.coeffs, (sol.alg.num_degree + 2):(length(sol.coeffs)))
        )
    )(x)
end

# """
# Auxiliary function used in nonlinear least squares
# """
# function make_rat_fun(p, q)
#     return let p = p, q = q
#         (y, x, a) -> begin
#             num = view(a, 1:(p + 1))
#             den = vcat(one(eltype(x)), view(a, (p + 2):(p + q + 1)))

#             @inbounds @simd ivdep for i in eachindex(y)
#                 y[i] = evalpoly(x[1, i], num) / evalpoly(x[1, i], den) - x[2, i]
#             end

#             return y
#         end
#     end
# end

# """
# # Carry out a nonlinear least squares of rational polynomials

# Find the polynomial coefficients that best approximate
# the points given by `x` and `y`.
# """
# function rational_fit(x, y, p, q, args...; kwargs...)
#     coefs0 = linear_rational_fit(x, y, p, q)
#     sol = nonlinear_fit(
#         make_rat_fun(p, q), stack((x, y); dims = 1), coefs0, args...;
#         resid_prototype = Vector{eltype(coefs0)}(undef, length(x)), iip = Val(true), kwargs...
#     )
#     return sol.u
# end
