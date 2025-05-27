# struct NonlinearFunctionWrapper{F, T}
#     f::F
#     target::T
# end

# (nlf::NonlinearFunctionWrapper{F, Nothing})(p, X) where {F} = nlf.f(X, p)
# (nlf::NonlinearFunctionWrapper{F, T})(p, X) where {F, T} = nlf.f(X, p) .- nlf.target
# function (nlf::NonlinearFunctionWrapper{F, Nothing})(dp, p, X) where {F}
#     nlf.f(dp, X, p)
#     return dp
# end
# function (nlf::NonlinearFunctionWrapper{F, T})(dp, p, X) where {F, T}
#     nlf.f(dp, X, p)
#     dp .-= nlf.target
#     return dp
# end

# @doc doc"""
#     nonlinear_fit(
#         f::F, data, p0, alg = nothing;
#         target = nothing, iip = nothing, solve_kwargs = (;), kwargs...
#     ) where {F}

# Nonlinear least squares fitting of data. This function is a wrapper around
# [`NonlinearLeastSquaresProblem`](@ref) and [`solve`](@ref) from the
# [`NonlinearSolve.jl`](https://github.com/SciML/NonlinearSolve.jl) package. We are fitting
# the function `f` to the data `data` using the initial guess `p0`.

# ```math
# \begin{equation}
#     \underset{p}{\text{argmin}} ~ \| f(\text{data}, p) - \text{target} \|_2
# \end{equation}
# ```

# If `target` is not provided, then it is treated as a zero vector.

# ## Parameters

#  * `f` the function that should be fitted
#  * `data` the data to be fitted
#  * `p0` the initial guess for the parameters
#  * `alg` the algorithm to be used for solving the problem. If not provided, the default
#    poly-algorithm is used.

# ## Keyword arguments

#  * `iip` determines if the function is inplace or not. If not provided, we automatically
#    determine this.
#  * `target` the target value to be fitted. If not provided, it is treated as a zero vector.
#  * `kwargs` are passed directly to the [`NonlinearLeastSquaresProblem`](@ref) constructor.
#  * `solve_kwargs` are passed directly to the `solve` function. See the
#    [documentation](https://docs.sciml.ai/NonlinearSolve/stable/basics/solve/) for more
#    details.

# !!! warning

#     The return is type unstable if `iip` is not provided.

# ## Return values

# A NonlinearSolution object is returned.
# """
# function nonlinear_fit(
#         f::F,
#         data,
#         p0,
#         alg = nothing;
#         iip::Union{Nothing, Val{true}, Val{false}} = nothing,
#         target = nothing,
#         solve_kwargs = (;),
#         kwargs...
# ) where {F}
#     iip === nothing && (iip = any(==(3), SciMLBase.numargs(f)))

#     return solve(
#         NonlinearLeastSquaresProblem(
#             NonlinearFunction{SciMLBase._unwrap_val(iip)}(
#                 NonlinearFunctionWrapper(f, target);
#                 kwargs...
#             ),
#             p0,
#             data
#         ),
#         alg;
#         solve_kwargs...
#     )
# end

# """
#    a = secant_nls_fit(x, y, fun, a0[[, eps,] maxiter])

# Secant/Gauss-Newton nonlinear least squares. DOESN'T NEED A DERIVATIVE FUNCTION. Given vectors `x` and `y`, the tries to fit parameters `a` to
# a function `f` using least squares approximation:

#  ``y = f(x, a₁, ..., aₙ)``

# For more general approximations, see [`gauss_newton_fit`](@ref).

# ### Arguments:

#  * `x` Vector with x values
#  * `y` Vector with y values
#  * `fun` a function that is called as `fun(x, a)` where `a` is a vector of parameters.
#  * `a0` Vector with the initial guesses of the parameters
#  * `eps` Maximum residpal for convergence
#  * `maxiter` Maximum number of iterations for convergence

# ## Return value

# A vector with the convrged array. If no convergence is achieved, the function throws an error.

# ## Specification of the fitting function

# The function that should be fitted shoud be specified by Julia funcion with the following signature:

# ```julia
# fun(x::T, a::AbstractVector{T}) where {T<:Number}
# ```

# ## Initial approximation (guess)

# If the initial approximation is not good enough, divergence is possible.

# **Careful** with parameters close to 0. The initial guess should never be 0.0 because the
# initial value of the parameter is used as reference value for computing resiudpals.

# ## Convergence criteria

# The argumento `maxiter` specifies the maximum number of iterations that should be carried
# out. At each iteration,

# ``aₖⁿ⁺¹ = aₖⁿ + δₖ``

# Convergence is achieved when

# ``|δᵏ / aₖ⁰| < ε``

# ## Example
# ```julia
# x = 1.0:10.0
# a = [3.0, 2.0, 1.0]
# y = a[1] + a[2]*x + a[3]*x^2
# fun(x, a) = a[1] + a[2]*x + a[3]*x^2

# a = secant_nls_fit(x, y, fun, [0.5, 0.5, 0.5], 1e-8, 30)
# ```
# """
# function secant_nls_fit(
#         x::AbstractVector{T}, y::AbstractVector{T}, fun, aguess::AbstractVector{T},
#         eps = 1e-8, maxiter = 200) where {T <: Number}
#     # TODO: once NonlinearSolve.jl has derivative-free NLLS methods, we can remove this
#     P = length(x) # Number of points
#     N = length(aguess) # Number of parameters

#     xi = zero(T)
#     df = zeros(T, N)
#     a = zeros(T, N)
#     for i in 1:N
#         a[i] = aguess[i]
#         if a[i] == 0
#             a[i] = 0.01
#         end
#     end

#     δ = a .* (one(T) / 20)
#     f1 = zeros(T, P)
#     a .+= δ
#     A = zeros(T, N, N)
#     b = zeros(T, N)

#     δref = abs.(a)
#     maxerr = zero(T)
#     for iter in 1:maxiter
#         A .= zero(T)
#         b .= zero(T)
#         f1 .= fun.(x, Ref(a))
#         for i in 1:P
#             xi = x[i]
#             yi = y[i]
#             f = f1[i] - yi
#             for k in 1:N
#                 a[k] -= δ[k]
#                 df[k] = (f1[i] - fun(xi, a)) / δ[k]
#                 a[k] += δ[k]
#             end
#             # Assemble LHS
#             for k in 1:N
#                 for j in 1:N
#                     A[j, k] += df[j] * df[k]
#                 end
#             end
#             # Assemble RHS
#             for j in 1:N
#                 b[j] -= f * df[j]
#             end
#         end
#         δ = A \ b
#         a .+= δ
#         # Verify convergence:
#         maxerr = maximum(abs, δ ./ δref)
#         if maxerr < eps
#             return (a)
#         end
#     end

#     error("gauss_newton_fit failed to converge in $maxiter iterations with relative residpal of $maxerr !")

#     return (a)
# end
