# using LinearAlgebra: diagm, qr!, eigvals

# struct ExpSumFit{T <: Real, S <: Number} <: AbstractLeastSquares
#     k::T          # constant
#     p::Vector{S}  # coefficients
#     λ::Vector{S}  # rate
#     τ::Vector{S}  # -1/rate
# end

# struct ExpSumFitInit{T <: Real}
#     n::Int                 # number of exponentials in the sum
#     m::Int                 # split integration intervals
#     Y::Matrix{T}           #
#     S::Matrix{T}           #
#     A::Vector{T}           #
#     Ā::Matrix{T}           # companion matrix
#     Xc::Matrix{Complex{T}} #
#     Xr::Matrix{T}          #
#     coeff::Matrix{T}       # numerical integration coefficients
# end

# """
#         expsum_fit(
#             x::Union{T,AbstractVector{T}}, y::AbstractVector{T}, n::Int; m::Int = 1,
#             withconst::Bool = true, init::Union{Nothing,ExpSumFitInit} = nothing
#         ) where {T <: Real}

# Fits the sum of `n` exponentials and a constant.

# ```math
#     y = k + p_1 e^{\\lambda_1 t} + p_2 e^{\\lambda_2 t} + \\ldots + p_n e^{\\lambda_n t}
# ```

# If the keyword `withconst` is set to `false`, the constant is not fitted but set `k=0`.

# Uses numerical integration with `m` strips, where the default `m=1` uses linear
# interpolation. `m=2` and higher require uniform interval and usually lead to better
# accuracy.

# Passing `init` preallocates most needed memory, and can be initialized with
# [`expsum_init`](@ref).

# Returns a `ExpSumFit` struct containing a constant `k` and vectors `p`, `λ`, `τ`, where
# `τ = -1/λ`.

# The algorithm is from
# [Matlab code of Juan Gonzales Burgos](https://github.com/juangburgos/FitSumExponentials).
# """
# function expsum_fit(x::AbstractVector{T}, y::AbstractVector{T},
#         n::Int; m::Int = 1, withconst::Bool = true,
#         init::Union{Nothing, ExpSumFitInit} = nothing) where {T <: Real}
#     n > 0 || throw(ArgumentError("number of exponent terms should be a positive integer"))
#     length(x) == length(y) ||
#         throw(DimensionMismatch("input vectors should be equal in length"))
#     m ≥ 2 && !all(isapprox.(diff(x), x[2] - x[1])) &&
#         error("m=$m requires uniformly spaced x")
#     if isnothing(init)
#         init = expsum_init(x, n, m = m, withconst = withconst)
#     else
#         (init.n == n && init.m == m && size(init.Xc)[1] == length(x)) ||
#             error("Init does not match parameters")
#     end
#     sc = expsum_scale!(x, y)
#     expfit_solve(init, x, y, sc)
# end

# function expsum_fit(dx::T, y::AbstractVector{T}, n::Int; kwargs...) where {T <: Real}
#     x = collect(range(0, step = dx, length = length(y)))
#     expsum_fit(x, y, n; kwargs...)
# end

# function expsum_scale!(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
#     sc = (x = maximum(abs, x), y = maximum(abs, y))
#     x ./= sc.x
#     y ./= sc.y
#     return sc
# end

# """
#     expsum_init(x::AbstractVector{T}, n::Int; m::Int = 1, withconst::Bool = true) where {T <: Real}

# Initialize most of the memory needed for [`expsum_fit`](@ref).
# """
# function expsum_init(
#         x::AbstractVector{T}, n::Int; m::Int = 1, withconst::Bool = true) where {T <: Real}
#     len = length(x)
#     nY, mY = 1 + (len - 1) ÷ m, 2n + withconst
#     return ExpSumFitInit(n, m, zeros(T, nY, mY), zeros(T, nY, n - 1),
#         zeros(T, mY),
#         diagm(-1 => ones(T, n - 1)), # companion matrix
#         zeros(Complex{T}, len, n + withconst),
#         zeros(T, len, n + withconst),
#         T.(calc_integral_rules(1:n, m = m)))
# end

# function expfit_solve(init::ExpSumFitInit, x, y, sc::NamedTuple)
#     n, Y, n = init.n, init.Y, init.n
#     expsum_fill_Y!(init, x, y)
#     λ, τ = expfit_solve_λ(init, y)
#     if isreal(λ)
#         X = init.Xr
#     else
#         X = init.Xc
#     end
#     expsum_fill_X!(init, x, λ, X)
#     qrX = qr!(X)
#     p = qrX \ y
#     if isreal(p)
#         p = real(p)
#     end
#     x .*= sc.x
#     y .*= sc.y
#     withconst = size(Y)[2] == 2n + 1
#     if withconst
#         return ExpSumFit(real(p[end]) * sc.y, p[1:n] * sc.y, λ / sc.x, τ * sc.x)
#     else
#         return ExpSumFit(0.0, p * sc.y, λ / sc.x, τ * sc.x)
#     end
# end

# """
#     cumints!(init::ExpSumFitInit, x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}

# Calculates in-place `init.n` cumulative integrals using coeffients `init.coeff`.
# """
# function cumints!(
#         init::ExpSumFitInit, x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
#     n, m, Y, S, coeff = init.n, init.m, init.Y, init.S, init.coeff
#     n > 0 || error("Number of exponential terms should be a positive integer")
#     m > 0 || error("Order of interpolating polynomial should be a positive integer")
#     len = size(init.Y)[1]
#     S .= 0.0
#     for j in 1:n
#         Y[1, j] = 0.0
#         for i in 2:len
#             dx = x[m * (i - 2) + 2] - x[m * (i - 2) + 1]
#             s = 0.0
#             for k in 1:(m + 1)
#                 s += coeff[j, k] * y[m * (i - 2) + k]
#             end
#             Y[i, j] = Y[i - 1, j] + dx^j * s
#         end
#     end
#     for j in 2:n
#         for i in 2:len
#             S[i, j - 1] += S[i - 1, j - 1] + Y[i, j - 1] # at this stage Y[:,j-1] is calculated
#         end
#         for k in 1:(j - 1)
#             f = factorial(k)
#             for i in 2:len
#                 dx = x[m * (i - 2) + 2] - x[m * (i - 2) + 1]
#                 Y[i, j] += S[i - 1, j - k] / f * (m * dx)^k
#             end
#         end
#     end
#     nothing
# end

# """
#     calc_integral_rules(n; m = 2)

# Determine coefficients of the rules cumulative integrals of order `n` using
# [method of undetermined coefficients](https://en.wikipedia.org/wiki/Simpson%27s_rule#Undetermined_coefficients).
# Interpolation order is `m`.

# * `n=1`, `m=1` [Trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule)
# * `n=1`, `m=2` [Simpson's 1/3 rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_1/3_rule)
# * `n=1`, `m=3` [Simpson's second (3/8) rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule)
# """
# function calc_integral_rules(n::Int; m::Int = 2)
#     n ≥ 1 || throw(ArgumentError("n=$n should be positive integer"))
#     n = big(n)
#     # evaluate m-th order polynomial terms at points x=0:m
#     polyvals = [x^n for x in 0:Rational(m), n in 0:m]
#     # evaluate m-th order polynomial terms integrated cumulatively n-times
#     integralvals = [m^(n + i) * factorial(i) // factorial(n + i) for i in 0:m]'
#     coeff = integralvals / polyvals
#     return coeff
# end

# function calc_integral_rules(ns::Union{UnitRange{Int}, Vector{Int}}; m::Int = 2)
#     reduce(vcat, [calc_integral_rules(n, m = m) for n in ns])
# end

# """
#     expsum_fill_Y!(init::ExpSumFitInit, x, y)

# Fills matrix `Y` with
# ```math
#     [\\int^1 \\int^2 \\ldots \\int^n x^0 x^1 \\ldots x^{m-n-1}]
# ```
# where `∫ⁱ` means the ith cumulative integral of `y` vs `x`, and `m` is the
# number of columns in matrix `Y`.
# """
# function expsum_fill_Y!(
#         init::ExpSumFitInit, x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
#     n, m, Y = init.n, init.m, init.Y
#     nY, mY = size(Y)
#     cumints!(init, x, y)
#     for j in 1:mY
#         Y[1, j] = 0.0
#     end
#     for j in 0:(mY - n - 1)
#         for i in 1:nY
#             Y[i, j + n + 1] = x[m * (i - 1) + 1]^j
#         end
#     end
#     nothing
# end

# """
#     expsum_fill_X!(init::ExpSumFitInit, x, λ, X)

# Fills matrix `init.X`
# ```math
#     [e^{x\\lambda_1} e^{x\\lambda_2} \\ldots e^{x\\lambda_n}]
# ```
# or if fitting is done with constant
# ```math
#     [e^{x\\lambda_1} e^{x\\lambda_2} \\ldots e^{x\\lambda_n} 1]
# ```
# """
# function expsum_fill_X!(init::ExpSumFitInit, x::AbstractVector{T}, λ::AbstractVector{S},
#         X::AbstractMatrix{S}) where {T <: Real, S <: Number}
#     n = init.n
#     nX, mX = size(X)
#     for j in 1:n
#         for i in 1:nX
#             X[i, j] = exp(x[i] * λ[j])
#         end
#     end
#     for j in (n + 1):mX
#         for i in 1:nX
#             X[i, j] = 1.0
#         end
#     end
#     nothing
# end

# function expfit_solve_λ(init::ExpSumFitInit, y)
#     n, m, Y, A, Ā = init.n, init.m, init.Y, init.A, init.Ā
#     qrY = qr!(Y)
#     A .= qrY \ view(y, 1:m:length(y))
#     Ā[1, 1:n] = A[1:n]
#     λ = eigvals(Ā)
#     if isreal(λ)
#         λ = real(λ)
#     end
#     τ = -1 ./ λ
#     return λ, τ
# end

# """
#     (sol::ExpSumFit)(x)

# Calculate the sum of exponentials using solution `sol` at points `x`.
# """
# function (f::ExpSumFit)(x)
#     k, p, λ = f.k, f.p, f.λ
#     y = k .+ sum(exp.(x * λ') .* p', dims = 2)
#     y = real.(y[:])
#     return y
# end
