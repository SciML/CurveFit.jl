
# """
# Original Kings law (1910) represents the relation between voltage and velocity
# in a hotwire anemometer. The law is given by:

# `E^2 = A + B * U^0.5`

# This function estimates `A` and `B`.
# """
# linear_king_fit(E, U) = __linear_fit_internal(sqrt, U, abs2, E)

# """
# Type that represents a Linear (original) King's law
# """
# struct LinearKingFit{T <: Number} <: AbstractLeastSquares
#     coefs::NTuple{2, T}

#     LinearKingFit(A::T, B::T) where {T <: Number} = new{T}((A, B))
#     LinearKingFit(coefs::NTuple{2, T}) where {T <: Number} = new{T}(coefs)
#     function LinearKingFit(E::AbstractVector{T}, U::AbstractVector{T}) where {T <: Number}
#         return new{T}(linear_king_fit(E, U))
#     end
# end

# (f::LinearKingFit)(E) = ((E .* E .- f.coefs[1]) ./ f.coefs[2]) .^ 2

# """
# Equation that computes the error of the modified King's law
# """
# function kingfun!(y, x, a)
#     @simd ivdep for i in eachindex(y)
#         y[i] = a[1] + a[2] * x[2, i]^a[3] - x[1, i]^2
#     end
#     return y
# end

# """
# Uses nonlinear least squares to fit the modified King's law:

# `E^2 = A + B * U^n`

# The Original (linear) King's law is used to estimate `A` and `B` when `n = 1/2`.
# This initial value is used as an initial guess for fitting the nonlinear modified King's law
# using the function `nonlinear_fit`.
# """
# function king_fit(E, U, args...; kwargs...)
#     a, b = linear_king_fit(E, U)
#     T = promote_type(typeof(a), typeof(b))

#     sol = nonlinear_fit(
#         kingfun!, stack((E, U); dims=1), [T(a), T(b), T(0.5)], args...;
#         resid_prototype = Vector{T}(undef, length(E)), iip=Val(true), kwargs...
#     )
#     return (sol.u[1], sol.u[2], sol.u[3])
# end

# """
# Type that represents the modified King's law
# """
# struct KingFit{T <: Number} <: AbstractLeastSquares
#     coefs::NTuple{3, T}

#     KingFit(A::T, B::T, n::T) where {T <: Number} = new{T}((A, B, n))
#     KingFit(coefs::NTuple{3, T}) where {T <: Number} = new{T}(coefs)
#     KingFit(A::T, B::T) where {T <: Number} = new{T}((A, B, one(T) / 2))
#     KingFit(coefs::NTuple{2, T}) where {T <: Number} = new{T}(coefs[1], coefs[2])
# end

# KingFit(args...; kwargs...) = KingFit(king_fit(args...; kwargs...))

# (f::KingFit)(E) = ((E * E - f.coefs[1]) / f.coefs[2])^(1 / f.coefs[3])
