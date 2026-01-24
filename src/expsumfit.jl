"""
    __calc_integral_rules(::Type{T}, n; m = 2) where {T}

Determine coefficients of the rules cumulative integrals of order `n` using
[method of undetermined coefficients](https://en.wikipedia.org/wiki/Simpson%27s_rule#Undetermined_coefficients).
Interpolation order is `m`.

* `n=1`, `m=1` [Trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule)
* `n=1`, `m=2` [Simpson's 1/3 rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_1/3_rule)
* `n=1`, `m=3` [Simpson's second (3/8) rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule)
"""
function __calc_integral_rules(::Type{T}, ns::AbstractVector{Int}; m::Int = 2) where {T}
    # evaluate m-th order polynomial terms at points x = 0:m
    polyvals = Matrix{FastRational{Int}}(undef, m + 1, m + 1)
    @inbounds for i in 0:m
        @simd ivdep for j in 0:m
            polyvals[j + 1, i + 1] = FastRational{Int}(i, 1, Val(true))^j
        end
    end
    polyvals_factorized = lu!(polyvals)

    integralvals = Matrix{FastRational{Int128}}(undef, length(ns), m + 1)
    result = similar(integralvals, T)
    @inbounds for (i, n) in enumerate(ns)
        __calc_integral_rules!(polyvals_factorized, view(integralvals, i, :), Int128(n), m)
        @simd ivdep for j in 1:(m + 1)
            result[i, j] = T(integralvals[i, j])
        end
    end
    return result
end

function __calc_integral_rules!(
        polyvals_factorized, integralvals::AbstractVector, n::nT, m::Integer
    ) where {nT <: Integer}
    @assert n ≥ 1 "n=$n should be positive integer"
    @assert n + m ≤ 33 "m + n=$n + $m should be less than or equal to 33"

    # evaluate m-th order polynomial terms integrated cumulatively n-times
    num = m^(n - 1)
    @inbounds @simd ivdep for i in 0:m
        num *= m
        integralvals[i + 1] = num * FastRational{nT}(
            factorial(nT(i)), factorial(n + i), Val(true)
        )
    end

    ldiv!(polyvals_factorized, integralvals)
    return nothing
end

function __expsum_scale!(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Real}
    scx = maximum(abs, x)
    scy = maximum(abs, y)
    @inbounds @simd ivdep for i in eachindex(x, y)
        x[i] /= scx
        y[i] /= scy
    end
    return (; x = scx, y = scy)
end

function __cumulative_integrals!(
        Y::AbstractMatrix{T}, S::AbstractMatrix{T}, coeff::AbstractMatrix{T},
        x::AbstractVector{T}, y::AbstractVector{T}, n::Integer, m::Integer
    ) where {T <: Real}
    @assert n > 0 "n=$n should be a positive integer"
    @assert m > 0 "m=$m should be a positive integer"

    fill!(S, false)

    @inbounds for j in 1:n
        Y[1, j] = false
        for i in 2:size(Y, 1)
            dx = x[m * (i - 2) + 2] - x[m * (i - 2) + 1]
            s = zero(T)
            @simd ivdep for k in 1:(m + 1)
                s += coeff[j, k] * y[m * (i - 2) + k]
            end
            Y[i, j] = Y[i - 1, j] + dx^j * s
        end
    end

    @inbounds for j in 2:n
        for i in 2:size(Y, 1)
            S[i, j - 1] += S[i - 1, j - 1] + Y[i, j - 1] # at this stage Y[:,j-1] is calculated
        end
        for k in 1:(j - 1)
            f = factorial(k)
            @simd ivdep for i in 2:size(Y, 1)
                dx = x[m * (i - 2) + 2] - x[m * (i - 2) + 1]
                Y[i, j] += S[i - 1, j - k] / f * (m * dx)^k
            end
        end
    end

    return nothing
end

function __expsum_fill_X!(
        x::AbstractVector, λ::AbstractVector, X::AbstractMatrix, n::Integer
    )
    @inbounds for j in 1:n
        @simd ivdep for i in axes(X, 1)
            X[i, j] = exp(x[i] * λ[j])
        end
    end
    fill!(view(X, :, (n + 1):size(X, 2)), true)
    return nothing
end

function __expsum_fill_Y!(
        Y::AbstractMatrix{T}, S::AbstractMatrix{T}, coeff::AbstractMatrix{T},
        x::AbstractVector{T}, y::AbstractVector{T}, n::Integer, m::Integer
    ) where {T <: Real}
    __cumulative_integrals!(Y, S, coeff, x, y, n, m)

    fill!(view(Y, 1, :), false)

    @inbounds for j in 0:(size(Y, 2) - n - 1)
        @simd ivdep for i in axes(Y, 1)
            Y[i, j + n + 1] = x[m * (i - 1) + 1]^j
        end
    end
    return nothing
end

function __expsum_solve_λ(Y, A, Ā, y, n, m)
    qrY = qr!(Y)
    ldiv!(A, qrY, view(y, 1:m:length(y)))
    copyto!(view(Ā, 1, 1:n), view(A, 1:n))
    return eigvals!(Ā)
end

# Common Solve Interface

@concrete struct ExpSumFitCache <: AbstractCurveFitCache
    prob <: CurveFitProblem
    alg
    kwargs
    Y <: AbstractMatrix
    S <: AbstractMatrix
    A <: AbstractVector
    Ā <: AbstractMatrix
    Xc <: AbstractMatrix{<:Complex}
    Xr <: AbstractMatrix
    coeff <: AbstractMatrix
end

function CommonSolve.init(
        prob::CurveFitProblem, alg::ExpSumFitAlgorithm; kwargs...
    )
    @assert !is_nonlinear_problem(prob) "Exponential sum fitting only works with linear \
                                         problems"
    sigma_not_supported(prob)

    T = eltype(prob.x)

    len = length(prob.x)
    nY, mY = 1 + (len - 1) ÷ alg.m, 2 * alg.n + alg.withconst
    return ExpSumFitCache(
        prob,
        alg,
        kwargs,
        similar(prob.x, nY, mY),
        similar(prob.x, nY, alg.n - 1),
        similar(prob.x, mY),
        diagm(-1 => ones(T, alg.n - 1)),
        similar(prob.x, Complex{T}, len, alg.n + alg.withconst),
        similar(prob.x, len, alg.n + alg.withconst),
        __calc_integral_rules(T, 1:(alg.n); alg.m)
    )
end

# TODO: allocations in this function aren't optimized
function CommonSolve.solve!(cache::ExpSumFitCache)
    sc = __expsum_scale!(cache.prob.x, cache.prob.y)

    __expsum_fill_Y!(
        cache.Y, cache.S, cache.coeff, cache.prob.x, cache.prob.y, cache.alg.n, cache.alg.m
    )

    λ = __expsum_solve_λ(
        cache.Y, cache.A, cache.Ā, cache.prob.y, cache.alg.n, cache.alg.m
    )

    X = isreal(λ) ? cache.Xr : cache.Xc
    __expsum_fill_X!(cache.prob.x, λ, X, cache.alg.n)

    qrX = qr!(X)
    p = qrX \ cache.prob.y
    isreal(p) && (p = real(p))

    for i in eachindex(cache.prob.x, cache.prob.y)
        cache.prob.x[i] *= sc.x
        cache.prob.y[i] *= sc.y
    end

    for i in eachindex(p)
        p[i] *= sc.y
    end
    for i in eachindex(λ)
        λ[i] /= sc.x
    end

    withconst = size(cache.Y, 2) == 2 * cache.alg.n + 1
    T = promote_type(eltype(p), eltype(λ))

    if withconst
        k = T[real(p[end])]
        p = p[1:(cache.alg.n)]
    else
        k = T[0]
    end
    backing = (; k, p, λ)

    y_pred = k .+ sum(exp.(cache.prob.x * λ') .* p'; dims = 2)
    resid = cache.prob.y .- y_pred

    return CurveFitSolution(
        cache.alg, NamedArrayPartition(backing), resid, cache.prob, ReturnCode.Success, nothing
    )
end

function (sol::CurveFitSolution{<:ExpSumFitAlgorithm})(x)
    (; k, p, λ) = sol.u
    y = k .+ sum(exp.(x * λ') .* p'; dims = 2)
    return real.(vec(y))
end
