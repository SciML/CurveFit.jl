abstract type AbstractCurveFitProblem end

abstract type AbstractCurveFitAlgorithm end

abstract type AbstractCurveFitSolution end

abstract type AbstractCurveFitCache end

# TODO: print more information about the cache
Base.show(io::IO, ::MIME"text/plain", x::T) where {T <: AbstractCurveFitCache} = print(io, nameof(T), "()")

# Core Problem Types
"""
    CurveFitProblem(x, y; nlfunc=nothing, u0=nothing, sigma=nothing, lb=nothing, ub=nothing)

Represents a curve fitting problem where `x` and `y` are the data points to fit.

Certain algorithms may require an initial guess `u0` for the coefficients to fit. See
specific solver documentation for more details.

Weights can be passed through `sigma`, which should be an array with the same
dimensions as `y`. As with `curve_fit()` from scipy, the elements should be the
standard deviation of `y`. Note that currently `sigma` is not supported for all
kinds of fits, check the problem or algorithm docstring to see if sigma is
supported.

Lower and upper bounds on the parameters can be passed through `lb` and `ub`,
which should be arrays with the same length as `u0`. Note that currently
`bounds` is not supported for all kinds of fits, check the problem or algorithm
docstring to see if bounds are supported.

See also [`NonlinearCurveFitProblem`](@ref).
"""
@concrete struct CurveFitProblem <: AbstractCurveFitProblem
    x <: AbstractArray
    y <: Union{AbstractArray, Nothing}
    sigma <: Union{AbstractArray, Nothing}
    nlfunc <: Union{Nothing, NonlinearFunction}
    u0 <: Union{Nothing, AbstractArray}
    lb <: Union{Nothing, AbstractArray}
    ub <: Union{Nothing, AbstractArray}
end

function SciMLBase.isinplace(prob::CurveFitProblem)
    !is_nonlinear_problem(prob) && return false
    return SciMLBase.isinplace(prob.nlfunc)
end

is_nonlinear_problem(prob::CurveFitProblem) = prob.nlfunc !== nothing

function sigma_not_supported(prob::CurveFitProblem)
    @assert isnothing(prob.sigma) "Passing weights (sigma) is not supported for this algorithm"
    return
end

function bounds_not_supported(prob::CurveFitProblem)
    @assert isnothing(prob.lb) && isnothing(prob.ub) "Passing bounds (lb/ub) is not supported for this algorithm"
    return
end

function CurveFitProblem(x, y; nlfunc = nothing, u0 = nothing, sigma = nothing, lb = nothing, ub = nothing)
    if nlfunc === nothing
        @assert ndims(x) == ndims(y) == 1 "x and y must be 1-dimensional arrays for linear \
                                       problems (`nlfunc` is `nothing`)"
    end

    return CurveFitProblem(x, y, sigma, nlfunc, u0, lb, ub)
end

@doc doc"""
    NonlinearCurveFitProblem(f, u0, x, y=nothing, sigma=nothing; lb=nothing, ub=nothing)

Nonlinear curve fitting problem where `f` is a nonlinear function to fit, `u0` is the
initial guess for the coefficients, `x` and `y` are the data points to fit, and
`sigma` is the standard deviation associated with `y`. The following
optimization problem is solved:

```math
\argmin_u ~ \left\| f(u, x) - y \right\|_2
```

If `y` is `nothing`, then it is treated as a zero vector and `f` is expected to
return residuals directly. `f` is a generic Julia function or ideally a
`NonlinearFunction` from
[`SciMLBase.jl`](https://github.com/SciML/SciMLBase.jl).

Lower and upper bounds on the parameters can be passed through `lb` and `ub` keyword
arguments. These should be arrays with the same length as `u0`.

## Function Signature

The model function `f` should have the signature `f(params, x)` where:
- `params` is a vector of parameters to be fitted
- `x` is the input data (can be a vector or matrix)

The function should return predictions for all input data points. For vectorized operations
over arrays, use Julia's broadcasting syntax with the `@.` macro:

```julia
# Vectorized function using @.
fn(a, x) = @. a[1] + a[2] * x^a[3]
```

For users who prefer to define scalar functions (e.g., those migrating from LsqFit.jl),
use the [`ScalarModel`](@ref) wrapper:

```julia
# Scalar function (operates on single x value)
fn_scalar(a, x) = a[1] + a[2] * x^a[3]
prob = NonlinearCurveFitProblem(ScalarModel(fn_scalar), u0, x, y)
```

See also [`ScalarModel`](@ref).
"""
function NonlinearCurveFitProblem(f::NonlinearFunction, u0, x, y = nothing, sigma = nothing; lb = nothing, ub = nothing)
    return CurveFitProblem(x, y; nlfunc = f, u0 = u0, sigma, lb, ub)
end
function NonlinearCurveFitProblem(f::F, u0, x, y = nothing, sigma = nothing; lb = nothing, ub = nothing) where {F}
    return NonlinearCurveFitProblem(NonlinearFunction(f), u0, x, y, sigma; lb, ub)
end

"""
    ScalarModel(f)

Wraps a scalar function `f(params, x_i)` that operates on a single data point `x_i`
into a vectorized form suitable for CurveFit.jl.

This is useful for users migrating from LsqFit.jl or those who prefer defining
scalar model functions without explicit broadcasting via the `@.` macro.

## Why is `@.` Typically Required?

In CurveFit.jl, model functions receive the entire data array `x` at once and must
return predictions for all data points. This vectorized design enables:
- Better GPU performance when using GPU arrays
- More efficient compilation with tools like Reactant.jl
- User control over array types

## Using ScalarModel

Instead of writing a vectorized function:
```julia
# Vectorized function (uses @.)
fn(a, x) = @. a[1] + a[2] * x^a[3]
prob = NonlinearCurveFitProblem(fn, u0, x, y)
```

You can write a simpler scalar function:
```julia
# Scalar function (no @. needed)
fn_scalar(a, x) = a[1] + a[2] * x^a[3]
prob = NonlinearCurveFitProblem(ScalarModel(fn_scalar), u0, x, y)
```

## Migration from LsqFit.jl

For users coming from LsqFit.jl, note that the parameter order is reversed:
- LsqFit.jl: `model(x, p)` (data first, then parameters)
- CurveFit.jl: `model(p, x)` (parameters first, then data)

Example migration:
```julia
# LsqFit.jl style (does NOT work directly)
lsqfit_model(x, p) = p[1] * exp(-x * p[2])

# CurveFit.jl with ScalarModel
curvefit_model(p, x) = p[1] * exp(-x * p[2])
prob = NonlinearCurveFitProblem(ScalarModel(curvefit_model), u0, x, y)
```
"""
struct ScalarModel{F}
    f::F
end

# When called with array data, broadcast over the data
(sm::ScalarModel)(out, params, x::AbstractArray) = out .= sm.f.(Ref(params), x)
(sm::ScalarModel)(params, x::AbstractArray) = sm.f.(Ref(params), x)
# When called with scalar data (for single-point evaluation), call directly
(sm::ScalarModel)(params, x::Number) = sm.f(params, x)

# Algorithms
@concrete struct LinearCurveFitAlgorithm <: AbstractCurveFitAlgorithm
    xfun <: Function
    yfun <: Function
    yfun_inverse <: Function
end

"""
    LinearCurveFitAlgorithm(;
        xfun = identity, yfun = identity, yfun_inverse = inverse(yfun)
    )

Represents a linear curve fitting algorithm where `x` and `y` are the data
points to fit. If the [`CurveFitProblem`](@ref) being solved has a `sigma` then
it will be used as weights.
We want to solve for `a` and `b` such that:

```math
f_y(y) = a f_x(x) + b
```

where ``f_x`` corresponds to `xfun` and ``f_y`` corresponds to `yfun`.
Note that this is a general problem specification of a curve fitting problem which can
be converted to a linear fit in a specific function space by choosing appropriate
`xfun` and `yfun`. The `yfun_inverse` is used to convert the fitted values back to the
original space (can be specified by defining `InverseFunctions.inverse`).

When `yfun` is not `identity`, `yfun_inverse` is applied to the fitted intercept before
storing it, so that `sol.u = (a, b)` directly satisfies the original-space formula
(e.g. `y = b * exp(a*x)` for [`ExpCurveFitAlgorithm`](@ref)).

This algorithm does not support bounds constraints (`lb`/`ub`).
"""
function LinearCurveFitAlgorithm(;
        xfun = identity, yfun = identity, yfun_inverse = inverse(yfun)
    )
    return LinearCurveFitAlgorithm(xfun, yfun, yfun_inverse)
end

@doc doc"""
    LogCurveFitAlgorithm()

Represents a log curve fitting algorithm where `x` and `y` are the data points
to fit. If the [`CurveFitProblem`](@ref) being solved has a `sigma` then
it will be used as weights. This algorithm does not support bounds constraints (`lb`/`ub`).
We want to solve for `a` and `b` such that:

```math
y = a \log(x) + b
```
"""
LogCurveFitAlgorithm() = LinearCurveFitAlgorithm(; xfun = log, yfun = identity)

@doc doc"""
    PowerCurveFitAlgorithm()

Represents a power curve fitting algorithm where `x` and `y` are the data points
to fit. This algorithm does not support passing weights through `sigma` in
[`CurveFitProblem`](@ref). This algorithm does not support bounds constraints (`lb`/`ub`).
We want to solve for `a` and `b` such that:

```math
y = b x^a
```

This is equivalent to a linear fit in log-log space, i.e.,

```math
\log(y) = a \log(x) + \log(b)
```

Because the fit is performed in log space, it minimizes *relative* (rather than
absolute) error: a fixed fractional deviation counts the same at small and large
`y`. This is a different estimator from a nonlinear least-squares fit of the same
model (which minimizes absolute error), so the coefficients and their
uncertainties will generally differ. [`vcov`](@ref) accounts for this when
computing the uncertainty of `(a, b)`.
"""
PowerCurveFitAlgorithm() = LinearCurveFitAlgorithm(; xfun = log, yfun = log)

@doc doc"""
    ExpCurveFitAlgorithm()

Represents an exponential curve fitting algorithm where `x` and `y` are the data points to
fit. This algorithm does not support passing weights through `sigma` in
[`CurveFitProblem`](@ref). This algorithm does not support bounds constraints (`lb`/`ub`).
We want to solve for `a` and `b` such that:

```math
y = b \exp(a x)
```

This is equivalent to a linear fit in log-linear space, i.e.,

```math
\log(y) = a x + \log(b)
```

Because the fit is performed in log space, it minimizes *relative* (rather than
absolute) error: a fixed fractional deviation counts the same at small and large
`y`. This is a different estimator from a nonlinear least-squares fit of the same
model (which minimizes absolute error), so the coefficients and their
uncertainties will generally differ. [`vcov`](@ref) accounts for this when
computing the uncertainty of `(a, b)`.
"""
ExpCurveFitAlgorithm() = LinearCurveFitAlgorithm(; xfun = identity, yfun = log)

@doc doc"""
    KingCurveFitAlgorithm()

Represents a king curve fitting problem where `x` and `y` are the data points to
fit. This algorithm does not support passing weights through `sigma` in
[`CurveFitProblem`](@ref). This algorithm does not support bounds constraints (`lb`/`ub`).
We want to solve for `A` and `B` according to original King's law (1910) that represents
the relationship between voltage (E) and velocity (U) in a hotwire anemometer:

```math
E^2 = A + B U^{1/2}
```

or

```math
x^2 = A + B y^{1/2}
```
"""
struct KingCurveFitAlgorithm <: AbstractCurveFitAlgorithm end

@doc doc"""
    ModifiedKingCurveFitAlgorithm(alg::Union{Nothing, AbstractNonlinearAlgorithm} = nothing)

Similar to [`KingCurveFitAlgorithm`](@ref), but uses the modified King's law:

```math
E^2 = A + B U^n
```

where `n` is also a parameter.

This algorithm supports bounds constraints via `lb` and `ub` in
[`CurveFitProblem`](@ref).
"""
@kwdef @concrete struct ModifiedKingCurveFitAlgorithm <: AbstractCurveFitAlgorithm
    alg <: Union{Nothing, AbstractNonlinearAlgorithm} = nothing
end

"""
    PolynomialFitAlgorithm(degree::Int)
    PolynomialFitAlgorithm(;
        degree::Int,
        linsolve_algorithm::Union{Nothing, AbstractLinearAlgorithm} = nothing
    )

Represents a polynomial fitting algorithm of degree `degree`. Only applicable to
[`LinearCurveFitAlgorithm`](@ref)s. This algorithm does not support passing
weights through `sigma` in [`CurveFitProblem`](@ref). This algorithm does not support
bounds constraints (`lb`/`ub`).

!!! tip

    For ill-conditioned problems, it is recommended to use linear solvers like
    `QRFactorization`. Alternatively, pass in
    `assumptions = OperatorAssumptions(false; condition = OperatorsCondition.<condition>)`
    to `solve`/`init`.
"""
@kwdef @concrete struct PolynomialFitAlgorithm <: AbstractCurveFitAlgorithm
    degree::Int
    linsolve_algorithm <: Union{Nothing, AbstractLinearAlgorithm} = nothing
end

PolynomialFitAlgorithm(degree::Int) = PolynomialFitAlgorithm(degree, nothing)

"""
    RationalPolynomialFitAlgorithm(num_degree::Int, den_degree::Int)
    RationalPolynomialFitAlgorithm(;
        num_degree::Int, den_degree::Int, alg = nothing
    )

Represents a rational polynomial fitting algorithm with numerator degree `num_degree`
and denominator degree `den_degree`. The internal polynomial fitting algorithm is
determined by the `alg` keyword argument. If `alg` is `nothing` or a
`AbstractNonlinearAlgorithm` (like solvers from NonlinearSolve.jl), it will use a
nonlinear curve fitting approach. If `alg` is a `AbstractLinearAlgorithm`, it will use
linear least squares fitting. This algorithm does not support passing weights
through `sigma` in [`CurveFitProblem`](@ref).

## Linear Rational Polynomial Fitting

In this case the following curve fit is done:

```math
y = \\frac{p(x)}{q(x)}
```

where `p(x)` is a polynomial of degree `num_degree` and `q(x)` is a polynomial of degree
`den_degree`. The linear case is solved by doing a least squares fit on:

```math
y q(x) = p(x)
```

where the zero order term of `q(x)` is assumed to be 1.

## Nonlinear Rational Polynomial Fitting

If an `u0` is not provided to the problem, then we will use linear least squares for an
initial guess.

The nonlinear variant of this algorithm supports bounds constraints via `lb` and `ub` in
[`CurveFitProblem`](@ref).
"""
@kwdef @concrete struct RationalPolynomialFitAlgorithm <: AbstractCurveFitAlgorithm
    num_degree::Int
    den_degree::Int
    alg <: Union{Nothing, AbstractLinearAlgorithm, AbstractNonlinearAlgorithm} = nothing
end

function RationalPolynomialFitAlgorithm(num_degree::Int, den_degree::Int)
    return RationalPolynomialFitAlgorithm(num_degree, den_degree, nothing)
end

@doc doc"""
    ExpSumFitAlgorithm(; n::Int, m::Int = 1, withconst::Bool = true)

Fits the sum of `n` exponentials and a constant. This algorithm does not support
passing weights through `sigma` in [`CurveFitProblem`](@ref). This algorithm does not
support bounds constraints (`lb`/`ub`).

```math
y = k + p_1 e^{λ_1 t} + p_2 e^{λ_2 t} + ⋯ + p_n e^{λ_n t}
```

If the keyword `withconst` is set to `false`, the constant is not fitted but set `k=0`.

Uses numerical integration with `m` strips, where the default `m=1` uses linear
interpolation. `m=2` and higher require uniform interval and usually lead to better
accuracy.

This algorithm is from
[Matlab code of Juan Gonzales Burgos](https://github.com/juangburgos/FitSumExponentials).
"""
@kwdef @concrete struct ExpSumFitAlgorithm <: AbstractCurveFitAlgorithm
    n::Int
    m::Int = 1
    withconst::Bool = true
end

## Internal types for dispatch
struct __FallbackLinearFitAlgorithm <: AbstractCurveFitAlgorithm end
@concrete struct __FallbackNonlinearFitAlgorithm <: AbstractCurveFitAlgorithm
    alg <: Union{Nothing, AbstractNonlinearAlgorithm}
end

# Solution Types
"""
    CurveFitSolution(alg, coeffs, resid, prob, retcode, original=nothing)

Represents the solution to a curve fitting problem. This is a callable struct and
can be used to evaluate the solution at a point. Exact evaluation mechanism depends on the
algorithm used to solve the problem.

!!! note "Aliasing semantics"
    The `x`/`y`/`sigma` arrays in `sol.prob` alias the solver cache's internal
    arrays rather than copies of them. For nonlinear fits those internal arrays
    are, by default, copies of the problems input arrays. A subsequent `reinit!`
    of the cache mutates them in place, so a solution returned beforehand will
    reflect the new data. Copy the `sol.prob` fields if you need a stable
    snapshot.

    Whether the cache aliases or copies the input `x`/`y`/`sigma`/`u0` arrays can
    be controlled by passing a [`CurveFitAliasSpecifier`](@ref) as the `alias`
    keyword to `init`/`solve`. By default `x`/`y`/`sigma` are aliased and `u0` is
    copied for the linear fits, while the nonlinear fits copy all of them so that
    `reinit!` doesn't overwrite the users input.
"""
@concrete struct CurveFitSolution <: AbstractCurveFitSolution
    alg <: AbstractCurveFitAlgorithm
    u
    resid
    prob <: CurveFitProblem
    retcode::ReturnCode.T
    original
end

function CurveFitSolution(alg, coeffs, resid, prob, retcode)
    return CurveFitSolution(alg, coeffs, resid, prob, retcode, nothing)
end

function Base.show(io::IO, ::MIME"text/plain", sol::CurveFitSolution)
    alg = @something(sol.alg, sol.original.alg) |> typeof |> nameof

    println(io, "retcode: ", sol.retcode)
    if is_nonlinear_problem(sol.prob)
        println(io, "f: ", sol.prob.nlfunc.f)
    end
    println(io, "alg: $(alg)")

    mean_resid = sum(sol.resid) / length(sol.resid)
    println(io, "residuals mean: ", mean_resid)
    print(io, "u: $(sol.u)")

    return nothing
end

SciMLBase.successful_retcode(sol::CurveFitSolution) = SciMLBase.successful_retcode(sol.retcode)


"""
    CurveFitAliasSpecifier(; alias_x = nothing, alias_y = nothing, alias_sigma = nothing, alias_u0 = nothing, alias = nothing)

Holds information on what variables to alias when solving a curve fit problem.
Conforms to the [`SciMLBase.AbstractAliasSpecifier`](@extref) interface.

When a keyword argument is `nothing`, the default behaviour is used.

### Keywords

* `alias_x::Union{Bool, Nothing}`: alias the `x` array.
* `alias_y::Union{Bool, Nothing}`: alias the `y` array.
* `alias_sigma::Union{Bool, Nothing}`: alias the `sigma` array.
* `alias_u0::Union{Bool, Nothing}`: alias the `u0` array.
* `alias::Union{Bool, Nothing}`: sets all fields of the `CurveFitAliasSpecifier` to `alias`.
"""
struct CurveFitAliasSpecifier <: SciMLBase.AbstractAliasSpecifier
    alias_x::Union{Bool, Nothing}
    alias_y::Union{Bool, Nothing}
    alias_sigma::Union{Bool, Nothing}
    alias_u0::Union{Bool, Nothing}

    function CurveFitAliasSpecifier(;
            alias_x = nothing, alias_y = nothing,
            alias_sigma = nothing, alias_u0 = nothing,
            alias = nothing
        )
        return if isnothing(alias)
            new(alias_x, alias_y, alias_sigma, alias_u0)
        elseif alias
            new(true, true, true, true)
        elseif !alias
            new(false, false, false, false)
        end
    end
end

_maybe_alias(::Nothing, ::Union{Bool, Nothing}, ::Bool) = nothing

function _maybe_alias(x::AbstractArray, flag::Union{Bool, Nothing}, default::Bool)
    return @something(flag, default) ? x : copyto!(similar(x), x)
end

# Return a `CurveFitProblem` whose input arrays are aliased or copied according
# to `alias` (falling back to `default_*` when `alias = nothing`).
function _alias_inputs(
        prob::CurveFitProblem, alias::CurveFitAliasSpecifier;
        default_x::Bool = true, default_y::Bool = true,
        default_sigma::Bool = true, default_u0::Bool = false
    )
    return CurveFitProblem(
        _maybe_alias(prob.x, alias.alias_x, default_x),
        _maybe_alias(prob.y, alias.alias_y, default_y),
        _maybe_alias(prob.sigma, alias.alias_sigma, default_sigma),
        prob.nlfunc,
        _maybe_alias(prob.u0, alias.alias_u0, default_u0),
        prob.lb,
        prob.ub
    )
end

# Common Solve Interface
function CommonSolve.init(prob::AbstractCurveFitProblem; kwargs...)
    return init(
        prob,
        is_nonlinear_problem(prob) ? __FallbackNonlinearFitAlgorithm(nothing) :
            error("Default algorithm is not defined for linear problems");
        kwargs...
    )
end

function CommonSolve.init(
        prob::AbstractCurveFitProblem, alg::AbstractNonlinearAlgorithm; kwargs...
    )
    @assert is_nonlinear_problem(prob) "Nonlinear algorithm can only be used with \
                                       nonlinear problems"
    return init(prob, __FallbackNonlinearFitAlgorithm(alg); kwargs...)
end

function CommonSolve.solve!(::AbstractCurveFitCache)
    error("solve!() must be implemented by a concrete subtype of `AbstractCurveFitCache`")
end
