abstract type AbstractCurveFitProblem end

abstract type AbstractCurveFitAlgorithm end

abstract type AbstractCurveFitSolution end

# Core Problem Types
@doc doc"""
    CurveFitProblem(x, y; nlfunc=nothing, u0=nothing)

Represents a curve fitting problem where `x` and `y` are the data points to fit.

Certain algorithms may require an initial guess `u0` for the coefficients to fit. See
specific solver documentation for more details.

See also [`NonlinearCurveFitProblem`](@ref).
"""
@concrete struct CurveFitProblem <: AbstractCurveFitProblem
    x <: AbstractArray
    y <: Union{AbstractArray, Nothing}
    nlfunc <: Union{Nothing, NonlinearFunction}
    u0 <: Union{Nothing, AbstractArray}
end

function SciMLBase.isinplace(prob::CurveFitProblem)
    !is_nonlinear_problem(prob) && return false
    return SciMLBase.isinplace(prob.nlfunc)
end

is_nonlinear_problem(prob::CurveFitProblem) = prob.nlfunc !== nothing

function CurveFitProblem(x, y; nlfunc = nothing, u0 = nothing)
    if nlfunc === nothing
        @assert ndims(x)==ndims(y)==1 "x and y must be 1-dimensional arrays for linear \
                                       problems (`nlfunc` is `nothing`)"
    end

    return CurveFitProblem(x, y, nlfunc, u0)
end

@doc doc"""
    NonlinearCurveFitProblem(f, u0, x, y)

Nonlinear curve fitting problem where `f` is a nonlinear function to fit, `u0` is the
initial guess for the coefficients, and `x` and `y` are the data points to fit. The
following optimization problem is solved:

```math
\begin{equation}
    \underset{u}{\text{argmin}} ~ \| f(u, x) - y \|_2
\end{equation}
```

If `y` is `nothing`, then it is treated as a zero vector. `f` is a generic Julia function or
ideally a `NonlinearFunction` from [`SciMLBase.jl`](https://github.com/SciML/SciMLBase.jl).
"""
function NonlinearCurveFitProblem(f::NonlinearFunction, u0, x, y = nothing)
    return CurveFitProblem(x, y; nlfunc = f, u0 = u0)
end
function NonlinearCurveFitProblem(f::F, u0, x, y = nothing) where {F}
    return NonlinearCurveFitProblem(NonlinearFunction(f), u0, x, y)
end

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

Represents a linear curve fitting algorithm where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` such that:

```math
yfun(y) = a * xfun(x) + b
```

Note that this is a general problem specification of a curve fitting problem which can
be converted to a linear fit in a specific function space by choosing appropriate
`xfun` and `yfun`. The `yfun_inverse` is used to convert the fitted values back to the
original space (can be specified by defining `InverseFunctions.inverse`).
"""
function LinearCurveFitAlgorithm(;
        xfun = identity, yfun = identity, yfun_inverse = inverse(yfun)
)
    return LinearCurveFitAlgorithm(xfun, yfun, yfun_inverse)
end

@doc doc"""
    LogCurveFitAlgorithm()

Represents a log curve fitting algorithm where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` such that:

```math
y = a * log(x) + b
```
"""
LogCurveFitAlgorithm() = LinearCurveFitAlgorithm(; xfun = log, yfun = identity)

@doc doc"""
    PowerCurveFitAlgorithm()

Represents a power curve fitting algorithm where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` such that:

```math
y = b * x^a
```

This is equivalent to a linear fit in log-log space, i.e.,

```math
log(y) = a * log(x) + log(b)
"""
PowerCurveFitAlgorithm() = LinearCurveFitAlgorithm(; xfun = log, yfun = log)

@doc doc"""
    ExpCurveFitAlgorithm()

Represents an exponential curve fitting algorithm where `x` and `y` are the data points to
fit. We want to solve for `a` and `b` such that:

```math
y = b * exp(a * x)
```

This is equivalent to a linear fit in log-linear space, i.e.,

```math
log(y) = a * x + log(b)
```
"""
ExpCurveFitAlgorithm() = LinearCurveFitAlgorithm(; xfun = identity, yfun = log)

@doc doc"""
    KingCurveFitAlgorithm()

Represents a king curve fitting problem where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` according to original King's law (1910) that represents
the relationship between voltage (E) and velocity (U) in a hotwire anemometer:

```math
E^2 = A + B * U^0.5
```

or

```math
x^2 = A + B * y^0.5
```
"""
KingCurveFitAlgorithm() = LinearCurveFitAlgorithm(; xfun = abs2, yfun = sqrt)

@doc doc"""
    PolynomialFitAlgorithm(degree::Int)
    PolynomialFitAlgorithm(;
        degree::Int,
        linsolve_algorithm::Union{Nothing, AbstractLinearAlgorithm} = nothing
    )

Represents a polynomial fitting algorithm of degree `degree`. Only applicable to
[`LinearCurveFitAlgorithm`](@ref)s.

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
linear least squares fitting.

## Linear Rational Polynomial Fitting

In this case the following curve fit is done:

```math
y = \frac{p(x)}{q(x)}
```

where `p(x)` is a polynomial of degree `num_degree` and `q(x)` is a polynomial of degree
`den_degree`. The linear case is solved by doing a least squares fit on:

```math
y * q(x) = p(x)
```

where the zero order term of `q(x)` is assumed to be 1.

## Nonlinear Rational Polynomial Fitting

If an `u0` is not provided to the problem, then we will use linear least squares for an
initial guess.
"""
@kwdef @concrete struct RationalPolynomialFitAlgorithm <: AbstractCurveFitAlgorithm
    num_degree::Int
    den_degree::Int
    alg <: Union{Nothing, AbstractLinearAlgorithm, AbstractNonlinearAlgorithm} = nothing
end

function RationalPolynomialFitAlgorithm(num_degree::Int, den_degree::Int)
    return RationalPolynomialFitAlgorithm(num_degree, den_degree, nothing)
end

## Internal types for dispatch
struct __FallbackLinearFitAlgorithm <: AbstractCurveFitAlgorithm end
@concrete struct __FallbackNonlinearFitAlgorithm <: AbstractCurveFitAlgorithm
    alg <: Union{Nothing, AbstractNonlinearAlgorithm}
end

# Solution Types
"""
    CurveFitSolution(alg, coeffs, prob)

Represents the solution to a curve fitting problem. This is a callable struct and
can be used to evaluate the solution at a point. Exact evaluation mechanism depends on the
algorithm used to solve the problem.
"""
@concrete struct CurveFitSolution <: AbstractCurveFitSolution
    alg <: AbstractCurveFitAlgorithm
    u
    prob <: CurveFitProblem
    retcode::ReturnCode.T
    original
end

function CurveFitSolution(alg, coeffs, prob, retcode)
    return CurveFitSolution(alg, coeffs, prob, retcode, nothing)
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
