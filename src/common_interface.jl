abstract type AbstractCurveFitProblem end

abstract type AbstractCurveFitAlgorithm end

abstract type AbstractCurveFitSolution end

# Core Problem Types
## Linear curve fitting
@doc doc"""
    LinearCurveFitProblem(x, y; xfun=identity, yfun=identity, yfun_inverse=inverse(yfun))

Represents a linear curve fitting problem where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` such that:

```math
yfun(y) = a * xfun(x) + b
```
"""
@concrete struct LinearCurveFitProblem <: AbstractCurveFitProblem
    x <: AbstractVector
    y <: AbstractVector
    xfun <: Function
    yfun <: Function
    yfun_inverse <: Function
end

function LinearCurveFitProblem(
        x, y; xfun = identity, yfun = identity, yfun_inverse = inverse(yfun)
)
    return LinearCurveFitProblem(x, y, xfun, yfun, yfun_inverse)
end

### Helpful aliases
@doc doc"""
    LogCurveFitProblem(x, y)

Represents a log curve fitting problem where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` such that:

```math
y = a * log(x) + b
```
"""
LogCurveFitProblem(x, y) = LinearCurveFitProblem(x, y; xfun = log, yfun = identity)

@doc doc"""
    PowerCurveFitProblem(x, y)
Represents a power curve fitting problem where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` such that:

```math
y = b * x^a
```

This is equivalent to a linear fit in log-log space, i.e.,

```math
log(y) = a * log(x) + log(b)
"""
PowerCurveFitProblem(x, y) = LinearCurveFitProblem(x, y; xfun = log, yfun = log)

@doc doc"""
    ExpCurveFitProblem(x, y)

Represents an exponential curve fitting problem where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` such that:

```math
y = b * exp(a * x)
```

This is equivalent to a linear fit in log-linear space, i.e.,

```math
log(y) = a * x + log(b)
```
"""
ExpCurveFitProblem(x, y) = LinearCurveFitProblem(x, y; xfun = identity, yfun = log)

## Nonlinear curve fitting

# Algorithms
@doc doc"""
    PolynomialFitAlgorithm(degree::Int)
    PolynomialFitAlgorithm(;
        degree::Int,
        linsolve_algorithm::Union{Nothing, AbstractLinearAlgorithm} = nothing
    )

Represents a polynomial fitting algorithm of degree `degree`. Only applicable to
[`LinearCurveFitProblem`](@ref)s.

!!! tip

    For ill-conditioned problems, it is recommended to use a linear solver algorithm
    such as `QRFactorization`.
"""
@kwdef @concrete struct PolynomialFitAlgorithm <: AbstractCurveFitAlgorithm
    degree::Int
    linsolve_algorithm <: Union{Nothing, AbstractLinearAlgorithm} = nothing
end

PolynomialFitAlgorithm(degree::Int) = PolynomialFitAlgorithm(degree, nothing)

# Solution Types
"""
    LinearCurveFitSolution(alg, coeffs, prob)

Represents the solution to a linear curve fitting problem. This is a callable struct and
can be used to evaluate the solution at a point. Exact evaluation mechanism depends on the
algorithm used to solve the problem.
"""
@concrete struct LinearCurveFitSolution <: AbstractCurveFitSolution
    alg <: Union{Nothing, AbstractCurveFitAlgorithm}
    coeffs
    prob <: LinearCurveFitProblem
end

# Common Solve Interface
function CommonSolve.solve(prob::AbstractCurveFitProblem; kwargs...)
    return solve(prob, nothing; kwargs...)
end
