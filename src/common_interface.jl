abstract type AbstractCurveFitProblem end

abstract type AbstractCurveFitAlgorithm end

abstract type AbstractCurveFitSolution end

# Core Problem Types
@doc doc"""
    CurveFitProblem(
        x, y; xfun=identity, yfun=identity, yfun_inverse=inverse(yfun), nlfunc=nothing
    )

Represents a curve fitting problem where `x` and `y` are the data points to fit. It is not
recommende to use this constructor directly, instead use one of the specialized
curve fitting problem constructors like [`LinearCurveFitProblem`](@ref),
[`LogCurveFitProblem`](@ref), [`PowerCurveFitProblem`](@ref),
[`ExpCurveFitProblem`](@ref) or [`NonlinearCurveFitProblem`](@ref).
"""
@concrete struct CurveFitProblem <: AbstractCurveFitProblem
    x <: AbstractArray
    y <: AbstractArray
    xfun <: Union{Nothing, Function}
    yfun <: Union{Nothing, Function}
    yfun_inverse <: Union{Nothing, Function}
    nlfunc <: Union{Nothing, NonlinearFunction}
end

is_nonlinear_problem(prob::CurveFitProblem) = prob.nlfunc !== nothing

function CurveFitProblem(
        x, y; xfun = nothing, yfun = nothing, nlfunc = nothing
)
    if nlfunc !== nothing
        @assert xfun===nothing "Nonlinear function must have xfun = identity"
    else
        @assert xfun!==nothing "xfun must be provided for linear problems \
                                (`nlfunc` is `nothing`)"
        @assert yfun!==nothing "yfun must be provided for linear problems \
                                (`nlfunc` is `nothing`)"
        @assert ndims(x)==ndims(y)==1 "x and y must be 1-dimensional arrays for linear \
                                       problems (`nlfunc` is `nothing`)"
    end

    return CurveFitProblem(
        x, y, xfun, yfun, yfun === nothing ? nothing : inverse(yfun), nlfunc
    )
end

### Helpful aliases
"""
    LinearCurveFitProblem(x, y; xfun = identity, yfun = identity)

Represents a linear curve fitting problem where `x` and `y` are the data points to fit.
We want to solve for `a` and `b` such that:

```math
yfun(y) = a * xfun(x) + b
```

Note that this is a general problem specification of a curve fitting problem which can
be converted to a linear fit in a specific function space by choosing appropriate
`xfun` and `yfun`. The `yfun_inverse` is used to convert the fitted values back to the
original space (can be specified by defining `InverseFunctions.inverse`)
"""
function LinearCurveFitProblem(x, y; xfun = identity, yfun = identity)
    return CurveFitProblem(x, y; xfun, yfun)
end

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

"""
    RationalPolynomialFitAlgorithm(num_degree::Int, den_degree::Int)
    RationalPolynomialFitAlgorithm(;
        num_degree::Int, den_degree::Int, alg = nothing
    )

Represents a rational polynomial fitting algorithm with numerator degree `num_degree`
and denominator degree `den_degree`. The internal polynomial fitting algorithm is
determined by the `alg` keyword argument. If `alg` is `nothing` or a
`AbstractNonlinearAlgorithm` (like solvers from NonlinearSolve.jl), it will use a
nonlinear curve fitting approach. If `alg` is a `PolynomialFitAlgorithm`, it will use
a polynomial fitting approach.
"""
@kwdef @concrete struct RationalPolynomialFitAlgorithm <: AbstractCurveFitAlgorithm
    num_degree::Int
    den_degree::Int
    alg <: Union{Nothing, PolynomialFitAlgorithm, AbstractNonlinearAlgorithm} = nothing
end

function RationalPolynomialFitAlgorithm(num_degree::Int, den_degree::Int)
    return RationalPolynomialFitAlgorithm(num_degree, den_degree, nothing)
end

## Internal types for dispatch
struct __FallbackLinearFitAlgorithm <: AbstractCurveFitAlgorithm end
struct __FallbackNonlinearFitAlgorithm <: AbstractCurveFitAlgorithm end

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
    prob <: CurveFitProblem
end

# Common Solve Interface
function CommonSolve.solve(prob::AbstractCurveFitProblem; kwargs...)
    return solve(
        prob,
        is_nonlinear_problem(prob) ? __FallbackNonlinearFitAlgorithm() : __FallbackLinearFitAlgorithm();
        kwargs...
    )
end
