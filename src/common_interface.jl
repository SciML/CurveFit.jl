abstract type AbstractCurveFitProblem end

abstract type AbstractCurveFitAlgorithm end

abstract type AbstractCurveFitSolution end

# Core Problem Types
@doc doc"""
    CurveFitProblem(
        x, y;
        xfun=identity, yfun=identity, yfun_inverse=inverse(yfun), nlfunc=nothing,
        u0=nothing
    )

Represents a curve fitting problem where `x` and `y` are the data points to fit. It is not
recommende to use this constructor directly, instead use one of the specialized
curve fitting problem constructors like [`LinearCurveFitProblem`](@ref),
[`LogCurveFitProblem`](@ref), [`PowerCurveFitProblem`](@ref),
[`ExpCurveFitProblem`](@ref) or [`NonlinearCurveFitProblem`](@ref).

Certain algorithms may require an initial guess `u0` for the coefficients to fit. See
specific solver documentation for more details.
"""
@concrete struct CurveFitProblem <: AbstractCurveFitProblem
    x <: AbstractArray
    y <: Union{AbstractArray, Nothing}
    xfun <: Union{Nothing, Function}
    yfun <: Union{Nothing, Function}
    yfun_inverse <: Union{Nothing, Function}
    nlfunc <: Union{Nothing, NonlinearFunction}
    u0 <: Union{Nothing, AbstractArray}
end

is_nonlinear_problem(prob::CurveFitProblem) = prob.nlfunc !== nothing

function CurveFitProblem(
        x, y; xfun = nothing, yfun = nothing, nlfunc = nothing, u0 = nothing
)
    if nlfunc !== nothing
        @assert xfun===nothing "Nonlinear function must have xfun = identity"
        @assert yfun===nothing "Nonlinear function must have yfun = identity"
    else
        @assert y isa AbstractArray "y must be an array for linear problems"
        @assert ndims(x)==ndims(y)==1 "x and y must be 1-dimensional arrays for linear \
                                       problems (`nlfunc` is `nothing`)"
        xfun === nothing && (xfun = identity)
        yfun === nothing && (yfun = identity)
    end

    return CurveFitProblem(
        x, y, xfun, yfun, yfun === nothing ? nothing : inverse(yfun), nlfunc, u0
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
original space (can be specified by defining `InverseFunctions.inverse`).
"""
function LinearCurveFitProblem(x, y; xfun = identity, yfun = identity, u0 = nothing)
    return CurveFitProblem(x, y; xfun, yfun, u0)
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
function NonlinearCurveFitProblem(f::NonlinearFunction, u0, x, y)
    return CurveFitProblem(x, y; nlfunc = f, u0 = u0)
end
function NonlinearCurveFitProblem(f::F, u0, x, y) where {F}
    return NonlinearCurveFitProblem(NonlinearFunction(f), x, y, u0)
end

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
struct __FallbackNonlinearFitAlgorithm <: AbstractCurveFitAlgorithm end

# Solution Types
"""
    CurveFitSolution(alg, coeffs, prob)

Represents the solution to a curve fitting problem. This is a callable struct and
can be used to evaluate the solution at a point. Exact evaluation mechanism depends on the
algorithm used to solve the problem.
"""
@concrete struct CurveFitSolution <: AbstractCurveFitSolution
    alg <: AbstractCurveFitAlgorithm
    coeffs
    prob <: CurveFitProblem
    retcode::ReturnCode.T
end

# Common Solve Interface
function CommonSolve.solve(prob::AbstractCurveFitProblem; kwargs...)
    return solve(
        prob,
        is_nonlinear_problem(prob) ? __FallbackNonlinearFitAlgorithm() :
        __FallbackLinearFitAlgorithm();
        kwargs...
    )
end
