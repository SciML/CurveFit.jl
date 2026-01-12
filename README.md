# CurveFit.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/CurveFit/stable/)

[![codecov](https://codecov.io/gh/SciML/CurveFit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/CurveFit.jl)
[![Build Status](https://github.com/SciML/CurveFit.jl/workflows/CI/badge.svg)](https://github.com/SciML/CurveFit.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

CurveFit.jl is a high-performance curve fitting library for Julia that provides linear, polynomial, special function, and nonlinear least squares fitting algorithms. It is part of the [SciML](https://sciml.ai/) ecosystem and implements the common `solve` interface from [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl).

## Installation

To install CurveFit.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("CurveFit")
```

## Features

CurveFit.jl provides the following fitting algorithms:

- **Linear fitting**: `LinearCurveFitAlgorithm` - General linear fits with customizable transformations
- **Log fitting**: `LogCurveFitAlgorithm` - Fits `y = a*log(x) + b`
- **Power fitting**: `PowerCurveFitAlgorithm` - Fits `y = b*x^a`
- **Exponential fitting**: `ExpCurveFitAlgorithm` - Fits `y = b*exp(a*x)`
- **Polynomial fitting**: `PolynomialFitAlgorithm` - Fits polynomials of arbitrary degree
- **Rational polynomial fitting**: `RationalPolynomialFitAlgorithm` - Fits rational functions p(x)/q(x)
- **Sum of exponentials**: `ExpSumFitAlgorithm` - Fits `y = k + p1*exp(λ1*t) + p2*exp(λ2*t) + ...`
- **King's law**: `KingCurveFitAlgorithm`, `ModifiedKingCurveFitAlgorithm` - For hotwire anemometry
- **Nonlinear least squares**: Via `NonlinearCurveFitProblem` with any user-defined function

## Quick Start

### Linear Fit

```julia
using CurveFit

x = 1.0:10.0
y = @. 1.0 + 2.0 * x  # y = 1 + 2x

prob = CurveFitProblem(x, y)
sol = solve(prob, LinearCurveFitAlgorithm())

sol.u  # (2.0, 1.0) - coefficients (slope, intercept)
sol(5.0)  # Evaluate fitted curve at x=5
```

### Polynomial Fit

```julia
using CurveFit

x = 1.0:10.0
y = @. 1.0 + 2.0*x + 0.5*x^2

prob = CurveFitProblem(x, y)
sol = solve(prob, PolynomialFitAlgorithm(degree=2))

sol.u  # [1.0, 2.0, 0.5] - polynomial coefficients
sol(3.0)  # Evaluate at x=3
```

### Nonlinear Fit

```julia
using CurveFit

x = 1.0:10.0
fn(a, x) = @. a[1] + a[2] * x^a[3]  # Nonlinear model
y = fn([3.0, 2.0, 0.7], x)  # Generate data

prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], x, y)  # Initial guess
sol = solve(prob)

sol.u  # Fitted parameters ≈ [3.0, 2.0, 0.7]
sol(5.0)  # Evaluate at x=5
```

### Statistical Analysis

CurveFit.jl provides statistical functions compatible with [StatsAPI.jl](https://github.com/JuliaStats/StatsAPI.jl):

```julia
using CurveFit

# ... after fitting ...
coef(sol)        # Fitted coefficients
residuals(sol)   # Residuals
fitted(sol)      # Fitted values
rss(sol)         # Residual sum of squares
stderror(sol)    # Standard errors
confint(sol)     # Confidence intervals
```

## Documentation

For more details, see the [documentation](https://docs.sciml.ai/CurveFit/stable/).
