# API Reference

This page provides a comprehensive reference for all exported functions and types in CurveFit.jl.

## Problem Types

```@docs
CurveFitProblem
NonlinearCurveFitProblem
```

## Algorithms

### Linear and Polynomial Fitting

```@docs
LinearCurveFitAlgorithm
PolynomialFitAlgorithm
```

### Logarithmic and Power Law Fitting

```@docs
LogCurveFitAlgorithm
PowerCurveFitAlgorithm
```

### Exponential Fitting

```@docs
ExpCurveFitAlgorithm
ExpSumFitAlgorithm
```

### Rational Polynomial Fitting

```@docs
RationalPolynomialFitAlgorithm
```

### Specialized Algorithms

```@docs
KingCurveFitAlgorithm
ModifiedKingCurveFitAlgorithm
```

## Solutions

```@docs
CurveFitSolution
```

## Common Interface Functions

```@docs
solve
solve!
init
```

## Utilities

### Model Evaluation

The solution object can be called as a function to evaluate the fitted model:

```julia
sol = solve(prob)
y_pred = sol(x_new)  # Evaluate at new points
```

### Algorithm-Specific Features

Some algorithms provide additional functionality:

- **Linear algorithms**: Direct solution via QR decomposition
- **Polynomial algorithms**: Stable evaluation using orthogonal polynomials
- **Rational algorithms**: Padé approximation techniques
- **Nonlinear algorithms**: Interface to NonlinearSolve.jl

## Algorithm Selection Guide

| Data Type | Recommended Algorithm | Notes |
|-----------|----------------------|-------|
| Linear relationship | `LinearCurveFitAlgorithm()` | Fast, exact solution |
| Polynomial | `PolynomialFitAlgorithm(n)` | Specify degree `n` |
| Exponential growth/decay | `ExpCurveFitAlgorithm()` | Handles `y = a*exp(b*x)` |
| Logarithmic | `LogCurveFitAlgorithm()` | For `y = a*log(x) + b` |
| Power law | `PowerCurveFitAlgorithm()` | For `y = a*x^b` |
| Rational functions | `RationalPolynomialFitAlgorithm(m,n)` | Degrees m (num), n (den) |
| Multi-exponential | `ExpSumFitAlgorithm(n)` | Sum of n exponentials |
| Custom nonlinear | `NonlinearCurveFitProblem` | Full flexibility |

## Performance Tips

1. **Linear problems**: Always use `LinearCurveFitAlgorithm` when applicable
2. **Polynomial fitting**: Use orthogonal polynomials for numerical stability
3. **Initial guesses**: For nonlinear problems, good initial parameters are crucial
4. **Data scaling**: Normalize data to similar magnitudes before fitting

## Examples

### Basic Linear Fit

```julia
prob = CurveFitProblem(x, y, LinearCurveFitAlgorithm())
sol = solve(prob)
```

### Polynomial with Specific Degree

```julia
prob = CurveFitProblem(x, y, PolynomialFitAlgorithm(3))  # Cubic
sol = solve(prob)
```

### Custom Nonlinear Model

```julia
model(x, p) = @. p[1] * (1 - exp(-p[2] * x))
prob = NonlinearCurveFitProblem(model, x, y, p0)
sol = solve(prob)
```

## Index

```@index
```