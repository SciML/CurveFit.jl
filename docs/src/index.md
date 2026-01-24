# CurveFit.jl

CurveFit provides a unified and extensible interface for linear, nonlinear, and
specialized curve fitting in Julia. It offers built-in solvers for common linear
and special-function models, while general nonlinear curve fitting is handled
through nonlinear least squares methods from NonlinearSolve.jl.

Curve fitting problems are defined in a consistent problem–solver style, allowing
flexible solver selection and access to common statistical diagnostics such as
residuals, standard errors, and confidence intervals via the StatsAPI.jl interface.

## Installation

```julia
using Pkg
Pkg.add("CurveFit")
```

## Quick start

```@example quick_start
using CurveFit

# Sample data
x = 0:0.1:10
y = @. 2x + 1

# Create and solve the problem
prob = CurveFitProblem(x, y)
sol = solve(prob, LinearCurveFitAlgorithm())

# Check the fitted coefficients
println("Slope (a): ", sol.u[1])
println("Intercept (b): ", sol.u[2])

# Make predictions
println("Prediction at x=5: ", sol(5.0))
```

See [Getting started](@ref) for more step-by-step examples with common fits.
