# CurveFit.jl

Linear, special function, and nonlinear curve fitting in Julia.

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

## See next

- [Basic overview](@ref) – Full documentation of problem types, algorithms, and usage details.  
- [Getting started with CurveFit.jl](@ref) – Step-by-step examples to get started with common fits.  