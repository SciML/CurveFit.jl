# CurveFit.jl

Linear, special function, and nonlinear curve fitting in Julia.

## Installation

```julia
using Pkg
Pkg.add("CurveFit")
```

## Quick start

```@example quick start
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

- [Manual](manual.md) – Full documentation of problem types, algorithms, and usage details.  
- [Tutorial](tutorial.md) – Step-by-step examples to get started with common fits.  
- [API Reference](api.md) – Complete list of exported functions, algorithms, and types.