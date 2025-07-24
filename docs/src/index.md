# CurveFit.jl

Linear, special function, and nonlinear curve fitting in Julia.

## Installation

```julia
using Pkg
Pkg.add("CurveFit")
```

## Example

```julia
using CurveFit

# Linear fitting
x = 0:0.1:10
y = @. 2.5 * x + 3.0 + 0.1 * randn()

prob = CurveFitProblem(x, y, LinearCurveFitAlgorithm())
sol = solve(prob)
```

## API

See the [API](@ref) page for detailed documentation of all exported functions.