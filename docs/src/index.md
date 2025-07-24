# CurveFit.jl

[![Build Status](https://github.com/SciML/CurveFit.jl/workflows/CI/badge.svg)](https://github.com/SciML/CurveFit.jl/actions)
[![Coverage](https://codecov.io/gh/SciML/CurveFit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/CurveFit.jl)

Linear, special function, and nonlinear curve fitting in Julia with high-performance and robustness.

## Features

CurveFit.jl provides a unified interface for various curve fitting algorithms, following the SciML common interface:

- **Linear Curve Fitting**: Linear and polynomial regression
- **Logarithmic and Power Law Fitting**: Log and power transformations
- **Exponential Fitting**: Single and multi-exponential models
- **Rational Polynomial Fitting**: Rational function approximation
- **King's Method**: Specialized fitting for certain astrophysical profiles
- **General Nonlinear Fitting**: Interface to NonlinearSolve.jl

## Installation

To install CurveFit.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("CurveFit")
```

## Quick Example

```julia
using CurveFit

# Generate some example data
x = 0:0.1:10
y = @. 2.5 * x + 1.2 + 0.1 * randn()

# Create a linear fitting problem
prob = CurveFitProblem(x, y, LinearCurveFitAlgorithm())

# Solve the problem
sol = solve(prob)

# Access the fitted parameters
println("Fitted parameters: ", sol.u)
```

## Getting Help

- See the [Getting Started](@ref getting_started) guide for a tutorial introduction
- Check the [API Reference](@ref) for detailed function documentation
- Browse the [Tutorials](@ref) for specific use cases

## Citation

If you use CurveFit.jl in your research, please cite:

```bibtex
@software{curvefit_jl,
  author = {Paulo José Saiz Jabardo and Avik Pal and contributors},
  title = {CurveFit.jl: Curve Fitting for Julia},
  url = {https://github.com/SciML/CurveFit.jl},
  version = {v1.0.0},
  year = {2024}
}
```

## Contributing

Contributions are welcome! Please see the [SciML contribution guidelines](https://sciml.ai/contributing.html) for more information.