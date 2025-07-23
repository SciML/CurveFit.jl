# CurveFit.jl

High-performance curve fitting library for Julia, providing easy-to-use methods for fitting data to various mathematical models.

## Features

- **Linear Regression**: Simple and multiple linear regression
- **Polynomial Fitting**: Fit data to polynomials of any degree
- **Exponential Models**: Exponential growth and decay fitting
- **Logarithmic Fitting**: Natural and base-10 logarithmic models
- **Power Law Fitting**: Power law relationships
- **Custom Models**: Framework for fitting arbitrary nonlinear models

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
x = 1:10
y = 2.5 .* x .+ 1.2 .+ 0.1 .* randn(10)

# Perform linear fit
result = linear_fit(x, y)

# Access the fitted parameters
a, b = result.params  # y = a*x + b
```

## Getting Help

- See the [Getting Started](@ref getting_started) guide for a tutorial introduction
- Check the [API Reference](@ref) for detailed function documentation
- Browse the [Tutorials](@ref) for specific use cases

## Citation

If you use CurveFit.jl in your research, please cite:

```bibtex
@software{curvefit_jl,
  author = {CurveFit Contributors},
  title = {CurveFit.jl: Curve Fitting for Julia},
  url = {https://github.com/YourOrg/CurveFit.jl},
  version = {v0.1.0},
  year = {2024}
}
```

## Contributing

Contributions are welcome! Please see our [contribution guidelines](https://github.com/YourOrg/CurveFit.jl/blob/main/CONTRIBUTING.md) for more information.