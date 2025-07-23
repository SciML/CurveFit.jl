# [Getting Started](@id getting_started)

This tutorial will guide you through the basic usage of CurveFit.jl, demonstrating the core ideas and functionality.

## Installation

First, install CurveFit.jl using the Julia package manager:

```julia
using Pkg
Pkg.add("CurveFit")
```

## Basic Linear Fitting

Let's start with a simple example of fitting a line to data:

```julia
using CurveFit
using Plots

# Generate synthetic data with noise
x = 1:0.5:10
y_true = 2.5 .* x .+ 3.0
y = y_true .+ randn(length(x))

# Perform linear fit
result = linear_fit(x, y)

# Extract parameters
a, b = result.params
println("Fitted line: y = $(a)x + $(b)")

# Plot the results
scatter(x, y, label="Data", markersize=4)
plot!(x, a .* x .+ b, label="Fitted line", linewidth=2)
plot!(x, y_true, label="True line", linewidth=2, linestyle=:dash)
```

## Polynomial Fitting

CurveFit.jl can fit polynomials of any degree:

```julia
# Generate data from a quadratic function
x = -2:0.1:2
y_true = 1.5 .* x.^2 .- 2.0 .* x .+ 0.5
y = y_true .+ 0.1 .* randn(length(x))

# Fit a second-degree polynomial
result = poly_fit(x, y, 2)

# Extract coefficients (from highest to lowest degree)
coeffs = result.params
println("Fitted polynomial: y = $(coeffs[1])x² + $(coeffs[2])x + $(coeffs[3])")

# Evaluate the fitted polynomial
y_fitted = result.eval(x)

# Plot
scatter(x, y, label="Data", markersize=3)
plot!(x, y_fitted, label="Fitted polynomial", linewidth=2)
```

## Exponential Fitting

For exponential relationships:

```julia
# Generate exponential data
x = 0:0.1:5
y_true = 2.0 .* exp.(0.5 .* x)
y = y_true .+ 0.1 .* y_true .* randn(length(x))  # Multiplicative noise

# Fit exponential model: y = a * exp(b * x)
result = exp_fit(x, y)

a, b = result.params
println("Fitted exponential: y = $(a) * exp($(b) * x)")

# Plot
scatter(x, y, label="Data", markersize=3)
plot!(x, result.eval(x), label="Fitted exponential", linewidth=2)
```

## Custom Nonlinear Models

For more complex models, use the general curve fitting function:

```julia
# Define a custom model
model(x, p) = p[1] .* sin.(p[2] .* x .+ p[3]) .+ p[4]

# Generate data
x = 0:0.1:10
params_true = [2.0, 1.5, 0.5, 1.0]
y = model(x, params_true) .+ 0.1 .* randn(length(x))

# Initial guess for parameters
p0 = [1.0, 1.0, 0.0, 0.0]

# Fit the model
result = curve_fit(model, x, y, p0)

println("Fitted parameters: ", result.params)
println("True parameters: ", params_true)

# Plot
scatter(x, y, label="Data", markersize=3)
plot!(x, result.eval(x), label="Fitted model", linewidth=2)
```

## Assessing Fit Quality

CurveFit.jl provides various metrics to assess the quality of your fit:

```julia
# Using the linear fit from earlier
result = linear_fit(x, y)

# R-squared value
println("R² = ", result.r_squared)

# Root mean square error
println("RMSE = ", result.rmse)

# Residuals
residuals = result.residuals
println("Mean residual = ", mean(residuals))
println("Std of residuals = ", std(residuals))

# Plot residuals
scatter(x, residuals, label="Residuals", markersize=4)
hline!([0], label="", color=:black, linestyle=:dash)
xlabel!("x")
ylabel!("Residual")
```

## Weighted Fitting

When your data points have different uncertainties:

```julia
# Data with varying uncertainty
x = 1:10
y = 2 .* x .+ 1 .+ randn(10)
weights = 1 ./ (0.1 .+ 0.1 .* x)  # Higher uncertainty for larger x

# Weighted linear fit
result = linear_fit(x, y; weights=weights)

println("Weighted fit: y = $(result.params[1])x + $(result.params[2])")
```

## Next Steps

- Explore the [Linear Regression](@ref) tutorial for more details on linear models
- Learn about [Polynomial Fitting](@ref) for higher-order fits
- Check out [Nonlinear Curve Fitting](@ref) for complex models
- See the [API Reference](@ref) for all available functions