# Curve Fitting Problems

## Overview

A curve fitting problem involves finding the parameters of a mathematical model that best describes a set of data points. CurveFit.jl provides a unified interface for solving various types of curve fitting problems.

## Problem Types

### Linear Problems

Linear curve fitting problems have the form:

```math
y = a_1 f_1(x) + a_2 f_2(x) + ... + a_n f_n(x)
```

where the parameters `a_i` appear linearly. Common examples include:
- Simple linear regression: `y = ax + b`
- Polynomial fitting: `y = a_n x^n + ... + a_1 x + a_0`

### Nonlinear Problems

Nonlinear problems involve parameters that appear nonlinearly:

```math
y = f(x, p)
```

where `p` is a vector of parameters. Examples include:
- Exponential: `y = a * exp(b * x)`
- Power law: `y = a * x^b`
- Sinusoidal: `y = a * sin(b * x + c) + d`

## Problem Formulation

### Least Squares Formulation

CurveFit.jl solves curve fitting problems by minimizing the sum of squared residuals:

```math
\min_p \sum_{i=1}^n w_i (y_i - f(x_i, p))^2
```

where:
- `(x_i, y_i)` are the data points
- `f(x, p)` is the model function
- `p` are the parameters to be optimized
- `w_i` are optional weights

### Weighted Least Squares

When data points have different uncertainties, weighted least squares can be used:

```julia
# Example with measurement uncertainties
x = [1, 2, 3, 4, 5]
y = [2.1, 4.2, 5.8, 8.1, 10.2]
σ = [0.1, 0.2, 0.15, 0.3, 0.25]  # Standard deviations

# Weights are inverse variances
weights = 1 ./ σ.^2

result = linear_fit(x, y; weights=weights)
```

## Creating Curve Fitting Problems

### Using Built-in Functions

For common curve types, use the specialized functions:

```julia
# Linear fit
result = linear_fit(x, y)

# Polynomial fit
result = poly_fit(x, y, degree)

# Exponential fit
result = exp_fit(x, y)
```

### Custom Models

For arbitrary models, define your own function:

```julia
# Define the model
function my_model(x, p)
    return p[1] .* exp.(-p[2] .* x) .* cos.(p[3] .* x .+ p[4])
end

# Initial parameter guess
p0 = [1.0, 0.1, 2.0, 0.0]

# Solve
result = curve_fit(my_model, x, y, p0)
```

## Multi-dimensional Problems

CurveFit.jl also supports fitting with multiple independent variables:

```julia
# 2D surface fitting
function surface_model(xy, p)
    x, y = xy[1, :], xy[2, :]
    return p[1] .+ p[2] .* x .+ p[3] .* y .+ p[4] .* x .* y
end

# Data: rows are dimensions, columns are data points
xy_data = [x_vals'; y_vals']
z_data = z_vals

result = curve_fit(surface_model, xy_data, z_data, p0)
```

## Constraints and Bounds

### Parameter Bounds

Constrain parameters to specific ranges:

```julia
# Fit with bounds
lower = [0.0, -Inf, 0.0]  # Lower bounds for each parameter
upper = [10.0, Inf, 2π]    # Upper bounds

result = curve_fit(model, x, y, p0; lower=lower, upper=upper)
```

### Linear Constraints

For problems with linear constraints on parameters:

```julia
# Example: p[1] + p[2] = 1
A = [1 1 0]
b = [1.0]

result = curve_fit(model, x, y, p0; A_eq=A, b_eq=b)
```

## Choosing Initial Parameters

Good initial parameters are crucial for nonlinear fitting:

```julia
# Strategy 1: Use domain knowledge
# For exponential decay, estimate from data
y_0 = y[1]
y_end = y[end]
decay_estimate = -log(y_end / y_0) / (x[end] - x[1])
p0 = [y_0, decay_estimate]

# Strategy 2: Use linearization
# For y = a * x^b, take logarithms
log_y = log.(y)
log_x = log.(x)
linear_result = linear_fit(log_x, log_y)
b_estimate = linear_result.params[1]
a_estimate = exp(linear_result.params[2])
p0 = [a_estimate, b_estimate]

# Strategy 3: Grid search
# Try multiple starting points
best_result = nothing
best_rmse = Inf

for a in [0.1, 1.0, 10.0], b in [-1.0, 0.0, 1.0]
    p0 = [a, b]
    try
        result = curve_fit(model, x, y, p0)
        if result.rmse < best_rmse
            best_rmse = result.rmse
            best_result = result
        end
    catch
        continue
    end
end
```

## Troubleshooting Common Issues

### Convergence Problems

If the fitting doesn't converge:

1. **Check initial parameters**: Are they reasonable?
2. **Scale your data**: Normalize to similar magnitudes
3. **Simplify the model**: Start with fewer parameters
4. **Increase iterations**: Use `maxiter` option

```julia
# Example with options
result = curve_fit(model, x, y, p0; 
    maxiter = 1000,
    tol = 1e-8,
    verbose = true
)
```

### Numerical Issues

For poorly conditioned problems:

```julia
# Normalize data
x_norm = (x .- mean(x)) ./ std(x)
y_norm = (y .- mean(y)) ./ std(y)

# Fit normalized data
result_norm = curve_fit(model_norm, x_norm, y_norm, p0)

# Transform parameters back
# (transformation depends on your model)
```

## See Also

- [Fitting Functions](@ref) - Available fitting functions
- [Solutions and Results](@ref) - Understanding fit results
- [Goodness of Fit](@ref) - Evaluating fit quality