# Fitting Functions

CurveFit.jl provides a variety of built-in fitting functions for common curve types, as well as a general framework for custom models.

## Built-in Fitting Functions

### Linear Fitting

```@docs
linear_fit
```

The `linear_fit` function performs simple linear regression:

```julia
# Fit y = ax + b
result = linear_fit(x, y)
a, b = result.params

# With weights
result = linear_fit(x, y; weights=w)

# Force through origin (b = 0)
result = linear_fit(x, y; intercept=false)
```

### Polynomial Fitting

```@docs
poly_fit
```

Fit data to a polynomial of specified degree:

```julia
# Quadratic fit: y = a₂x² + a₁x + a₀
result = poly_fit(x, y, 2)
coeffs = result.params  # [a₂, a₁, a₀]

# Higher-order polynomial
result = poly_fit(x, y, 5)

# Weighted polynomial fit
result = poly_fit(x, y, 3; weights=w)
```

### Exponential Fitting

```@docs
exp_fit
```

Fit exponential models of the form `y = a * exp(b * x)`:

```julia
# Basic exponential fit
result = exp_fit(x, y)
a, b = result.params

# Exponential with offset: y = a * exp(b * x) + c
result = exp_fit(x, y; offset=true)
a, b, c = result.params

# Constrained exponential (e.g., decay with b < 0)
result = exp_fit(x, y; upper=[Inf, 0.0])
```

### Logarithmic Fitting

```@docs
log_fit
```

Fit logarithmic models:

```julia
# Natural logarithm: y = a * ln(x) + b
result = log_fit(x, y)

# Base-10 logarithm: y = a * log₁₀(x) + b
result = log_fit(x, y; base=10)

# General base: y = a * log_b(x) + c
result = log_fit(x, y; base=2)
```

### Power Law Fitting

```@docs
power_fit
```

Fit power law relationships `y = a * x^b`:

```julia
# Basic power law
result = power_fit(x, y)
a, b = result.params

# Power law with offset: y = a * x^b + c
result = power_fit(x, y; offset=true)

# Constrained exponent (e.g., b > 0)
result = power_fit(x, y; lower=[-Inf, 0.0])
```

## General Curve Fitting

```@docs
curve_fit
```

For arbitrary nonlinear models:

```julia
# Define your model
model(x, p) = p[1] .* sin.(p[2] .* x .+ p[3]) .+ p[4]

# Initial parameters
p0 = [1.0, 1.0, 0.0, 0.0]

# Fit
result = curve_fit(model, x, y, p0)
```

### Advanced Options

All fitting functions support various options:

```julia
result = curve_fit(model, x, y, p0;
    # Optimization options
    maxiter = 1000,          # Maximum iterations
    tol = 1e-8,             # Convergence tolerance
    
    # Parameter bounds
    lower = [-Inf, 0.0, -π, -Inf],
    upper = [Inf, Inf, π, Inf],
    
    # Weights
    weights = w,
    
    # Algorithm selection
    algorithm = :lm,         # :lm (Levenberg-Marquardt) or :tr (Trust Region)
    
    # Verbosity
    verbose = true,          # Print progress
    
    # Jacobian
    autodiff = true,         # Use automatic differentiation
    jacobian = my_jacobian,  # Or provide custom Jacobian
)
```

## Model Functions Best Practices

### Writing Efficient Model Functions

```julia
# Good: Vectorized operations
model_good(x, p) = p[1] .* exp.(p[2] .* x) .+ p[3]

# Bad: Loop-based (slower)
function model_bad(x, p)
    y = similar(x)
    for i in 1:length(x)
        y[i] = p[1] * exp(p[2] * x[i]) + p[3]
    end
    return y
end
```

### Handling Multiple Variables

For models with multiple independent variables:

```julia
# 2D Gaussian
function gaussian_2d(xy, p)
    x, y = xy[1, :], xy[2, :]
    A, x0, y0, σx, σy = p
    return A .* exp.(
        -((x .- x0).^2 ./ (2σx^2) .+ (y .- y0).^2 ./ (2σy^2))
    )
end

# Usage
xy_data = [x_coords'; y_coords']
result = curve_fit(gaussian_2d, xy_data, z_data, p0)
```

### Providing Jacobians

For better performance with complex models:

```julia
# Model function
model(x, p) = p[1] .* (1 .- exp.(-p[2] .* x))

# Jacobian function
function model_jacobian(x, p)
    J = zeros(length(x), length(p))
    J[:, 1] = 1 .- exp.(-p[2] .* x)
    J[:, 2] = p[1] .* x .* exp.(-p[2] .* x)
    return J
end

# Use with curve_fit
result = curve_fit(model, x, y, p0; jacobian=model_jacobian)
```

## Specialized Fitting Scenarios

### Piecewise Functions

```julia
# Piecewise linear function
function piecewise_linear(x, p)
    # p = [a1, b1, a2, b2, x_break]
    y = similar(x)
    for i in 1:length(x)
        if x[i] < p[5]
            y[i] = p[1] * x[i] + p[2]
        else
            y[i] = p[3] * x[i] + p[4]
        end
    end
    return y
end

# Ensure continuity at break point
function constrained_piecewise(x, p)
    # p = [a1, b1, a2, x_break]
    # b2 calculated to ensure continuity
    b2 = p[1] * p[4] + p[2] - p[3] * p[4]
    return piecewise_linear(x, [p[1], p[2], p[3], b2, p[4]])
end
```

### Periodic Functions

```julia
# Fourier series (first 3 terms)
function fourier_series(x, p)
    # p = [a0, a1, b1, a2, b2, a3, b3, ω]
    y = p[1] / 2
    for n in 1:3
        y .+= p[2n] .* cos.(n * p[8] .* x) .+ p[2n+1] .* sin.(n * p[8] .* x)
    end
    return y
end
```

### Robust Fitting

For data with outliers:

```julia
# Huber loss function for robust fitting
function huber_loss(residuals, δ=1.0)
    loss = 0.0
    for r in residuals
        if abs(r) <= δ
            loss += r^2 / 2
        else
            loss += δ * (abs(r) - δ/2)
        end
    end
    return loss
end

# Custom robust fitting
function robust_fit(model, x, y, p0; δ=1.0)
    # Minimize Huber loss instead of squared residuals
    objective(p) = huber_loss(y .- model(x, p), δ)
    
    # Use optimization package
    result = optimize(objective, p0)
    return result
end
```

## Performance Tips

1. **Use in-place operations** when possible:
   ```julia
   function model!(y, x, p)
       @. y = p[1] * exp(p[2] * x) + p[3]
   end
   ```

2. **Pre-allocate arrays** for temporary calculations:
   ```julia
   function efficient_model(x, p, cache=similar(x))
       @. cache = exp(p[2] * x)
       return p[1] .* cache .+ p[3]
   end
   ```

3. **Use StaticArrays** for small parameter vectors:
   ```julia
   using StaticArrays
   p0 = @SVector [1.0, 2.0, 3.0]
   ```

## See Also

- [Curve Fitting Problems](@ref) - Problem formulation
- [Solutions and Results](@ref) - Understanding outputs
- [API Reference](@ref) - Complete function documentation