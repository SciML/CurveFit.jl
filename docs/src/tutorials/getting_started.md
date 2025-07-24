# [Getting Started](@id getting_started)

This tutorial will guide you through the basic usage of CurveFit.jl, demonstrating the core concepts and common workflows.

## Installation

First, install CurveFit.jl using the Julia package manager:

```julia
using Pkg
Pkg.add("CurveFit")
```

## Basic Concepts

CurveFit.jl follows the SciML common interface pattern. The basic workflow is:

1. Create a `CurveFitProblem` with your data and algorithm choice
2. Call `solve` to obtain a `CurveFitSolution`
3. Extract and use the fitted parameters from the solution

## Linear Fitting Example

Let's start with a simple linear regression:

```julia
using CurveFit
using Plots

# Generate synthetic data with noise
x = 0:0.1:10
y_true = @. 2.5 * x + 3.0
y = y_true + 0.5 * randn(length(x))

# Create and solve the fitting problem
prob = CurveFitProblem(x, y, LinearCurveFitAlgorithm())
sol = solve(prob)

# Extract fitted parameters
println("Fitted parameters: ", sol.u)
println("Converged: ", sol.retcode)

# Plot the results
scatter(x, y, label="Data", markersize=3, alpha=0.6)
plot!(x, sol(x), label="Fitted line", linewidth=2)
plot!(x, y_true, label="True line", linewidth=2, linestyle=:dash)
```

## Polynomial Fitting

CurveFit.jl can fit polynomials of any degree:

```julia
# Generate data from a quadratic function
x = -2:0.1:2
y_true = @. 1.5 * x^2 - 2.0 * x + 0.5
y = y_true + 0.1 * randn(length(x))

# Fit a second-degree polynomial
prob = CurveFitProblem(x, y, PolynomialFitAlgorithm(2))
sol = solve(prob)

println("Polynomial coefficients: ", sol.u)

# Evaluate and plot
y_fitted = sol(x)
scatter(x, y, label="Data", markersize=3)
plot!(x, y_fitted, label="Fitted polynomial", linewidth=2)
```

## Exponential Fitting

For exponential relationships:

```julia
# Generate exponential data
x = 0:0.1:5
y_true = @. 2.0 * exp(0.5 * x)
y = y_true + 0.1 * y_true .* randn(length(x))  # Multiplicative noise

# Fit exponential model
prob = CurveFitProblem(x, y, ExpCurveFitAlgorithm())
sol = solve(prob)

println("Exponential parameters: ", sol.u)

# Plot
scatter(x, y, label="Data", markersize=3)
plot!(x, sol(x), label="Fitted exponential", linewidth=2)
```

## Logarithmic and Power Law Fitting

```julia
# Logarithmic fit
x_log = 1:0.1:10
y_log = @. 2.0 * log(x_log) + 1.0 + 0.1 * randn()

prob_log = CurveFitProblem(x_log, y_log, LogCurveFitAlgorithm())
sol_log = solve(prob_log)

# Power law fit
x_power = 1:0.1:10
y_power = @. 3.0 * x_power^0.7 + 0.2 * randn()

prob_power = CurveFitProblem(x_power, y_power, PowerCurveFitAlgorithm())
sol_power = solve(prob_power)

# Plot both
p1 = scatter(x_log, y_log, label="Data")
plot!(p1, x_log, sol_log(x_log), label="Log fit")

p2 = scatter(x_power, y_power, label="Data")
plot!(p2, x_power, sol_power(x_power), label="Power fit")

plot(p1, p2, layout=(1,2))
```

## Rational Polynomial Fitting

For more complex functions, rational polynomials can provide excellent approximations:

```julia
# Generate data from a rational function
x = -2:0.1:2
y_true = @. (1 + 2*x) / (1 - x + x^2)
y = y_true + 0.05 * randn(length(x))

# Fit with rational polynomial (degree 1 numerator, degree 2 denominator)
prob = CurveFitProblem(x, y, RationalPolynomialFitAlgorithm(1, 2))
sol = solve(prob)

println("Rational polynomial coefficients: ", sol.u)

scatter(x, y, label="Data", markersize=3)
plot!(x, sol(x), label="Rational fit", linewidth=2)
```

## General Nonlinear Fitting

For custom models, use the NonlinearCurveFitProblem:

```julia
# Define a custom model function
model(x, p) = @. p[1] * sin(p[2] * x + p[3]) + p[4]

# Generate data
x = 0:0.1:10
p_true = [2.0, 1.5, 0.5, 1.0]
y = model(x, p_true) + 0.1 * randn(length(x))

# Initial guess for parameters
p0 = [1.0, 1.0, 0.0, 0.0]

# Create nonlinear fitting problem
prob = NonlinearCurveFitProblem(model, x, y, p0)
sol = solve(prob)

println("Fitted parameters: ", sol.u)
println("True parameters: ", p_true)

# Plot
scatter(x, y, label="Data", markersize=3)
plot!(x, sol(x), label="Fitted model", linewidth=2)
```

## Working with Solutions

The solution object provides various useful information:

```julia
# From any fitting solution
println("Converged: ", sol.retcode == ReturnCode.Success)
println("Parameters: ", sol.u)

# Evaluate at new points
x_new = range(minimum(x), maximum(x), length=100)
y_pred = sol(x_new)

# For some algorithms, access residuals
if hasproperty(sol, :resid)
    println("Residual norm: ", norm(sol.resid))
end
```

## Tips for Success

1. **Choose the right algorithm**: Match the algorithm to your model type
2. **Scale your data**: Normalize data for better numerical stability
3. **Check convergence**: Always verify `sol.retcode`
4. **Visualize results**: Plot your fit to check for reasonableness
5. **Try multiple algorithms**: Some problems benefit from specific methods

## Next Steps

- Explore [Linear and Polynomial Fitting](@ref) for detailed regression techniques
- Learn about [Exponential and Special Functions](@ref) for specialized models
- Check out [Nonlinear Curve Fitting](@ref) for complex custom models
- See the [API Reference](@ref) for all available algorithms and options