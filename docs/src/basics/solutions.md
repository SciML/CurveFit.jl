# Solutions and Results

When you perform a curve fit with CurveFit.jl, the result is returned as a `FitResult` object containing comprehensive information about the fit.

## Understanding FitResult

The `FitResult` type contains all information about a completed fit:

```julia
result = linear_fit(x, y)

# Access fitted parameters
params = result.params

# Check convergence
if result.converged
    println("Fit converged in $(result.iterations) iterations")
else
    println("Fit did not converge!")
end
```

## Key Fields

### Parameters and Uncertainties

```julia
# Fitted parameters
a, b = result.params  # For linear fit: y = ax + b

# Parameter standard errors
param_errors = sqrt.(diag(result.covariance))
println("a = $(a) ± $(param_errors[1])")
println("b = $(b) ± $(param_errors[2])")

# Full covariance matrix
cov_matrix = result.covariance

# Correlation between parameters
correlation = cov_matrix[1,2] / (param_errors[1] * param_errors[2])
```

### Residuals and Predictions

```julia
# Residuals (observed - predicted)
residuals = result.residuals

# Make predictions
x_new = range(minimum(x), maximum(x), length=100)
y_pred = result.eval(x_new)

# Or use the predict function
y_pred = predict(result, x_new)

# Prediction intervals
y_pred, lower, upper = predict(result, x_new; interval=0.95)
```

### Goodness of Fit Metrics

```julia
# Root Mean Square Error
println("RMSE: $(result.rmse)")

# R-squared (coefficient of determination)
println("R²: $(result.r_squared)")

# Chi-squared statistic
println("χ²: $(result.chi_squared)")

# Reduced chi-squared (χ² per degree of freedom)
dof = length(x) - length(result.params)
reduced_chi2 = result.chi_squared / dof
println("Reduced χ²: $(reduced_chi2)")

# Information criteria
println("AIC: $(result.aic)")  # Akaike Information Criterion
println("BIC: $(result.bic)")  # Bayesian Information Criterion
```

## Confidence Intervals

### Parameter Confidence Intervals

```julia
# 95% confidence intervals for parameters
ci = confidence_intervals(result, 0.95)

for i in 1:length(result.params)
    println("Parameter $i: $(result.params[i]) ∈ [$(ci[i,1]), $(ci[i,2])]")
end
```

### Prediction Bands

```julia
# Confidence band (uncertainty in the mean)
y_fit, y_lower_conf, y_upper_conf = confidence_band(result, x_new, 0.95)

# Prediction band (uncertainty in individual predictions)
y_fit, y_lower_pred, y_upper_pred = prediction_band(result, x_new, 0.95)

# Plot with bands
plot(x_new, y_fit, label="Fit", ribbon=(y_fit - y_lower_conf, y_upper_conf - y_fit))
scatter!(x, y, label="Data")
```

## Residual Analysis

### Basic Residual Plots

```julia
using Plots

# Residuals vs fitted values
scatter(result.eval(x), result.residuals,
    xlabel="Fitted values", ylabel="Residuals",
    title="Residual Plot")
hline!([0], color=:red, linestyle=:dash, label="")

# Q-Q plot for normality check
using StatsPlots
qqplot(Normal(), result.residuals,
    title="Q-Q Plot of Residuals")
```

### Advanced Diagnostics

```julia
# Standardized residuals
std_residuals = standardized_residuals(result)

# Cook's distance (influence of each point)
cooks_d = cooks_distance(result)
influential = findall(cooks_d .> 4/length(x))
println("Influential points: $influential")

# Leverage values
leverage_vals = leverage(result)
high_leverage = findall(leverage_vals .> 2*length(result.params)/length(x))
```

## Model Comparison

When comparing different models:

```julia
# Fit multiple models
result_linear = linear_fit(x, y)
result_quad = poly_fit(x, y, 2)
result_exp = exp_fit(x, y)

# Compare using information criteria
models = ["Linear", "Quadratic", "Exponential"]
aics = [result_linear.aic, result_quad.aic, result_exp.aic]
bics = [result_linear.bic, result_quad.bic, result_exp.bic]

# Lower AIC/BIC is better
best_aic = argmin(aics)
println("Best model by AIC: $(models[best_aic])")

# Likelihood ratio test for nested models
lr_statistic = result_linear.chi_squared - result_quad.chi_squared
p_value = 1 - cdf(Chisq(1), lr_statistic)
println("Linear vs Quadratic p-value: $p_value")
```

## Extracting Detailed Information

### Jacobian and Sensitivity

```julia
# Jacobian matrix at solution
J = result.jacobian

# Parameter sensitivity
# How much does each parameter affect the predictions?
sensitivity = sum(abs.(J), dims=1) / size(J, 1)
```

### Correlation Analysis

```julia
# Parameter correlation matrix
param_corr = correlation_matrix(result)

# Visualize correlations
using Plots
heatmap(param_corr, 
    title="Parameter Correlations",
    color=:RdBu,
    clims=(-1, 1))
```

## Working with Weighted Fits

For weighted least squares:

```julia
# Perform weighted fit
weights = 1 ./ measurement_errors.^2
result = linear_fit(x, y; weights=weights)

# Weighted residuals
weighted_residuals = sqrt.(weights) .* result.residuals

# Weighted RMSE
wrmse = sqrt(mean(weighted_residuals.^2))
```

## Saving and Loading Results

```julia
using JLD2

# Save fit result
@save "fit_result.jld2" result

# Load fit result
@load "fit_result.jld2" result

# Export to DataFrame
using DataFrames

df_params = DataFrame(
    parameter = 1:length(result.params),
    value = result.params,
    std_error = sqrt.(diag(result.covariance))
)

df_fit = DataFrame(
    x = x,
    y_observed = y,
    y_fitted = result.eval(x),
    residual = result.residuals
)
```

## Custom Result Processing

```julia
# Create a summary function
function fit_summary(result::FitResult)
    println("=== Fit Summary ===")
    println("Converged: $(result.converged)")
    println("Parameters: $(result.params)")
    println("RMSE: $(round(result.rmse, digits=4))")
    println("R²: $(round(result.r_squared, digits=4))")
    
    # Parameter table
    errors = sqrt.(diag(result.covariance))
    for (i, (p, e)) in enumerate(zip(result.params, errors))
        println("  p[$i] = $(round(p, digits=4)) ± $(round(e, digits=4))")
    end
end

# Use it
fit_summary(result)
```

## Troubleshooting Non-Convergence

If `result.converged == false`:

```julia
# Check why fit didn't converge
if !result.converged
    println("Fit failed after $(result.iterations) iterations")
    
    # Examine final residuals
    println("Final RMSE: $(result.rmse)")
    println("Max residual: $(maximum(abs.(result.residuals)))")
    
    # Try with different options
    result_retry = curve_fit(model, x, y, p0;
        maxiter = 5000,
        tol = 1e-6,
        verbose = true
    )
end
```

## See Also

- [Goodness of Fit](@ref) - Detailed metrics explanation
- [API Reference](@ref) - Complete field descriptions
- [Advanced Examples](@ref) - Complex usage scenarios