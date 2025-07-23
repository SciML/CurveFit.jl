# Nonlinear Curve Fitting

This tutorial covers advanced nonlinear curve fitting scenarios, demonstrating how to handle complex models, multiple datasets, and challenging optimization problems.

## Basic Nonlinear Least Squares

Let's start with a classic nonlinear least squares problem:

```julia
using CurveFit
using Plots

# Generate data from a damped oscillator
t = 0:0.1:10
A_true = 2.0
ω_true = 2π
φ_true = π/4
τ_true = 3.0

y_true = A_true .* exp.(-t/τ_true) .* sin.(ω_true .* t .+ φ_true)
y = y_true .+ 0.1 .* randn(length(t))

# Define the model
damped_oscillator(t, p) = p[1] .* exp.(-t/p[4]) .* sin.(p[2] .* t .+ p[3])

# Initial guess
p0 = [1.5, 6.0, 0.0, 2.0]  # [A, ω, φ, τ]

# Fit the model
result = curve_fit(damped_oscillator, t, y, p0)

println("True parameters: [A=$A_true, ω=$ω_true, φ=$φ_true, τ=$τ_true]")
println("Fitted parameters: $(result.params)")
```

## Handling Multiple Minima

Some models have multiple local minima. Here's how to handle them:

```julia
# Model with multiple minima
multi_modal(x, p) = p[1] .* sin.(p[2] .* x) .+ p[3] .* cos.(p[4] .* x)

# Generate synthetic data
x = range(0, 10, length=100)
p_true = [2.0, 1.5, 1.0, 3.0]
y = multi_modal(x, p_true) .+ 0.1 .* randn(length(x))

# Try multiple starting points
n_trials = 20
results = []

for trial in 1:n_trials
    # Random initial guess
    p0 = 4 .* rand(4)
    
    try
        result = curve_fit(multi_modal, x, y, p0; verbose=false)
        push!(results, result)
    catch
        # Skip failed fits
        continue
    end
end

# Find best result
best_idx = argmin([r.rmse for r in results])
best_result = results[best_idx]

println("Best fit RMSE: $(best_result.rmse)")
println("Best parameters: $(best_result.params)")
```

## Constrained Nonlinear Fitting

Many physical models require constraints:

```julia
# Michaelis-Menten kinetics: v = Vmax * [S] / (Km + [S])
# Constraints: Vmax > 0, Km > 0
michaelis_menten(S, p) = p[1] .* S ./ (p[2] .+ S)

# Generate data
S = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
Vmax_true = 100.0
Km_true = 2.5
v = michaelis_menten(S, [Vmax_true, Km_true]) .+ 2 .* randn(length(S))

# Fit with constraints
p0 = [50.0, 1.0]
lower = [0.0, 0.0]  # Both parameters must be positive
upper = [Inf, Inf]

result = curve_fit(michaelis_menten, S, v, p0; lower=lower, upper=upper)

println("Vmax = $(result.params[1]) (true: $Vmax_true)")
println("Km = $(result.params[2]) (true: $Km_true)")

# Plot with confidence intervals
S_fine = range(0, maximum(S), length=100)
v_pred, v_lower, v_upper = predict(result, S_fine; interval=0.95)

plot(S_fine, v_pred, ribbon=(v_pred - v_lower, v_upper - v_pred),
     label="Fit with 95% CI", alpha=0.3)
scatter!(S, v, label="Data", markersize=5)
```

## Global Optimization Strategies

For difficult optimization problems:

```julia
# Complex model with many parameters
function complex_model(x, p)
    result = p[1]
    for i in 2:2:length(p)-1
        result .+= p[i] .* exp.(-((x .- p[i+1]) ./ 0.5).^2)
    end
    return result
end

# Generate data (sum of 3 Gaussians)
x = range(-2, 8, length=200)
p_true = [0.5, 2.0, 1.0, 3.0, 3.5, 1.5, 5.0]
y = complex_model(x, p_true) .+ 0.05 .* randn(length(x))

# Global optimization with simulated annealing-like approach
function global_fit(model, x, y, n_params; n_restarts=50, temp_schedule=exp)
    best_result = nothing
    best_rmse = Inf
    
    for restart in 1:n_restarts
        # Temperature for this restart
        temp = temp_schedule(-restart/10)
        
        # Random perturbation based on temperature
        p0 = randn(n_params) .* temp .+ randn(n_params)
        
        try
            result = curve_fit(model, x, y, p0; verbose=false, maxiter=500)
            if result.rmse < best_rmse
                best_rmse = result.rmse
                best_result = result
            end
        catch
            continue
        end
    end
    
    return best_result
end

result = global_fit(complex_model, x, y, 7)
```

## Robust Nonlinear Fitting

When dealing with outliers:

```julia
# Generate data with outliers
x = range(0, 10, length=50)
y_clean = 2.0 .* exp.(-0.5 .* x) .+ 1.0
y = copy(y_clean)

# Add outliers
outlier_idx = [10, 25, 40]
y[outlier_idx] .+= 3 .* randn(length(outlier_idx))

# Standard fit (sensitive to outliers)
model(x, p) = p[1] .* exp.(p[2] .* x) .+ p[3]
result_standard = curve_fit(model, x, y, [1.0, -0.3, 0.5])

# Robust fit using iteratively reweighted least squares
function robust_curve_fit(model, x, y, p0; n_iter=5)
    current_p = p0
    weights = ones(length(x))
    
    for iter in 1:n_iter
        # Fit with current weights
        result = curve_fit(model, x, y, current_p; weights=weights)
        current_p = result.params
        
        # Update weights based on residuals
        residuals = y .- model(x, current_p)
        mad = median(abs.(residuals .- median(residuals)))
        
        # Tukey's biweight
        c = 4.685 * mad
        weights = (abs.(residuals) .< c) .* (1 .- (residuals ./ c).^2).^2
    end
    
    return curve_fit(model, x, y, current_p; weights=weights)
end

result_robust = robust_curve_fit(model, x, y, [1.0, -0.3, 0.5])

# Compare results
plot(x, y_clean, label="True function", linewidth=2)
scatter!(x, y, label="Data with outliers", markersize=4)
scatter!(x[outlier_idx], y[outlier_idx], label="Outliers", color=:red, markersize=6)
plot!(x, model(x, result_standard.params), label="Standard fit", linewidth=2)
plot!(x, model(x, result_robust.params), label="Robust fit", linewidth=2, linestyle=:dash)
```

## Multi-Response Fitting

Fitting multiple related datasets simultaneously:

```julia
# Model for enzyme kinetics at different temperatures
function enzyme_model(S, T, p)
    # p = [Vmax_ref, Km_ref, Ea_V, Ea_K]
    R = 8.314  # Gas constant
    T_ref = 298.15  # Reference temperature (25°C)
    
    Vmax = p[1] * exp(-p[3]/R * (1/T - 1/T_ref))
    Km = p[2] * exp(-p[4]/R * (1/T - 1/T_ref))
    
    return Vmax .* S ./ (Km .+ S)
end

# Generate data at multiple temperatures
temperatures = [293.15, 298.15, 303.15, 308.15]  # 20, 25, 30, 35°C
S_values = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

# True parameters
p_true = [100.0, 2.0, 20000.0, 10000.0]

# Collect all data
all_S = Float64[]
all_T = Float64[]
all_v = Float64[]

for T in temperatures
    for S in S_values
        push!(all_S, S)
        push!(all_T, T)
        v = enzyme_model(S, T, p_true) + 2*randn()
        push!(all_v, v)
    end
end

# Define wrapper for curve_fit
function multi_temp_model(data, p)
    S, T = data[1, :], data[2, :]
    return [enzyme_model(S[i], T[i], p) for i in 1:length(S)]
end

# Fit all data simultaneously
data = [all_S'; all_T']
p0 = [80.0, 3.0, 15000.0, 8000.0]
result = curve_fit(multi_temp_model, data, all_v, p0)

println("Fitted parameters: $(result.params)")
println("True parameters: $p_true")

# Visualize fits
colors = [:blue, :red, :green, :orange]
for (i, T) in enumerate(temperatures)
    idx = all_T .== T
    scatter!(S_values, all_v[idx], label="Data $(T-273.15)°C", color=colors[i])
    
    S_fine = range(0, 10, length=100)
    v_pred = [enzyme_model(s, T, result.params) for s in S_fine]
    plot!(S_fine, v_pred, label="Fit $(T-273.15)°C", color=colors[i])
end
xlabel!("Substrate concentration [S]")
ylabel!("Reaction rate v")
```

## Advanced Jacobian Handling

For better performance with complex models:

```julia
# Model with analytical Jacobian
function exp_decay_model(t, p)
    return p[1] .* exp.(p[2] .* t) .+ p[3]
end

function exp_decay_jacobian(t, p)
    n = length(t)
    J = zeros(n, 3)
    
    exp_term = exp.(p[2] .* t)
    J[:, 1] = exp_term
    J[:, 2] = p[1] .* t .* exp_term
    J[:, 3] .= 1.0
    
    return J
end

# Compare performance
using BenchmarkTools

t = range(0, 5, length=100)
y = 2.0 .* exp.(-0.5 .* t) .+ 1.0 .+ 0.1 .* randn(length(t))
p0 = [1.5, -0.3, 0.5]

# With automatic differentiation (default)
@btime curve_fit($exp_decay_model, $t, $y, $p0)

# With analytical Jacobian
@btime curve_fit($exp_decay_model, $t, $y, $p0; jacobian=$exp_decay_jacobian)
```

## Model Selection and Validation

```julia
# Compare different growth models
models = Dict(
    "Exponential" => (x, p) -> p[1] .* exp.(p[2] .* x),
    "Logistic" => (x, p) -> p[1] ./ (1 .+ exp.(-p[2] .* (x .- p[3]))),
    "Gompertz" => (x, p) -> p[1] .* exp.(-exp.(p[2] .- p[3] .* x)),
    "Power" => (x, p) -> p[1] .* x .^ p[2]
)

# Initial guesses for each model
p0_dict = Dict(
    "Exponential" => [1.0, 0.1],
    "Logistic" => [10.0, 1.0, 5.0],
    "Gompertz" => [10.0, 1.0, 0.5],
    "Power" => [1.0, 0.5]
)

# Fit all models
results = Dict{String, FitResult}()
for (name, model) in models
    try
        results[name] = curve_fit(model, x, y, p0_dict[name])
    catch
        println("Failed to fit $name model")
    end
end

# Model comparison table
println("\nModel Comparison:")
println("Model        RMSE    R²      AIC     BIC")
println("-" * 45)
for (name, result) in results
    @printf("%-12s %.3f   %.3f   %.1f   %.1f\n", 
            name, result.rmse, result.r_squared, result.aic, result.bic)
end

# Cross-validation
function cross_validate(model, x, y, p0; k=5)
    n = length(x)
    fold_size = n ÷ k
    rmse_cv = Float64[]
    
    for fold in 1:k
        # Create train/test split
        test_idx = (fold-1)*fold_size+1 : min(fold*fold_size, n)
        train_idx = setdiff(1:n, test_idx)
        
        # Fit on training data
        result_train = curve_fit(model, x[train_idx], y[train_idx], p0)
        
        # Evaluate on test data
        y_pred = model(x[test_idx], result_train.params)
        rmse_fold = sqrt(mean((y[test_idx] .- y_pred).^2))
        push!(rmse_cv, rmse_fold)
    end
    
    return mean(rmse_cv), std(rmse_cv)
end

# Perform cross-validation
best_model = "Exponential"
best_cv_rmse, cv_std = cross_validate(
    models[best_model], x, y, p0_dict[best_model]
)
println("\nCross-validation RMSE for $best_model: $best_cv_rmse ± $cv_std")
```

## Tips for Successful Nonlinear Fitting

1. **Choose good initial parameters**: Use domain knowledge or preliminary analysis
2. **Scale your variables**: Normalize to similar magnitudes
3. **Check parameter identifiability**: Ensure parameters are not redundant
4. **Validate assumptions**: Check residual patterns and normality
5. **Use appropriate constraints**: Incorporate physical limits
6. **Consider global optimization**: For complex landscapes
7. **Perform sensitivity analysis**: Understand parameter importance

## See Also

- [Getting Started](@ref getting_started) - Basic usage
- [Advanced Examples](@ref) - More complex scenarios
- [API Reference](@ref) - Complete function documentation