# Frequently Asked Questions

## General Questions

### What is the difference between curve fitting and interpolation?

**Curve fitting** finds a smooth function that approximates your data, minimizing the overall error. The fitted curve typically doesn't pass through all data points exactly.

**Interpolation** finds a function that passes through all data points exactly. It's used when you trust your data completely and want to estimate values between known points.

Use curve fitting when:
- Your data has measurement errors
- You want to understand underlying relationships
- You need to extrapolate beyond your data range
- You want a simple model with few parameters

### How do I choose between different model types?

Consider:
1. **Physical/theoretical basis**: Does the model match the underlying process?
2. **Visual inspection**: Plot your data and see what shape it suggests
3. **Model comparison metrics**: Use AIC/BIC for statistical comparison
4. **Residual analysis**: Check if residuals are randomly distributed
5. **Parsimony**: Prefer simpler models when performance is similar

Example decision tree:
```julia
# Linear relationship?
result_linear = linear_fit(x, y)
if result_linear.r_squared > 0.95
    return result_linear
end

# Exponential growth/decay?
if all(y .> 0)  # Exponential requires positive y
    result_exp = exp_fit(x, y)
    if result_exp.aic < result_linear.aic
        return result_exp
    end
end

# Try polynomial
result_poly2 = poly_fit(x, y, 2)
result_poly3 = poly_fit(x, y, 3)

# Choose based on BIC (penalizes complexity)
# ... etc
```

### Why doesn't my fit converge?

Common causes and solutions:

1. **Poor initial parameters**
   ```julia
   # Bad: Random guess
   p0 = rand(4)
   
   # Good: Educated guess based on data
   p0 = [maximum(y), mean(x), std(y), minimum(y)]
   ```

2. **Model-data mismatch**
   ```julia
   # Check if model can represent data
   plot(x, y)  # Visual inspection
   ```

3. **Scaling issues**
   ```julia
   # Normalize data
   x_norm = (x .- mean(x)) ./ std(x)
   y_norm = (y .- mean(y)) ./ std(y)
   result = curve_fit(model_norm, x_norm, y_norm, p0)
   ```

4. **Numerical issues**
   ```julia
   # Increase tolerance and iterations
   result = curve_fit(model, x, y, p0; 
       tol=1e-6,      # Less strict
       maxiter=5000   # More iterations
   )
   ```

## Technical Questions

### How do I handle errors in both x and y?

Standard least squares assumes errors only in y. For errors in both variables:

```julia
# Total least squares / Orthogonal distance regression
function orthogonal_fit(x, y, σ_x, σ_y, model, p0)
    # Define augmented residual function
    function residuals(p_aug)
        n = length(x)
        p = p_aug[1:length(p0)]
        x_true = p_aug[length(p0)+1:end]
        
        res_y = (y .- model(x_true, p)) ./ σ_y
        res_x = (x .- x_true) ./ σ_x
        
        return vcat(res_y, res_x)
    end
    
    # Initial guess includes true x values
    p_aug_0 = vcat(p0, x)
    
    # Solve
    # ... implementation ...
end
```

### How do I fit a discontinuous function?

For piecewise functions:

```julia
# Piecewise linear with unknown break point
function piecewise_model(x, p)
    # p = [slope1, intercept1, slope2, intercept2, x_break]
    return @. ifelse(x < p[5], 
                     p[1] * x + p[2],
                     p[3] * x + p[4])
end

# Smooth approximation (better for optimization)
function smooth_piecewise(x, p, smoothness=0.1)
    # Use tanh for smooth transition
    transition = @. 0.5 * (1 + tanh((x - p[5]) / smoothness))
    left = @. p[1] * x + p[2]
    right = @. p[3] * x + p[4]
    return @. (1 - transition) * left + transition * right
end
```

### How do I enforce parameter relationships?

Use constraints or reparameterization:

```julia
# Method 1: Constraints
# Ensure p[1] + p[2] = 1
A = [1 1 0]
b = [1.0]
result = curve_fit(model, x, y, p0; A_eq=A, b_eq=b)

# Method 2: Reparameterization
# Instead of fitting p[1] and p[2] with p[1] + p[2] = 1,
# fit p[1] and compute p[2] = 1 - p[1]
function constrained_model(x, p_reduced)
    p1 = p_reduced[1]
    p2 = 1 - p1
    p3 = p_reduced[2]
    # Use p1, p2, p3 in your model...
end
```

### What's the difference between parameter errors and prediction errors?

**Parameter errors** quantify uncertainty in fitted parameters:
```julia
param_errors = sqrt.(diag(result.covariance))
println("Parameter 1: $(result.params[1]) ± $(param_errors[1])")
```

**Prediction errors** quantify uncertainty in model predictions:
```julia
# Confidence interval (uncertainty in mean prediction)
y_pred, y_lower_conf, y_upper_conf = predict(result, x_new; interval=0.95)

# Prediction interval (uncertainty in individual predictions)
y_pred, y_lower_pred, y_upper_pred = predict(result, x_new; 
    interval=0.95, prediction=true)
```

## Performance Questions

### How can I speed up fitting?

1. **Provide analytical Jacobian**:
   ```julia
   result = curve_fit(model, x, y, p0; jacobian=my_jacobian)
   ```

2. **Use StaticArrays for small problems**:
   ```julia
   using StaticArrays
   p0 = @SVector [1.0, 2.0, 3.0]
   ```

3. **Parallelize multiple fits**:
   ```julia
   using Distributed
   results = pmap(dataset -> curve_fit(model, dataset.x, dataset.y, p0), 
                  datasets)
   ```

### How much data do I need?

Rule of thumb: At least 10-20 data points per parameter, more for:
- Nonlinear models
- Noisy data  
- Correlated parameters

Check parameter uncertainties:
```julia
relative_errors = sqrt.(diag(result.covariance)) ./ abs.(result.params)
if any(relative_errors .> 0.5)
    @warn "Large parameter uncertainties - consider more data"
end
```

## Troubleshooting

### "Singular matrix" error

The fitting matrix is not invertible. Solutions:
1. Remove redundant parameters
2. Add regularization
3. Get more diverse data points
4. Check for perfect multicollinearity

### Fit looks good but metrics are poor

Check for:
1. Outliers skewing metrics
2. Heteroscedasticity (varying noise)
3. Wrong error distribution assumption

```julia
# Diagnose with residual plots
scatter(result.eval(x), result.residuals)
histogram(result.residuals)  # Should be roughly normal
```

### Different runs give different results

The problem has multiple local minima. Solutions:
1. Try multiple starting points
2. Use global optimization
3. Add constraints based on domain knowledge
4. Simplify the model

## Best Practices

1. **Always visualize**: Plot data and fit
2. **Check residuals**: Look for patterns
3. **Validate**: Use cross-validation or hold-out data
4. **Report uncertainties**: Include parameter errors
5. **Document model choice**: Explain why you chose this model
6. **Test robustness**: Try different initial conditions

## See Also

- [Getting Started](@ref getting_started) - Basic tutorial
- [Curve Fitting Problems](@ref) - Problem types
- [Solutions and Results](@ref) - Understanding outputs