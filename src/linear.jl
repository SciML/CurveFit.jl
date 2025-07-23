"""
    linear_fit(x, y; weights=nothing, intercept=true)

Perform linear regression to fit the model `y = ax + b`.

# Arguments
- `x`: Independent variable data
- `y`: Dependent variable data
- `weights`: Optional weights for weighted least squares
- `intercept`: Whether to include intercept term (default: true)

# Returns
- `FitResult`: Object containing fitted parameters, statistics, and methods

# Example
```julia
x = 1:10
y = 2.5 .* x .+ 1.2 .+ randn(10)
result = linear_fit(x, y)
a, b = result.params  # slope and intercept
```
"""
function linear_fit(x::AbstractVector, y::AbstractVector; 
                   weights=nothing, intercept=true)
    @assert length(x) == length(y) "x and y must have the same length"
    
    n = length(x)
    
    # Build design matrix
    if intercept
        X = [x ones(n)]
        n_params = 2
    else
        X = reshape(x, :, 1)
        n_params = 1
    end
    
    # Apply weights if provided
    if weights !== nothing
        @assert length(weights) == n "weights must have the same length as x"
        W = Diagonal(sqrt.(weights))
        X_weighted = W * X
        y_weighted = W * y
    else
        X_weighted = X
        y_weighted = y
    end
    
    # Solve using QR decomposition for numerical stability
    Q, R = qr(X_weighted)
    params = R \ (Q' * y_weighted)
    
    # Calculate residuals and statistics
    y_pred = X * params
    residuals = y .- y_pred
    
    # Calculate covariance matrix
    if weights !== nothing
        s2 = sum(weights .* residuals.^2) / (n - n_params)
        covariance = s2 * inv(X' * Diagonal(weights) * X)
    else
        s2 = sum(residuals.^2) / (n - n_params)
        covariance = s2 * inv(X' * X)
    end
    
    # Calculate R-squared
    if intercept
        ss_tot = sum((y .- mean(y)).^2)
        ss_res = sum(residuals.^2)
        r_squared = 1 - ss_res / ss_tot
    else
        # For no-intercept model, use uncentered R-squared
        r_squared = 1 - sum(residuals.^2) / sum(y.^2)
    end
    
    # Create evaluation function
    if intercept
        eval_fn = x_new -> params[1] .* x_new .+ params[2]
    else
        eval_fn = x_new -> params[1] .* x_new
    end
    
    return FitResult(
        params = params,
        covariance = covariance,
        residuals = residuals,
        jacobian = X,
        converged = true,
        iterations = 1,  # Direct solution
        rmse = sqrt(mean(residuals.^2)),
        r_squared = r_squared,
        chi_squared = sum(residuals.^2),
        aic = n * log(2π) + n * log(s2) + n + 2 * n_params,
        bic = n * log(2π) + n * log(s2) + n + n_params * log(n),
        eval = eval_fn
    )
end

"""
    linear_fit!(x, y; kwargs...)

In-place version of `linear_fit`. Modifies input arrays during computation.
"""
function linear_fit!(x::AbstractVector, y::AbstractVector; kwargs...)
    # For linear fitting, in-place doesn't provide significant benefits
    # but we provide it for API consistency
    return linear_fit(x, y; kwargs...)
end