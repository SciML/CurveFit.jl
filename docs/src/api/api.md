# API Reference

This page provides a comprehensive reference for all exported functions and types in CurveFit.jl.

## Fitting Functions

### Linear Fitting

```@docs
linear_fit
linear_fit!
```

### Polynomial Fitting

```@docs
poly_fit
poly_fit!
polyval
```

### Exponential Fitting

```@docs
exp_fit
exp_fit!
```

### Logarithmic Fitting

```@docs
log_fit
log_fit!
```

### Power Law Fitting

```@docs
power_fit
power_fit!
```

### General Curve Fitting

```@docs
curve_fit
curve_fit!
```

## Types

### FitResult

```@docs
FitResult
```

Fields:
- `params::Vector{Float64}`: Fitted parameters
- `covariance::Matrix{Float64}`: Parameter covariance matrix
- `residuals::Vector{Float64}`: Residuals (y - ŷ)
- `jacobian::Matrix{Float64}`: Jacobian at solution
- `converged::Bool`: Whether fit converged
- `iterations::Int`: Number of iterations
- `rmse::Float64`: Root mean square error
- `r_squared::Float64`: Coefficient of determination
- `chi_squared::Float64`: Chi-squared statistic
- `aic::Float64`: Akaike Information Criterion
- `bic::Float64`: Bayesian Information Criterion

### CurveFitOptions

```@docs
CurveFitOptions
```

Fields:
- `maxiter::Int`: Maximum iterations (default: 1000)
- `tol::Float64`: Convergence tolerance (default: 1e-8)
- `verbose::Bool`: Print progress (default: false)
- `algorithm::Symbol`: Algorithm choice (default: :lm)
- `autodiff::Bool`: Use automatic differentiation (default: true)

## Utility Functions

### Model Evaluation

```@docs
evaluate_model
predict
```

### Statistical Functions

```@docs
confidence_intervals
prediction_intervals
parameter_errors
correlation_matrix
```

### Goodness of Fit

```@docs
rmse
r_squared
chi_squared
reduced_chi_squared
aic
bic
```

### Residual Analysis

```@docs
residuals
standardized_residuals
studentized_residuals
cooks_distance
leverage
```

## Advanced Functions

### Weighted Fitting

```@docs
weighted_fit
set_weights
```

### Constrained Fitting

```@docs
constrained_fit
add_constraint
```

### Multi-Response Fitting

```@docs
multiresponse_fit
```

### Bootstrap Analysis

```@docs
bootstrap_fit
bootstrap_confidence_intervals
```

## Low-Level Interface

### Problem Types

```@docs
CurveFitProblem
LinearFitProblem
NonlinearFitProblem
```

### Solver Interface

```@docs
solve
solve!
```

### Jacobian Functions

```@docs
compute_jacobian
finite_difference_jacobian
forward_diff_jacobian
```

## Constants and Defaults

```@docs
DEFAULT_OPTIONS
SUPPORTED_ALGORITHMS
```

## Exceptions

```@docs
CurveFitException
ConvergenceException
SingularException
DimensionMismatch
```

## Index

```@index
Pages = ["api.md"]
```