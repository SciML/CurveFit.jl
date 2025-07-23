"""
    FitResult

Type containing the results of a curve fitting operation.

# Fields
- `params::Vector{Float64}`: Fitted parameters
- `covariance::Matrix{Float64}`: Parameter covariance matrix  
- `residuals::Vector{Float64}`: Residuals (y - ŷ)
- `jacobian::Matrix{Float64}`: Jacobian matrix at solution
- `converged::Bool`: Whether the fit converged
- `iterations::Int`: Number of iterations taken
- `rmse::Float64`: Root mean square error
- `r_squared::Float64`: Coefficient of determination (R²)
- `chi_squared::Float64`: Chi-squared statistic
- `aic::Float64`: Akaike Information Criterion
- `bic::Float64`: Bayesian Information Criterion
- `eval::Function`: Function to evaluate the fitted model

# Methods
- `predict(result, x)`: Make predictions at new x values
- `confidence_intervals(result, level)`: Get parameter confidence intervals
- `parameter_errors(result)`: Get parameter standard errors
"""
struct FitResult{T<:Real}
    params::Vector{T}
    covariance::Matrix{T}
    residuals::Vector{T}
    jacobian::Matrix{T}
    converged::Bool
    iterations::Int
    rmse::T
    r_squared::T
    chi_squared::T
    aic::T
    bic::T
    eval::Function
end

"""
    CurveFitOptions

Options for curve fitting algorithms.

# Fields
- `maxiter::Int`: Maximum number of iterations (default: 1000)
- `tol::Float64`: Convergence tolerance (default: 1e-8)
- `verbose::Bool`: Whether to print progress (default: false)
- `algorithm::Symbol`: Algorithm to use (:lm or :tr, default: :lm)
- `autodiff::Bool`: Use automatic differentiation (default: true)
"""
Base.@kwdef struct CurveFitOptions
    maxiter::Int = 1000
    tol::Float64 = 1e-8
    verbose::Bool = false
    algorithm::Symbol = :lm
    autodiff::Bool = true
end

"""
    CurveFitProblem

Abstract type for curve fitting problems.
"""
abstract type CurveFitProblem end

"""
    LinearFitProblem

Problem type for linear least squares fitting.
"""
struct LinearFitProblem{T<:Real} <: CurveFitProblem
    X::Matrix{T}
    y::Vector{T}
    weights::Union{Nothing, Vector{T}}
end

"""
    NonlinearFitProblem

Problem type for nonlinear least squares fitting.
"""
struct NonlinearFitProblem{T<:Real, F<:Function, J<:Union{Nothing, Function}} <: CurveFitProblem
    model::F
    x::Union{AbstractVector, AbstractMatrix}
    y::Vector{T}
    p0::Vector{T}
    weights::Union{Nothing, Vector{T}}
    jacobian::J
    lower::Union{Nothing, Vector{T}}
    upper::Union{Nothing, Vector{T}}
end

# Exception types
struct CurveFitException <: Exception
    msg::String
end

struct ConvergenceException <: Exception
    msg::String
    iterations::Int
    final_error::Float64
end

struct SingularException <: Exception
    msg::String
end

struct DimensionMismatch <: Exception
    msg::String
    expected::Int
    actual::Int
end