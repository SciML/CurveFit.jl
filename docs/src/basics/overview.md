# Basic overview

CurveFit provides a unified and extensible interface for linear, nonlinear, and
specialized curve fitting in Julia. It offers built-in solvers for common linear
and special-function models, while general nonlinear curve fitting is handled
through nonlinear least squares methods from NonlinearSolve.jl.

Curve fitting problems are defined in a consistent problem–solver style, allowing
flexible solver selection and access to common statistical diagnostics such as
residuals, standard errors, and confidence intervals via the StatsAPI.jl interface.
