```@meta
CurrentModule = CurveFit
```
# API Reference

This page contains the documenation for all public API exported by CurveFit.

## Problem Types

Problem types define the data, model, and optional initial guesses used for curve
fitting. 

```@docs
CurveFitProblem
NonlinearCurveFitProblem
```

## Algorithms

CurveFit provides algorithms for linear fitting, nonlinear fitting with
NonlinearSolve.jl and several specialized model with dedicated solvers.

```@docs
LinearCurveFitAlgorithm
PolynomialFitAlgorithm
LogCurveFitAlgorithm
PowerCurveFitAlgorithm
ExpCurveFitAlgorithm
RationalPolynomialFitAlgorithm
KingCurveFitAlgorithm
ModifiedKingCurveFitAlgorithm
ExpSumFitAlgorithm
LM_linsolve
```

## Solutions

Solvers return a CurveFitSolution, which stores fitted parameters, residuals,
convergence information, and the original problem definition. Solutions are 
callable and integrate with StatsAPI.jl for statistical analysis.

```@docs
CurveFitSolution
```

## Common Interface

CurveFit follows the CommonSolve.jl interface, enabling a consistent workflow for
initialization, caching, and solving. These methods allow CurveFit to interoperate
with other SciML tools and solver abstractions.

```@docs
CommonSolve.solve
CommonSolve.solve!
CommonSolve.solve!(::AbstractCurveFitCache)
CommonSolve.init
CommonSolve.init(::AbstractCurveFitProblem)
```

## StatsAPI Interface

CurveFit implements the StatsAPI.jl interface, allowing its users to
examine their curve fitting solutions in more detail.

```@docs
coef
residuals
dof
dof_residual
nobs
predict
fitted
mse
rss
isconverged
vcov
stderror
confint
```