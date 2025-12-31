```@meta
CurrentModule = CurveFit
```
# API Reference

This page contains the documenation for all public API exported by CurveFit.

## Problem Types

```@docs
CurveFitProblem
NonlinearCurveFitProblem
```

## Algorithms

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
```

## Solutions

```@docs
CurveFitSolution
```

## Common Interface

```@docs
CommonSolve.solve
CommonSolve.solve!
CommonSolve.solve!(::AbstractCurveFitCache)
CommonSolve.init
CommonSolve.init(::AbstractCurveFitProblem)
```

## StatsAPI Interface

```@docs
StatsAPI.coef
StatsAPI.residuals
StatsAPI.dof
StatsAPI.dof_residual
StatsAPI.nobs
StatsAPI.predict
StatsAPI.fitted
mse
StatsAPI.rss
StatsAPI.vcov
StatsAPI.stderror
StatsAPI.confint
```