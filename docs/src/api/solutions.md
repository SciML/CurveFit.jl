# Solutions
 
Curve fitting results are returned as `CurveFitSolution` objects. These objects store the
fitted coefficients, residuals, the problem solved and the algorithm used to solve it. It also
stores a return code which holds information about the success of the solver. For more information
on the `ReturnCode` type, see
[SciMLBase.jl](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/#retcodes).

```@docs
CurveFitSolution
```

## StatsAPI functions
`CurveFitSolution` objects can be treated as statistical models and CurveFit
defines various statistics methods for the `CurveFitSolution` type, many of
which extend StatsAPI.jl functions. Import these functions from StatsAPI.jl.

```@docs
StatsAPI.residuals(::CurveFit.CurveFitSolution)
CurveFit.mse
StatsAPI.dof(::CurveFit.CurveFitSolution)
StatsAPI.dof_residual(::CurveFit.CurveFitSolution)
StatsAPI.predict(::CurveFit.CurveFitSolution)
StatsAPI.coef(::CurveFit.CurveFitSolution)
StatsAPI.nobs(::CurveFit.CurveFitSolution)
StatsAPI.fitted(::CurveFit.CurveFitSolution)
StatsAPI.rss(::CurveFit.CurveFitSolution)
CurveFit.isconverged
StatsAPI.vcov(::CurveFit.CurveFitSolution)
StatsAPI.stderror(::CurveFit.CurveFitSolution)
CurveFit.margin_error
StatsAPI.confint(::CurveFit.CurveFitSolution)
```
