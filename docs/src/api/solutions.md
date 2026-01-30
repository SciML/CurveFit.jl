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
which extend StatsAPI.jl functions.

```@docs
residuals
mse
dof
dof_residual
predict
coef
nobs
fitted
rss
isconverged
vcov
stderror
margin_error
confint
```
