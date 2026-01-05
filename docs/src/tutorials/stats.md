# StatsAPI interface usage

This tutorial goes over basic functionality of the StatsAPI.jl interface
as implemeted by CurveFit. Solvers find coefficients of the model so that
it fits the data as best as it can. Statistical tools allow user to assess
how good and reliable their fitting result is. 

After solving a curve fitting problem (does not have to be a nonlinear problem),
you can use statistical functions on the `CurveFitSolution` object. 

## Setup

```@example setup
using CurveFit
using NonlinearSolve

x = collect(1.0:10.0)
θ_true = [3.0, 2.0, 1.5]

f(θ, x) = @. θ[1] + θ[2] * x + x^θ[3]
y = f(θ_true, x)

nonf = NonlinearFunction(f)
prob = NonlinearCurveFitProblem(nonf, [1.0, 1.0, 1.0], x, y)

sol = solve(prob, LevenbergMarquardt())

residuals(sol)
 
rss(sol)
mse(sol)

nobs(sol)
 
dof(sol)
dof_residual(sol)

predict(sol)
isconverged(sol)

stderror(sol)
vcov(sol)

confint(sol)
```

## Basic quantities

The residuals measure the difference between the fitted model and the data.

The residual sum of squares (rss) and mean squared error (mse) summarize the 
overall fit quality.

The number of observations corresponds to the size of the data set used, i.e
the number of data points.

In CurveFit, `dof` returns the number of the coefficients used in the model,
while `dof_residual` returns the difference between the number of observations
and degrees of freedom (`dof`). These quantities are used internally in
other calculations.

The `predict` gives a prediction using the fitted coefficients and new data.
If only the solution object is passed, original data will be used in calculation.
`isconverged` checks if the solver was successful in solving the problem.


## Parameter uncertainty

CurveFit exposes parameter uncertainty through standard errors and covariance
matrices. The covariance matrix estimates the joint uncertainty of the fitted 
coefficients under standard least squares assumptions. The diagonal entries 
correspond to the variance of each coefficient estimate, while the off-diagonal
entries quantify correlations between coefficient. Standard errors are obtained 
as the square roots of the diagonal variances.

## Confidence intervals

Point estimates alone do not convey how uncertain a fitted coefficient is.
Confidence intervals provide a range of values that are statistically
consistent with the observed data under standard modeling assumptions.

Confidence intervals are computed using the estimated standard error which 
is calculated as the square root of the diagonal entries of the covariance
matrix. The default confidence level is 95% and the associated parameter α
in that case is 0.05. Critial values of the normal distributions are 
calculated in correspondence with α. Upper and lower bounds of the confidence
intervals are then calculated as the product of the standard errors and the
critical values ± the fitted coefficient.

### Interpreting confidence intervals

A 95% confidence interval means that, under repeated experiments with the same
data-generating process, approximately 95% of such intervals would contain the
true coefficient value.

Narrow intervals indicate well-determined cefficients, while wide intervals
suggest that the data provide limited information about a coefficient.

Confidence intervals are commonly used to assess the reliability of the
fitted coefficients.
