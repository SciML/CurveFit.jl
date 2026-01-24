# StatsAPI interface usage

This tutorial goes over basic functionality of the StatsAPI.jl interface
as implemented by CurveFit. Solvers find coefficients of the model so that
it fits the data as best as it can. Statistical tools allow user to assess
how good and reliable their fitting result is. 

After solving a curve fitting problem (does not have to be a nonlinear problem),
you can use statistical functions on the `CurveFitSolution` object. See
[StatsAPI functions](@ref) to view all the implemented functions.

## Examples

```@example stats
using CurveFit

x = collect(1.0:10.0)
θ_true = [3.0, 2.0, 1.5]

f(θ, x) = @. θ[1] + θ[2] * x + x^θ[3]
y = f(θ_true, x)

prob = NonlinearCurveFitProblem(f, [1.0, 1.0, 1.0], x, y)
sol = solve(prob)

# Get the confidence intervals
confint(sol)
```

## Basic quantities

- [`residuals()`](@ref) measure the difference between the fitted model and the data.
- The residual sum of squares ([`rss()`](@ref)) and mean squared error
  ([`mse()`](@ref)) summarize the overall fit quality.
- The number of observations ([`nobs()`](@ref)) corresponds to the size of the
  data set used, i.e the number of data points.
- [`predict()`](@ref) gives a prediction using the fitted coefficients and new
  data. If only the solution object is passed, original data will be used in
  calculation.
- [`isconverged()`](@ref) checks if the solver was successful in solving the
  problem.

See [StatsAPI functions](@ref) to view all the implemented functions.

## Parameter uncertainty

CurveFit exposes parameter uncertainty through standard errors and covariance
matrices. The covariance matrix estimates the joint uncertainty of the fitted
coefficients under standard least squares assumptions. The diagonal entries
correspond to the variance of each coefficient estimate, while the off-diagonal
entries quantify correlations between coefficients. Standard errors are obtained
as the square roots of the diagonal variances.

## Confidence intervals

Point estimates alone do not convey how uncertain a fitted coefficient is.
Confidence intervals provide a range of values that are statistically
consistent with the observed data under standard modelling assumptions.

Confidence intervals can be computed with [`confint()`](@ref).

### Interpreting confidence intervals

A 95% confidence interval means that, under repeated experiments with the same
data-generating process, approximately 95% of such intervals would contain the
true coefficient value.

Narrow intervals indicate well-determined coefficients, while wide intervals
suggest that the data provide limited information about a coefficient.

Confidence intervals are commonly used to assess the reliability of the
fitted coefficients.
