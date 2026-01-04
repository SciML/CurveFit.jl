# Tutorial

This tutorial introduces the basic workflow of CurveFit by walking through
several simple examples. 

## Linear Fitting

Fit a linear function `y = a * x + b`:

```@example linear
using CurveFit

# Generate sample data: y = 2.5 * x + 3.0
x = collect(0:0.1:10)
y = @. 2.5 * x + 3.0

# Create the problem and solve
prob = CurveFitProblem(x, y)
sol = solve(prob, LinearCurveFitAlgorithm())

# Access the coefficients: sol.u = (a, b)
println("Slope (a): ", sol.u[1])
println("Intercept (b): ", sol.u[2])

# Evaluate the solution at a point
println("Prediction at x=5: ", sol(5.0))
```

## Polynomial Fitting

Fit a polynomial of a given degree:

```@example polynomial
using CurveFit

# Generate sample data: y = 1.0 + 2.0*x + 3.0*x^2
x = collect(range(1, stop=10, length=20))
y = @. 1.0 + 2.0 * x + 3.0 * x^2

# Create the problem and solve with degree 2 polynomial
prob = CurveFitProblem(x, y)
sol = solve(prob, PolynomialFitAlgorithm(degree=2))

# Access the coefficients: [c0, c1, c2] for c0 + c1*x + c2*x^2
println("Coefficients: ", sol.u)

# Evaluate the solution at a point
println("Prediction at x=5: ", sol(5.0))
```

## Exponential Fitting

Fit an exponential function `y = b * exp(a * x)`:

```@example exponential
using CurveFit

# Generate sample data: y = 2.0 * exp(0.3 * x)
x = collect(range(0, stop=5, length=20))
y = @. 2.0 * exp(0.3 * x)

prob = CurveFitProblem(x, y)
sol = solve(prob, ExpCurveFitAlgorithm())

# sol.u[1] = a (exponent coefficient)
# sol.u[2] = log(b)
println("Exponent coefficient (a): ", sol.u[1])
println("Scale factor (b): ", exp(sol.u[2]))
```

## Power Law Fitting

Fit a power law `y = b * x^a`:

```@example power
using CurveFit

# Generate sample data: y = 2.0 * x^0.8
x = collect(range(1, stop=10, length=20))
y = @. 2.0 * x^0.8

prob = CurveFitProblem(x, y)
sol = solve(prob, PowerCurveFitAlgorithm())

# sol.u[1] = a (exponent)
# sol.u[2] = log(b)
println("Exponent (a): ", sol.u[1])
println("Scale factor (b): ", exp(sol.u[2]))
```

## Nonlinear Curve Fitting

For arbitrary nonlinear functions, use `NonlinearCurveFitProblem`:

```@example nonlinear
using CurveFit

# Define a nonlinear function: y = a[1] + a[2] * x^a[3]
fn(a, x) = @. a[1] + a[2] * x^a[3]

# True parameters
true_params = [3.0, 2.0, 0.7]

# Generate sample data
x = collect(1.0:0.5:10.0)
y = fn(true_params, x)

# Create problem with initial guess for parameters
u0 = [0.5, 0.5, 0.5]
prob = NonlinearCurveFitProblem(fn, u0, x, y)
sol = solve(prob)

println("Fitted parameters: ", sol.u)
println("Prediction at x=5: ", sol(5.0))
```

## Sum of Exponentials

Fit a sum of exponentials: `y = k + p * exp(λ * t)`:

```@example expsum
using CurveFit

# Generate sample data: y = 2.0 + 3.0*exp(-0.5*t)
t = collect(range(0, stop=10, length=50))
y = @. 2.0 + 3.0 * exp(-0.5 * t)

prob = CurveFitProblem(t, y)
sol = solve(prob, ExpSumFitAlgorithm(n=1, withconst=true))

# Access fitted parameters (returned as arrays for n exponentials)
# k[] extracts the scalar from a 1-element array
println("Constant (k): ", sol.u.k[])
println("Amplitude (p): ", sol.u.p[])
println("Decay rate (λ): ", sol.u.λ[])
```

## Modified King fitting

Fit with the modified king law: `x^2 = a + b * y^n`

```@example king
using CurveFit

# Generate the data: x^2 = a + b * y^n
x = collect(10.0:20.0)
θ_ref = [5.0, 1.3, 2.5]
y = @. ((X^2 - θ_ref[2])/θ_ref[3])^(1/θ_ref[1])

# Works with and without an initial guess
prob = CurveFitProblem(x, y)
sol = solve(prob, ModifiedKingCurveFitAlgorithm())

println("Fitted parameters: ", sol.u)
```

## Rational polynomial fitting

Fit a rational function: `y = p(x)/q(x)`

```@example rational
using CurveFit

# Generate sample data: y = (1 + 2*x) / (1 + 0.5*x - 0.1*x^2)
x = collect(range(0, stop=5, length=30))
y = @. (1.0 + 2.0*x) / (1.0 + 0.5*x - 0.1*x^2)

# Create the problem and solve with numerator degree 1, denominator degree 2
prob = CurveFitProblem(x, y)
alg = RationalPolynomialFitAlgorithm(num_degree=1, den_degree=2)
sol = solve(prob, alg)

# Access fitted coefficients
# Numerator: p0 + p1*x
println("Numerator coefficients: ", sol.u[1:2])
# Denominator: 1 + q1*x + q2*x^2
println("Denominator coefficients: ", vcat(1.0, sol.u[3:4]))

# Evaluate the solution at a point
println("Prediction at x=2: ", sol(2.0))
```

## Using StatsAPI

```@example stats
using CurveFit
using NonlinearSolve

# Generate sample data: y = 3 + 2*x + x^1
X = collect(1.0:10.0)
θ_ref = [3.0, 2.0, 1.0]

function f(θ, X)
    @. θ[1] + θ[2]*X + X^θ[3]
end

Y = f(θ_ref, X)

# Wrap as NonlinearFunction and define problem
nonf = NonlinearFunction(f)
prob = NonlinearCurveFitProblem(nonf, [0.5, 0.5, 0.5], X, Y)

# Solve using Levenberg-Marquardt
sol = solve(prob, LevenbergMarquardt())

# --- StatsAPI functions ---

# Access fitted coefficients, [θ1, θ2, θ3]
println("Coefficients: ", coef(sol))          

# Residuals
println("Residuals: ", residuals(sol))

# Predicted values at a single point
println("Prediction at X=5: ", predict(sol, 5.0))

# Predicted values at all X
println("Fitted values: ", fitted(sol))

# Number of observations and degrees of freedom
println("Number of observations: ", nobs(sol))
println("Degrees of freedom: ", dof(sol))
println("Residual degrees of freedom: ", dof_residual(sol))

# Sum of squared residuals and mean squared error
println("RSS: ", rss(sol))
println("MSE: ", mse(sol))

# Check convergence
println("Converged: ", isconverged(sol))

# Covariance matrix and standard errors
vcov(sol)
stderror(sol)

# Confidence intervals
confint(sol)
```

## API

See the [API Reference](@ref) page for detailed documentation of all exported functions.
