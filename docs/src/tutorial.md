# Tutorial

This tutorial covers the supported problem types. This includes
linear and non linear curve fitting along with some special cases 
such as polynomial, exponential, power law and sum of exponentials
curve fitting.

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

## API

See the [API Reference](@ref) page for detailed documentation of all exported functions.
