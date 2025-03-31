# CurveFit

[![Build Status](https://github.com/pjabardo/CurveFit.jl/workflows/CI/badge.svg)](https://github.com/pjabardo/CurveFit.jl/actions)
[![Coverage](https://codecov.io/gh/pjabardo/CurveFit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pjabardo/CurveFit.jl)


A package that implements a few curve fitting functions.

## Linear Least squares

Linear least square is commonly used technique to find approximation to a discrete
set of data. Given the sets of points `x[i]` and `y[i]` and a list of functions `f_i(x)`
the least squares method finds coefficients `a[i]` such that

```
a[1]*f_1(x) + a[2]*f_2(x) + ... + a[n]*f_n(x)
```
minimizes the squares of the errors in relation to `y[i]`.

The basic function is implemented using QR decomposition: `A \ y`:
```
coefs = A \ y
```
where `A[:,i] = f_i(x)`. While usually `x` is a single variable, in general, if several
independent variables are required, the same procedure can be used, something along the line of: 
`A[:,i] = f_i(x1, x2, ..., xn)`.

Several typical cases are possible:
 * `linear_fit(x, y)` finds coefficients `a` and `b` for `y[i] = a + b*x[i]`
 * `power_fit(x, y)` finds coefficients `a` and `b` for `y[i] = a *x[i]^b`
 * `log_fit(x, y)` finds coefficients `a` and `b` for `y[i] = a + b*log(x[i])`
 * `exp_fit(x, y)` finds coefficients `a` and `b` for `y[i] = a*exp(b*x[i])`
 * `expsum_fit(x, y, 2, withconst = true)` finds coefficients `k`, `p`, and `λ` for `y[i] = k + p[1]*exp(λ[1]*x[i]) + p[2]*exp(λ[2]*x[i])`
 * `poly_fit(x, y, n)` finds coefficients `a[k]`  for 
   `y[i] = a[1] + a[2]*x[i] + a[3]*x[i]^2 + a[n+1]*x[i]^n`
 * `linear_king_fit(E, U)`, find coefficients `a` and `b` for `E[i]^2 = a + b*U^0.5`
 * `linear_rational_fit(x, y, p, q)` finds the coefficients for rational polynomials: `y[i] = (a[1] + a[2]*x[i] + ... + a[p+1]*x[i]^p) / (1 + a[p+1+1]*x[i] + ... + a[p+1+q]*x[i]^q)`

## Nonlinear least squares

Sometimes the fitting function is non linear with respect to the  fitting coefficients. In this case, given
an approximation of the coefficients, the fitting function is linearized around this 
approximation and linear least squares is used to calculate a correction to the approximate coefficients. This iteration is repeated until convergence is 
reached. The fitting function has the following form:
```
f(x_1, x_2, x_3, ..., x_n, a_1, a_2, ...,  a_p) = 0
```
where `xi` are the known data points and `ai` are the coefficients that 
should be fitted. 

When the model formula is not linear on the fitting coefficients, a nonlinear algorithm is necessary. This library implements a Newton-type algorithm that doesn't explicitly need derivatives. This is implemented in the function:

`coefs, converged, iter = nonlinear_fit(x, fun, a0, eps=1e-7, maxiter=200)`

In this function, `x` is an array where each column represents a different variable of the data set,
`fun` is a callable that returns the fitting error and should be callable with the following signature:

`residual = fun(x, a)`

where `x` is a vector representing a row of the argument array `x` and `a` is an estimate of the
fitting coefficients which should all be different from zero (to provide a scale). `eps` and `maxiter`
are convergence parameters.

The `nonlinear_fit` function is used to implement the following fitting functions.

 * `king_fit(E, U)` find coefficients `a`, `b` and `n` for `E[i]^2 = a + b*U^n`
 * `rational_fit` Just like `linear_rational_fit` but tries to improve the results using nonlinear least squares (`nonlinear_fit`)

## Example
```julia
using PyPlot
using CurveFit

x = 0.1:0.1:10.0
y =  @. 2.0 / x  * exp(-1.5 / x ) + randn()/100

fun(xy, a) = xy[2] - a[1] / xy[1] * exp( - a[2] / xy[1])
coefs, converged, iter = nonlinear_fit([x y], fun, [2.1 1.4], 1e-7, 200)
y_fit = @. coefs[1] / x * exp(-coefs[2] / x)

plot(x, y, "o", x, y_fit, "r-", linewidth=3)
```

## Generic interface

A generic interface was developed to have a common interface for all curve fitting possibilities and to make it easy to use the results:

`fit = curve_fit(::Type{T}, x, y...)`

where `T` is the curve fitting type.

The following cases are implemented:

 * `curve_fit(LinearFit, x, y)` 
 * `curve_fit(LogFit, x, y)`
 * `curve_fit(PowerFit, x, y)`
 * `curve_fit(ExpFit, x, y)`
 * `curve_fit(Polynomial, x, y, n=1)`
 * `curve_fit(LinearKingFit, E, U)`
 * `curve_fit(KingFit, E, U)`
 * `curve_fit(RationalPoly, x, y, p, q)`

The `curve_fit` generic function returns an object that can be use to compute estimates of the model with `apply_fit`. `call` is overloaded so that the object can be used as a function.



## Example
```julia
using PyPlot
using CurveFit

x = 0.0:0.02:2.0
y0 = @. 1 + x + x*x + randn()/10
fit = curve_fit(Polynomial, x, y0, 2)
y0b = fit.(x) 
plot(x, y0, "o", x, y0b, "r-", linewidth=3)
```


## King's law

In hotwire anemometry, a simple expression for the calibration curve of the probe 
is known as King's law, expressed as:
```
E^2 = A + B*sqrt(U)
```
where `E` is voltage on the anemometer bridge, `U` is the flow velocity.
The coefficients A and B are obtained from a calibration. The function
`linear_king_fit` estimates coefficients `A` and `B`.

A better approximation for the calibration curve is known as modified
King's law:
```
E^2 = A + B*U^n
```
Now, this is a nonlinear curve fit. The linear fit (`linear_king_fit`) is usually
a very good first guess for the coefficients (where `n=0.5`). This curve fit is 
implemented in function `king_fit`.







### Example
```julia
using PyPlot
using CurveFit

U = 1.0:20.0
E = @. sqrt(2 + 1 * U ^ 0.45) + randn()/60
e = range(minimum(E), maximum(E), length=50)

f1 = curve_fit(KingFit, E, U)
U1 = f1.(e)

f2 = curve_fit(Polynomial, E, U, 5)
U2 = f2.(e)

plot(U, E, "o", U1, e, "r-", U2, e, "g-", linewidth=3)
```
