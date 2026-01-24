# Problems

CurveFit defines two primary problem types: `CurveFitProblem` and `NonlinearCurveFitProblem`. 
These types encapsulate the data to be fitted, along with optional model definitions and
initial parameter guesses. The data and the initial guess are subtypes of `AbstractArray`.

A `CurveFitProblem` represents a general curve fitting problem and requires only input
data `x` and output data `y`. When no model is specified, both `x` and `y` must be
one-dimensional arrays.

`NonlinearCurveFitProblem` is a convenience constructor for defining nonlinear fitting
problems. It returns a `CurveFitProblem` object with a user-provided model function 
and an initial guess `u0` for the model parameters. The model is supplied as a standard 
Julia function and is internally wrapped as a `NonlinearFunction`. For details, see the 
documentation for [`NonlinearFunction`](https://docs.sciml.ai/NonlinearSolve/stable/basics/nonlinear_functions/).
If the output data `y` is not provided, it is treated as a zero vector.

## Constructors 

```@docs
CurveFitProblem
NonlinearCurveFitProblem
```

## Scalar models
For convenience when creating a model CurveFit provides a helper `ScalarModel`
type that allows defining scalar models.

```@docs
ScalarModel
```
