# Solvers

CurveFit provides built-in solvers for linear curve fitting problems. Nonlinear problems
are delegated to `NonlinearSolve.jl`. In addition, CurveFit includes specialized algorithms
for selected nonstandard models.

## Linear fitting

Linear curve fitting in CurveFit solves problems of the general form $ f_y(y) = a \cdot f_x(x) + b $
where `x` and `y` are the data points being fitted and `a` and `b` are the fit parameters.

Linear fits do not require an initial guess. The user must explicitly select a linear
algorithm, as no default is assumed.

Linear fitting algorithms are represented by the `LinearCurveFitAlgorithm` type, which
encapsulates transformations applied to the input and output data. The fields `xfun` and
`yfun` define transformations applied to `x` and `y`, respectively, while
`yfun_inverse` maps fitted parameters back to the original data space.

The default constructor corresponds to the standard linear model $ y = a \cdot x + b $
where both transformations are identity functions. Users may supply custom transformations
to define alternative linear relationships. The inverse transformation is computed using
`InverseFunctions.jl`.

These are the convenience constructors exported by CurveFit used for defining some linear
curve fitting algorithms:

```@docs
LinearCurveFitAlgorithm
LogCurveFitAlgorithm
PowerCurveFitAlgorithm
ExpCurveFitAlgorithm
KingCurveFitAlgorithm
```

## Nonlinear fitting

Nonlinear curve fitting problems are solved through NonlinearSolve.jl. The user defines 
a nonlinear problem using `NonlinearCurveFitProblem`, supplying a model function and an 
initial guess for the parameters.

During initialization, CurveFit detects nonlinear problems by checking if a model function
is provided and internally constructs a `NonlinearLeastSquaresProblem`, which is then passed 
to NonlinearSolve.jl for solution. For details, see the documentation for [NonlinearLeastSquaresProblem](https://docs.sciml.ai/NonlinearSolve/stable/basics/nonlinear_problem/).

By default, NonlinearSolve.jl automatically selects an appropriate nonlinear least-squares
algorithm. Advanced users may explicitly specify a solver, such as `LevenbergMarquardt` or
`GaussNewton`. Documentation for available solvers can be found [here](https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_least_squares_solvers/).

CurveFit provides a convenience constructor for choosing the linear solver passed to
the Levenberg-Marquardt algorithm from NonlinearSolve.jl.

```@docs
LM_linsolve
```

For more information see [Tutorials](@ref).

## Special functions

CurveFit provides specialized algorithms for fitting selected classes of nonlinear models
with known structure. These algorithms exploit problem-specific properties to improve
robustness and performance compared to fully generic nonlinear least-squares approaches.

### Modified King fitting

The Modified King fit algorithm is designed for the model function of the form:

```math
x^2 = a + b \cdot y^n
```

Unlike the linear King fit where the exponent is a constant (1/2), `n` in the exponent is 
a coefficient and thus the problem becomes nonlinear. Since it is nonlinear solving it is
done via NonlinearSolve.jl. This Modified King fit has a defined model and if the user tries 
to solve a problem with a specified model function using this algorithm an error will be thrown.
However, the user can pass an initial guess for the coefficients `a`, `b` and `n`. In case an
initial guess is not provided CurveFit will obatin it internally by using the related linear
King fit.

```@docs
ModifiedKingCurveFitAlgorithm
```

### Sum of exponentials fitting

The sum of exponentials fitting algorithm is defined for models of the form:

```math
y = k + \sum_{i=1}^{n} p_i \exp(\lambda_i x)
```

The number of exponential terms is specified by the user. This method exploits the algebraic 
structure of exponential sums to avoid solving a nonlinear problem. By computing discrete 
cumulative integrals of the data, the problem is transformed into a linear recurrence 
whose coefficients can be estimated using linear least squares. The exponential rates 
\(\lambda_i\) are then recovered as eigenvalues of a matrix derived from this recurrence.

Once the rates are known, the amplitudes \(p_i\) and the optional constant term \(k\)
are obtained by solving a linear least squares problem with fixed exponential basis
functions.

For more details see [here](https://math.stackexchange.com/questions/1428566/fit-sum-of-exponentials).

```@docs
ExpSumFitAlgorithm
```

### Polynomial fitting

#### Standard

Polynomial fitting solves problems of the form

```math
y = \sum_{k=0}^{n} a_k x^k
```

The polynomial degree is specified by the user. The problem is formulated as 
a `LinearProblem` using a Vandermonde matrix constructed from the input data
The polynomial coefficients are obtained by solving this linear system with
a user-selectable linear solver. 

Polynomial fitting is only applicable to linear curve fitting problems and does
not support user-provided initial guesses. It is recommended to use numerically
stable solvers such as QR-based factorizations. CurveFit allows the linear solver
to be customized via the `linsolve_algorithm` field of `PolynomialFitAlgorithm`.

```@docs
PolynomialFitAlgorithm
```

#### Rational

Rational polynomial fitting solves problems of the form

```math
y = \frac{p(x)}{q(x)}
```

\(p(x)\) and \(q(x)\) are polynomials of user-specified degrees. The constant term
of the denominator is assumed to be 1. Rational polynomial fitting implements both
linear and nonlinear methods. `RationalPolynomialFitAlgorithm` has three fields, 
the first two are degrees of the numerator and denominator and third is the `alg`
field which contains the algorithm used for solving the problem.

In case the user chooses a linear fit algorithm initial guesses are not supported.
The problem is transformed into $y \cdot q(x) = p(x)$. This linearized problem is
then solved by creating linear fit cache with `LinearProblem` which is then handled
by the internal solver for linear rational fit of CurveFit. For more detail on
`Linear problem` [see](https://docs.sciml.ai/LinearSolve/stable/basics/LinearProblem/).

If the user chooses a nonlinear algorithm, then the problem is solved via NonlinearSolve.jl.
In this case the user can provide an initial guess. If it is not provided, then CurveFit
will obtain one from the linear rational fit.

```@docs
RationalPolynomialFitAlgorithm
```

