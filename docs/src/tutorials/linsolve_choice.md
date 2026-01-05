# Advanced usage

This tutorial covers advanced options for improving curve fitting performance
by controlling the linear solvers used internally by nonlinear algorithms.
It is intended for users who want more control over numerical stability,
performance, or problem structure. NonlinearSolve.jl provides the user with the
option to choose the linear solver they want for the Levenber-Marquardt
algorithm.

## Why does linear solver choice matter

The Levenberg–Marquardt algorithm work by repeatedly solving linear
systems involving the Jacobian matrix. At each iteration, systems of
the form:

```math
    (J^T \cdot J + \lamda \cdot I)\delta = J^T \cdot r
```

must be solved. `J` is the Jacobian matrix and `r` the residual
vector.

The numerical properties of these systems depend strongly on the structure
and conditioning of the Jacobian. While NonlinearSolve.jl selects a reasonable 
default, explicitly choosing a linear solver can significantly improve robustness 
or performance in some cases.

Inside the `LevenbergMarquardt` constructor the keyword argument `linsolve` specifies
the linear solver used in the algorithm. By default it is set to nothing and in that 
case NonlinearSolve.jl will choose a linear solver on its own. Users can pass their
desired linear solver to the constructor in which case that will be the one used.
Below is a simple example of how to pass a specific linear solver.    

## Example

```@example qr
using CurveFit
using NonlinearSolve

X = collect(1.0:10.0)
θ_true = [3.0, 1.5, 0.5]

function f(θ, X)
    return @. θ[1] + θ[2] * exp(θ[3] * X)
end
Y = f(θ_true, X)

nonf = NonlinearFunction(f)
alg = LevenbergMarquardt(linsolve = QRFactorization())
 
prob = NonlinearCurveFitProblem(nonf, [2.1, 3.3, 0.1], X, Y)
sol = solve(prob, alg)
```

## Choosing the right solver

As a general guideline:

- QR-based factorization is recommended for ill-conditioned or
rank-deficient problems. It is robust and numerically-stable.
Slower than cholesky factorization.

- Cholesky-based factorization is efficient for well-conditioned
problems and for symmetric positive definite matrices. Faster
than qr factorization, but can be unstable.

- LU-based factorization falls somewhere in between QR and Cholesky
factorization. Can be used, although the other two are usually
better.

- Default solver is usually a decent choice, unless some structural
information about the problem is known which can be exploited.

## When to customize

Custom solver selection is most useful when:

- Fits fail to converge due to numerical instability

- Jacobians are poorly conditioned

- Problem structure (dense, sparse, symmetric) is known in advance

- Performance becomes critical for large problems

In most cases, the default choice is more than enough, but advanced
users can benefit from this additional level of control.


For more detail see:
- [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/)
- [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve/stable/)
