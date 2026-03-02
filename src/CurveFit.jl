module CurveFit

using CommonSolve: CommonSolve, init, solve!, solve
using ConcreteStructs: @concrete
using InverseFunctions: inverse
using Markdown: @doc_str
using PrecompileTools: @compile_workload, @setup_workload
using Setfield: @set, @set!

using RecursiveArrayTools: NamedArrayPartition
using FastRationals: FastRational
using LinearAlgebra: LinearAlgebra, eigvals!, diagm, qr!, lu!, ldiv!
using NonlinearSolveFirstOrder: NonlinearSolveFirstOrder
using NonlinearSolveBase: NonlinearSolveBase
using SciMLBase: SciMLBase, AbstractNonlinearAlgorithm, AbstractLinearAlgorithm, ReturnCode,
    NonlinearFunction, LinearProblem, NonlinearLeastSquaresProblem, reinit!
using DifferentiationInterface: DifferentiationInterface
using ADTypes: AutoForwardDiff
using Distributions: TDist, quantile
using StatsAPI: StatsAPI, coef, residuals, predict, fitted, nobs, dof, dof_residual,
    rss, vcov, stderror, confint

# Abstract base class for fitting data
abstract type AbstractApproxFit end

# Abstract class for least squares fitting of data
abstract type AbstractLeastSquares <: AbstractApproxFit end

include("common_interface.jl")

include("linfit.jl")
include("rationalfit.jl")
include("nonlinfit.jl")
include("king.jl")
include("expsumfit.jl")
include("stats.jl")

# Exported functions
export CurveFitProblem, NonlinearCurveFitProblem, ScalarModel

export LinearCurveFitAlgorithm, LogCurveFitAlgorithm, PowerCurveFitAlgorithm,
    ExpCurveFitAlgorithm, PolynomialFitAlgorithm
export RationalPolynomialFitAlgorithm
export KingCurveFitAlgorithm, ModifiedKingCurveFitAlgorithm
export ExpSumFitAlgorithm

export CurveFitSolution

export solve, solve!, init

export coef, residuals, predict, fitted, nobs, dof, dof_residual, rss, mse, vcov, stderror, margin_error, confint
export isconverged

@setup_workload begin
    x_range = 1.0:10.0
    x_array = collect(1.0:10.0)
    fn(a, x) = @. a[1] + a[2] * x^a[3]
    a0 = [3.0, 2.0, 0.7]
    y = fn(a0, x_range)

    @compile_workload begin
        prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], x_range, y)
        solve(prob)

        prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], x_array, y)
        solve(prob)
    end
end

end
