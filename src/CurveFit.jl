module CurveFit

using CommonSolve: CommonSolve, init, solve!, solve
using ConcreteStructs: @concrete
using InverseFunctions: inverse
using Markdown: @doc_str
using Setfield: @set!

using RecursiveArrayTools: NamedArrayPartition
using FastRationals: FastRational
using LinearAlgebra: LinearAlgebra, eigvals!, diagm, qr!, lu!, ldiv!
using LinearSolve: LinearSolve
using NonlinearSolve: NonlinearSolve
using SciMLBase: SciMLBase, AbstractNonlinearAlgorithm, AbstractLinearAlgorithm, ReturnCode,
                 NonlinearFunction, LinearProblem, NonlinearLeastSquaresProblem
using DifferentiationInterface: DifferentiationInterface
using ADTypes: AutoForwardDiff
using Distributions: TDist, quantile
import StatsAPI: coef, residuals, predict, fitted, nobs, dof, dof_residual, rss, vcov, stderror,
                 confint

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
export CurveFitProblem, NonlinearCurveFitProblem

export LinearCurveFitAlgorithm, LogCurveFitAlgorithm, PowerCurveFitAlgorithm,
       ExpCurveFitAlgorithm, PolynomialFitAlgorithm
export RationalPolynomialFitAlgorithm
export KingCurveFitAlgorithm, ModifiedKingCurveFitAlgorithm
export ExpSumFitAlgorithm
export LM_linsolve

export CurveFitSolution

export solve, solve!, init

export coef, residuals, predict, fitted, nobs, dof, dof_residual, rss, mse, vcov, stderror, confint
export isconverged

end