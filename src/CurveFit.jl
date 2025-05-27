module CurveFit

using CommonSolve: CommonSolve, init, solve!, solve
using ConcreteStructs: @concrete
using InverseFunctions: inverse
using Markdown: @doc_str

using LinearSolve: LinearSolve
using SciMLBase: AbstractLinearAlgorithm

using SciMLBase
using NonlinearSolve
using LinearAlgebra

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

# Exported functions
export LinearCurveFitProblem, LogCurveFitProblem, PowerCurveFitProblem, ExpCurveFitProblem
export PolynomialFitAlgorithm
export LinearCurveFitSolution

export solve, solve!, init

end
