module CurveFit

using LinearAlgebra
using Statistics

export curve_fit, linear_fit, poly_fit, exp_fit, log_fit, power_fit
export FitResult, CurveFitOptions, CurveFitException

include("types.jl")
include("linear.jl")
include("polynomial.jl")
include("exponential.jl")
include("logarithmic.jl")
include("power.jl")
include("general.jl")

end