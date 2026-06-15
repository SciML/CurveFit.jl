using CurveFit
using Test
using SciMLBase

@testset "Linear Weighted Least Squares" begin
    x = 1:10
    a0 = [2.0, 1.0]
    fn(x) = a0[1] * x + a0[2]
    y = fn.(x)

    # Add a large outlier
    y[5] += 20.0

    # Fit without sigma - outlier should affect result
    prob = CurveFitProblem(x, y)
    sol_unweighted = solve(prob, LinearCurveFitAlgorithm())

    # Fit with a high sigma on the outlier
    sigma = ones(length(y))
    sigma[5] = 100.0
    prob = CurveFitProblem(x, y; sigma)
    sol_weighted = solve(prob, LinearCurveFitAlgorithm())

    # Weighted fit should (hopefully) be closer to the true parameters
    err_no_weight = maximum(abs.(sol_unweighted.u .- a0))
    err_weighted = maximum(abs.(sol_weighted.u .- a0))

    @test err_weighted < err_no_weight
    @test SciMLBase.successful_retcode(sol_weighted)
end
