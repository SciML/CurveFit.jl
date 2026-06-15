using CurveFit
using Test
using SciMLBase

@testset "Nonlinear Weighted Least Squares" begin
    fn(a, x) = @. a[1] + a[2] * x
    x = collect(1.0:10.0)
    a0 = [1.0, 2.0]
    y = fn(a0, x)

    # Add a large outlier
    y[5] += 20.0

    # Fit without sigma - outlier should affect result
    prob_no_weight = NonlinearCurveFitProblem(fn, [0.5, 0.5], x, y)
    sol_no_weight = solve(prob_no_weight)

    # Fit with a high sigma on the outlier
    sigma = ones(length(y))
    sigma[5] = 100.0
    prob_weighted = NonlinearCurveFitProblem(fn, [0.5, 0.5], x, y, sigma)
    sol_weighted = solve(prob_weighted)

    # Weighted fit should (hopefully) be closer to the true parameters
    err_no_weight = maximum(abs.(sol_no_weight.u .- a0))
    err_weighted = maximum(abs.(sol_weighted.u .- a0))

    @test err_weighted < err_no_weight
    @test SciMLBase.successful_retcode(sol_weighted)
end
