using CurveFit
using CommonSolve: init, solve, solve!
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

@testset "Weighted residual-style problem (y = nothing)" begin
    x = collect(1.0:10.0)
    y = @. 2.0 + 3.0 * x + 0.1 * sin(x)
    sigma = collect(range(0.5, 2.0; length = length(x)))
    u0 = [1.0, 1.0]

    pred(u, x) = @. u[1] + u[2] * x
    sol_pred = solve(NonlinearCurveFitProblem(pred, u0, x, y, sigma))

    # Residual-style problems (y = nothing) must apply sigma the same way as
    # the equivalent prediction-style problem.
    resid_oop(u, x) = y .- pred(u, x)
    sol_oop = solve(NonlinearCurveFitProblem(resid_oop, u0, x, nothing, sigma))
    @test sol_oop.u ≈ sol_pred.u
    @test sol_oop.resid ≈ sol_pred.resid

    resid_iip(resid, u, x) = resid .= y .- pred(u, x)
    sol_iip = solve(NonlinearCurveFitProblem(resid_iip, u0, x, nothing, sigma))
    @test sol_iip.u ≈ sol_pred.u
    @test sol_iip.resid ≈ sol_pred.resid
end
