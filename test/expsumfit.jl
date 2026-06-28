using CurveFit
using Test

@testset "ExpSumFit" begin
    x = collect(0.02:0.02:1.5)
    y = @. 5 * exp(0.5 * x) + 4 * exp(-3 * x) + 2 * exp(-2 * x) - 3 * exp(0.15 * x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, ExpSumFitAlgorithm(; n = 4, m = 2, withconst = false))

    @test sol.u.λ ≈ [-3, -2, 0.15, 0.5] rtol = 1.0e-3
    @test sol.u.p ≈ [4, 2, -3, 5] rtol = 1.0e-2

    y = @. -1 + 5 * exp(0.5 * x) + 4 * exp(-3 * x) + 2 * exp(-2 * x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, ExpSumFitAlgorithm(; n = 3, m = 1, withconst = true))

    @test sol.u.λ ≈ [-3, -2, 0.5] rtol = 7.0e-4
    @test sol.u.p ≈ [4, 2, 5] rtol = 3.0e-3
    @test sol.u.k[] ≈ -1 rtol = 2.0e-3

    sol = solve(prob, ExpSumFitAlgorithm(; n = 3, m = 2, withconst = true))
    @test sol.u.λ ≈ [-3, -2, 0.5] rtol = 5.0e-7
    @test sol.u.p ≈ [4, 2, 5] rtol = 9.0e-6
    @test sol.u.k[] ≈ -1 rtol = 2.0e-6
    @test sol(x) ≈ y rtol = 2.0e-6
    @test sol(x[1]) ≈ y[1] rtol = 2.0e-6

    # decay curve
    fs, ts, ω₀, τ = 20.0e3, 0.2, 6283.2, 0.0322
    t = range(0, step = 1 / fs, stop = ts)
    y = @. 1.23 * exp(-t / τ) * cos(ω₀ * t)

    prob = CurveFitProblem(collect(range(0, step = 1 / fs, length = length(y))), y)
    sol = solve(prob, ExpSumFitAlgorithm(; n = 2, m = 2))

    y_fit = sol(t)
    @test isapprox(y, y_fit, rtol = 9.0e-3)
    sol = solve(prob, ExpSumFitAlgorithm(; n = 2, m = 4))
    y_fit = sol(t)
    @test isapprox(y, y_fit, rtol = 6.0e-4)
    sol = solve(prob, ExpSumFitAlgorithm(; n = 2, m = 6))
    y_fit = sol(t)
    @test isapprox(y, y_fit, rtol = 5.0e-5)
end
