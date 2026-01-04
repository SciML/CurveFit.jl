@testitem "ExpSumFit: Integration rules" begin
    # Trapezoidal rule
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 1:1; m = 1) == [1 // 2 1 // 2]

    # Simpson's first (1/3) rule
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 1:1; m = 2) ==
        [1 // 3 4 // 3 1 // 3]
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 2:2; m = 2) ==
        [2 // 3 4 // 3 0 // 1]
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 3:3; m = 2) ==
        [3 // 5 4 // 5 -1 // 15]

    # Simpson's second (3/8) rule
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 1:1; m = 3) ==
        [3 // 8 9 // 8 9 // 8 3 // 8]
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 2:2; m = 3) ==
        [39 // 40 27 // 10 27 // 40 3 // 20]
end

@testitem "ExpSumFit: Cumulative integrals" begin
    n = 4

    x = collect(0:0.4:10)
    y = @. 1 + sin(x)
    cumints_analytic = @. [(x + 1 - cos(x)) (1 / 2 * x^2 + x - sin(x)) (
        1 / 6 * x^3 +
            x^2 / 2 +
            cos(x) - 1
    ) (
        1 /
            24 *
            x^4 +
            x^3 /
            6 +
            sin(x) -
            x
    )]

    @testset "m = 1 & n = 4" begin
        coeffs = CurveFit.__calc_integral_rules(Float64, 1:4; m = 1)

        len = length(x)
        nY, mY = 1 + (len - 1) ÷ 1, 2 * n + 1

        Y = Matrix{Float64}(undef, nY, mY)
        S = Matrix{Float64}(undef, nY, 4 - 1)

        CurveFit.__cumulative_integrals!(Y, S, coeffs, x, y, 4, 1)

        @test cumints_analytic[:, 1] ≈ Y[:, 1] rtol = 3.0e-3
        @test cumints_analytic[:, 2] ≈ Y[:, 2] rtol = 3.0e-3
        @test cumints_analytic[:, 3] ≈ Y[:, 3] rtol = 4.0e-3
        @test cumints_analytic[:, 4] ≈ Y[:, 4] rtol = 4.0e-3
    end

    @testset "m = 2 & n = 4" begin
        coeffs = CurveFit.__calc_integral_rules(Float64, 1:4; m = 2)

        len = length(x)
        nY, mY = 1 + (len - 1) ÷ 2, 2 * n + 1

        Y = Matrix{Float64}(undef, nY, mY)
        S = Matrix{Float64}(undef, nY, 4 - 1)

        CurveFit.__cumulative_integrals!(Y, S, coeffs, x, y, 4, 2)

        @test cumints_analytic[1:2:end, 1] ≈ Y[:, 1] rtol = 3.0e-5
        @test cumints_analytic[1:2:end, 2] ≈ Y[:, 2] rtol = 4.0e-5
        @test cumints_analytic[1:2:end, 3] ≈ Y[:, 3] rtol = 5.0e-5
        @test cumints_analytic[1:2:end, 4] ≈ Y[:, 4] rtol = 6.0e-5
    end

    x = collect(range(0, stop = 2, length = 13))
    y = @. exp(x)
    cumints_analytic = @. [(exp(x) - 1) (exp(x) - 1 - x) (exp(x) - 1 - x - x^2 / 2) (
        exp(x) -
            1 - x -
            x^2 /
            2 -
            x^3 /
            6
    )]

    @testset "m = 2 & n = 4" begin
        coeffs = CurveFit.__calc_integral_rules(Float64, 1:4; m = 2)

        len = length(x)
        nY, mY = 1 + (len - 1) ÷ 2, 2 * n + 1

        Y = Matrix{Float64}(undef, nY, mY)
        S = Matrix{Float64}(undef, nY, 4 - 1)

        CurveFit.__cumulative_integrals!(Y, S, coeffs, x, y, 4, 2)

        @test cumints_analytic[1:2:end, 1] ≈ Y[:, 1] rtol = 5.0e-6
        @test cumints_analytic[1:2:end, 2] ≈ Y[:, 2] rtol = 3.0e-5
        @test cumints_analytic[1:2:end, 3] ≈ Y[:, 3] rtol = 3.0e-5
        @test cumints_analytic[1:2:end, 4] ≈ Y[:, 4] rtol = 4.0e-5
    end
end

@testitem "ExpSumFit" begin
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
