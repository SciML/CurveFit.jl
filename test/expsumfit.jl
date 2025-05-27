@testitem "ExpSumFit" begin
    x = collect(0.02:0.02:1.5)
    y = @. 5 * exp(0.5 * x) + 4 * exp(-3 * x) + 2 * exp(-2 * x) - 3 * exp(0.15 * x)
    sol = expsum_fit(x, y, 4, withconst = false)
    @test isapprox(sol.λ, [-3, -2, 0.15, 0.5], rtol = 1e-3)
    @test isapprox(sol.p, [4, 2, -3, 5], rtol = 1e-2)

    y = @. -1 + 5 * exp(0.5 * x) + 4 * exp(-3 * x) + 2 * exp(-2 * x)
    sol = expsum_fit(x, y, 3, withconst = true)
    @test isapprox(sol.λ, [-3, -2, 0.5], rtol = 7e-4)
    @test isapprox(sol.p, [4, 2, 5], rtol = 3e-3)
    @test isapprox(sol.k, -1, rtol = 2e-3)

    init = expsum_init(x, 3, m = 2, withconst = true)
    sol = expsum_fit(x, y, 3, m = 2, withconst = true, init = init)
    @test isapprox(sol.λ, [-3, -2, 0.5], rtol = 5e-7)
    @test isapprox(sol.p, [4, 2, 5], rtol = 9e-6)
    @test isapprox(sol.k, -1, rtol = 2e-6)

    # decay curve
    fs, ts, ω₀, τ = 20e3, 0.2, 6283.2, 0.0322
    t = range(0, step = 1 / fs, stop = ts)
    y = @. 1.23 * exp(-t / τ) * cos(ω₀ * t)
    sol = expsum_fit(1 / fs, y, 2, m = 2)
    y_fit = sol(t)
    @test isapprox(y, y_fit, rtol = 9e-3)
    sol = expsum_fit(1 / fs, y, 2, m = 4)
    y_fit = sol(t)
    @test isapprox(y, y_fit, rtol = 6e-4)
    sol = expsum_fit(1 / fs, y, 2, m = 6)
    y_fit = sol(t)
    @test isapprox(y, y_fit, rtol = 5e-5)

    ## Integration rules
    # Trapezoidal rule
    @test CurveFit.calc_integral_rules(1, m = 1) == [1 // 2 1 // 2]
    # Simpson's first (1/3) rule
    @test CurveFit.calc_integral_rules(1, m = 2) == [1 // 3 4 // 3 1 // 3]
    @test CurveFit.calc_integral_rules(2, m = 2) == [2 // 3 4 // 3 0 // 1]
    @test CurveFit.calc_integral_rules(3, m = 2) == [3 // 5 4 // 5 -1 // 15]
    # Simpson's second (3/8) rule
    @test CurveFit.calc_integral_rules(1, m = 3) == [3 // 8 9 // 8 9 // 8 3 // 8]
    @test CurveFit.calc_integral_rules(2, m = 3) == [39 // 40 27 // 10 27 // 40 3 // 20]

    ## Cumulative integrals
    x = collect(0:0.4:10)
    y = @. 1 + sin(x)
    cumints_analytic = @. [(x + 1 - cos(x)) (1 / 2 * x^2 + x - sin(x)) (1 / 6 * x^3 +
                                                                        x^2 / 2 +
                                                                        cos(x) - 1) (1 /
                                                                                     24 *
                                                                                     x^4 +
                                                                                     x^3 /
                                                                                     6 +
                                                                                     sin(x) -
                                                                                     x)]
    p1 = expsum_init(x, 4, m = 1)
    CurveFit.cumints!(p1, x, y)
    @test isapprox(cumints_analytic[:, 1], p1.Y[:, 1], rtol = 3e-3)
    @test isapprox(cumints_analytic[:, 2], p1.Y[:, 2], rtol = 3e-3)
    @test isapprox(cumints_analytic[:, 3], p1.Y[:, 3], rtol = 4e-3)
    @test isapprox(cumints_analytic[:, 4], p1.Y[:, 4], rtol = 4e-3)
    p2 = expsum_init(x, 4, m = 2)
    CurveFit.cumints!(p2, x, y)
    @test isapprox(cumints_analytic[1:2:end, 1], p2.Y[:, 1], rtol = 3e-5)
    @test isapprox(cumints_analytic[1:2:end, 2], p2.Y[:, 2], rtol = 4e-5)
    @test isapprox(cumints_analytic[1:2:end, 3], p2.Y[:, 3], rtol = 5e-5)
    @test isapprox(cumints_analytic[1:2:end, 4], p2.Y[:, 4], rtol = 6e-5)

    x = collect(range(0, stop = 2, length = 13))
    y = @. exp(x)
    cumints_analytic = @. [(exp(x) - 1) (exp(x) - 1 - x) (exp(x) - 1 - x - x^2 / 2) (exp(x) -
                                                                                     1 - x -
                                                                                     x^2 /
                                                                                     2 -
                                                                                     x^3 /
                                                                                     6)]
    p2 = expsum_init(x, 4, m = 2)
    CurveFit.cumints!(p2, x, y)
    @test isapprox(cumints_analytic[1:2:end, 1], p2.Y[:, 1], rtol = 5e-6)
    @test isapprox(cumints_analytic[1:2:end, 2], p2.Y[:, 2], rtol = 3e-5)
    @test isapprox(cumints_analytic[1:2:end, 3], p2.Y[:, 3], rtol = 3e-5)
    @test isapprox(cumints_analytic[1:2:end, 4], p2.Y[:, 4], rtol = 4e-5)
end
