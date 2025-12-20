@testitem "Documentation Examples: Linear Fit" begin
    # Example from docs/src/index.md
    x = collect(0:0.1:10)
    y = @. 2.5 * x + 3.0

    prob = CurveFitProblem(x, y)
    sol = solve(prob, LinearCurveFitAlgorithm())

    @test sol.u[1] ≈ 2.5
    @test sol.u[2] ≈ 3.0
    @test sol(5.0) ≈ 2.5 * 5.0 + 3.0
end

@testitem "Documentation Examples: Polynomial Fit" begin
    # Example from docs/src/index.md
    x = collect(range(1, stop = 10, length = 20))
    y = @. 1.0 + 2.0 * x + 3.0 * x^2

    prob = CurveFitProblem(x, y)
    sol = solve(prob, PolynomialFitAlgorithm(degree = 2))

    @test sol.u[1] ≈ 1.0
    @test sol.u[2] ≈ 2.0
    @test sol.u[3] ≈ 3.0
    @test sol(5.0) ≈ 1.0 + 2.0 * 5.0 + 3.0 * 5.0^2
end

@testitem "Documentation Examples: Exponential Fit" begin
    # Example from docs/src/index.md
    x = collect(range(0, stop = 5, length = 20))
    y = @. 2.0 * exp(0.3 * x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, ExpCurveFitAlgorithm())

    @test sol.u[1] ≈ 0.3
    @test exp(sol.u[2]) ≈ 2.0
end

@testitem "Documentation Examples: Power Law Fit" begin
    # Example from docs/src/index.md
    x = collect(range(1, stop = 10, length = 20))
    y = @. 2.0 * x^0.8

    prob = CurveFitProblem(x, y)
    sol = solve(prob, PowerCurveFitAlgorithm())

    @test sol.u[1] ≈ 0.8
    @test exp(sol.u[2]) ≈ 2.0
end

@testitem "Documentation Examples: Nonlinear Fit" begin
    # Example from docs/src/index.md
    fn(a, x) = @. a[1] + a[2] * x^a[3]

    true_params = [3.0, 2.0, 0.7]
    x = collect(1.0:0.5:10.0)
    y = fn(true_params, x)

    u0 = [0.5, 0.5, 0.5]
    prob = NonlinearCurveFitProblem(fn, u0, x, y)
    sol = solve(prob)

    @test sol.u ≈ true_params
end

@testitem "Documentation Examples: ExpSum Fit" begin
    # Example from docs/src/index.md
    t = collect(range(0, stop = 10, length = 50))
    y = @. 2.0 + 3.0 * exp(-0.5 * t)

    prob = CurveFitProblem(t, y)
    sol = solve(prob, ExpSumFitAlgorithm(n = 1, withconst = true))

    @test sol.u.k[] ≈ 2.0 rtol = 1e-2
    @test sol.u.p[] ≈ 3.0 rtol = 1e-2
    @test sol.u.λ[] ≈ -0.5 rtol = 1e-2
end
