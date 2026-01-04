@testitem "Linear Fit" begin
    x = range(1, stop = 10, length = 10)

    fn(x) = 1.0 + 2.0 * x
    y = fn.(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, LinearCurveFitAlgorithm())

    @test sol.u[1] ≈ 2.0
    @test sol.u[2] ≈ 1.0

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(val)
    end
end

@testitem "Log Fit" begin
    x = range(1, stop = 10, length = 10)

    fn(x) = 1.0 + 2.0 * log(x)
    y = fn.(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, LogCurveFitAlgorithm())

    @test sol.u[1] ≈ 2.0
    @test sol.u[2] ≈ 1.0

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(val)
    end
end

@testitem "Power Fit" begin
    x = range(1, stop = 10, length = 10)

    fn(x) = 2.0 * x .^ 0.8
    y = fn(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, PowerCurveFitAlgorithm())

    @test sol.u[1] ≈ 0.8
    @test sol.u[2] ≈ log(2.0)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(val)
    end
end

@testitem "Exp Fit" begin
    x = range(1, stop = 10, length = 10)

    fn(x) = 2.0 * exp.(0.8 * x)
    y = fn(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, ExpCurveFitAlgorithm())

    @test sol.u[1] ≈ 0.8
    @test sol.u[2] ≈ log(2.0)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(val)
    end
end

@testitem "Polynomial Fit" begin
    using LinearSolve

    x = range(1, stop = 10, length = 10)

    fn(x) = 1.0 + 2.0 * x + 3.0 * x^2 + 0.5 * x^3
    y = fn.(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, PolynomialFitAlgorithm(degree = 4))

    @test sol.u[1] ≈ 1.0
    @test sol.u[2] ≈ 2.0
    @test sol.u[3] ≈ 3.0
    @test sol.u[4] ≈ 0.5
    @test sol.u[5] ≈ 0.0 atol = 1.0e-8

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(val)
    end

    @testset "ill-conditioned" begin
        true_coeffs = [80.0, -5.0e-18, -7.0e-20, -1.0e-36]
        x1 = 1.0e10 .* (0:0.1:5)
        y1 = evalpoly.(x1, (true_coeffs,))

        prob = CurveFitProblem(x1, y1)
        sol = solve(prob, PolynomialFitAlgorithm(3, QRFactorization()))

        @test sol.u[1] ≈ true_coeffs[1] rtol = 1.0e-5
        @test sol.u[2] ≈ true_coeffs[2] rtol = 1.0e-5
        @test sol.u[3] ≈ true_coeffs[3] rtol = 1.0e-5
        @test sol.u[4] ≈ true_coeffs[4] rtol = 1.0e-5

        @testset for val in (0.0, 1.5, 4.5, 10.0)
            @test sol(val) ≈ evalpoly(val, true_coeffs)
        end

        sol = solve(
            prob,
            PolynomialFitAlgorithm(3);
            assumptions = OperatorAssumptions(
                false; condition = OperatorCondition.VeryIllConditioned
            )
        )

        @test sol.u[1] ≈ true_coeffs[1] rtol = 1.0e-5
        @test sol.u[2] ≈ true_coeffs[2] rtol = 1.0e-5
        @test sol.u[3] ≈ true_coeffs[3] rtol = 1.0e-5
        @test sol.u[4] ≈ true_coeffs[4] rtol = 1.0e-5

        @testset for val in (0.0, 1.5, 4.5, 10.0)
            @test sol(val) ≈ evalpoly(val, true_coeffs)
        end
    end
end
