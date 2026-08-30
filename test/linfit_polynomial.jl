using CurveFit
using CommonSolve: init, solve, solve!
using Test
using LinearSolve

@testset "Polynomial Fit" begin
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

    @test sol.resid ≈ y .- sol.(x)

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

        @testset "default algorithm" begin
            sol = solve(prob, PolynomialFitAlgorithm(3))

            @test sol.u[1] ≈ true_coeffs[1] rtol = 1.0e-5
            @test sol.u[2] ≈ true_coeffs[2] rtol = 1.0e-5
            @test sol.u[3] ≈ true_coeffs[3] rtol = 1.0e-5
            @test sol.u[4] ≈ true_coeffs[4] rtol = 1.0e-5
        end

        @testset "degree 4, x ~ 1e6" begin
            true_coeffs4 = [3.0, -2.0e-6, 1.0e-12, -5.0e-19, 2.0e-25]
            x2 = 1.0e6 .* (0:0.1:5)
            prob2 = CurveFitProblem(x2, evalpoly.(x2, (true_coeffs4,)))
            sol = solve(prob2, PolynomialFitAlgorithm(4))

            @testset for i in eachindex(true_coeffs4)
                @test sol.u[i] ≈ true_coeffs4[i] rtol = 1.0e-5
            end
        end
    end

    @testset "integer x" begin
        xi = 1:10
        yi = evalpoly.(xi, ([1.0, 2.0, 3.0],))
        sol = solve(CurveFitProblem(xi, yi), PolynomialFitAlgorithm(2))

        @test sol.u[1] ≈ 1.0
        @test sol.u[2] ≈ 2.0
        @test sol.u[3] ≈ 3.0
    end
end
