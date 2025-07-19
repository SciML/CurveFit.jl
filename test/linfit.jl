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

@testitem "Residuals Field" begin
    # Test that residuals are properly computed and stored
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = [2.1, 4.0, 5.9, 8.1, 9.9]  # y ≈ 2x with some noise

    prob = CurveFitProblem(x, y)
    sol = solve(prob, LinearCurveFitAlgorithm())

    # Test that resid field exists and is not nothing
    @test sol.resid !== nothing
    @test length(sol.resid) == length(y)
    
    # Test that residuals are computed correctly as y - fitted_values
    fitted_values = sol.(x)
    expected_resid = y .- fitted_values
    @test sol.resid ≈ expected_resid
    
    # Test with perfect fit (should have near-zero residuals)
    x_perfect = [1.0, 2.0, 3.0, 4.0]
    y_perfect = [3.0, 5.0, 7.0, 9.0]  # exactly y = 2x + 1
    
    prob_perfect = CurveFitProblem(x_perfect, y_perfect)
    sol_perfect = solve(prob_perfect, LinearCurveFitAlgorithm())
    
    @test sol_perfect.resid !== nothing
    @test maximum(abs.(sol_perfect.resid)) < 1e-14  # Should be essentially zero
    
    # Test polynomial fit residuals
    x_poly = [1.0, 2.0, 3.0, 4.0, 5.0]
    y_poly = [1.0, 4.0, 9.0, 16.0, 25.0]  # exactly x^2
    
    prob_poly = CurveFitProblem(x_poly, y_poly)
    sol_poly = solve(prob_poly, PolynomialFitAlgorithm(degree=2))
    
    @test sol_poly.resid !== nothing
    @test maximum(abs.(sol_poly.resid)) < 1e-14  # Should be essentially zero for perfect fit
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
    @test sol.u[5]≈0.0 atol=1e-8

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(val)
    end

    @testset "ill-conditioned" begin
        true_coeffs = [80.0, -5e-18, -7e-20, -1e-36]
        x1 = 1e10 .* (0:0.1:5)
        y1 = evalpoly.(x1, (true_coeffs,))

        prob = CurveFitProblem(x1, y1)
        sol = solve(prob, PolynomialFitAlgorithm(3, QRFactorization()))

        @test sol.u[1]≈true_coeffs[1] rtol=1e-5
        @test sol.u[2]≈true_coeffs[2] rtol=1e-5
        @test sol.u[3]≈true_coeffs[3] rtol=1e-5
        @test sol.u[4]≈true_coeffs[4] rtol=1e-5

        @testset for val in (0.0, 1.5, 4.5, 10.0)
            @test sol(val) ≈ evalpoly(val, true_coeffs)
        end

        sol = solve(prob,
            PolynomialFitAlgorithm(3);
            assumptions = OperatorAssumptions(
                false; condition = OperatorCondition.VeryIllConditioned
            )
        )

        @test sol.u[1]≈true_coeffs[1] rtol=1e-5
        @test sol.u[2]≈true_coeffs[2] rtol=1e-5
        @test sol.u[3]≈true_coeffs[3] rtol=1e-5
        @test sol.u[4]≈true_coeffs[4] rtol=1e-5

        @testset for val in (0.0, 1.5, 4.5, 10.0)
            @test sol(val) ≈ evalpoly(val, true_coeffs)
        end
    end
end
