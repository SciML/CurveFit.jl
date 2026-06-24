using CurveFit
using Test

@testset "Power Fit" begin
    x = range(1, stop = 10, length = 10)

    fn(x) = 2.0 * x .^ 0.8
    y = fn(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, PowerCurveFitAlgorithm())

    @test sol.u[1] ≈ 0.8
    @test sol.u[2] ≈ 2.0

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(val)
    end

    # Sigma not supported
    prob_sigma = CurveFitProblem(x, y; sigma = ones(length(y)))
    @test_throws AssertionError solve(prob_sigma, PowerCurveFitAlgorithm())
end
