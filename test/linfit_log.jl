using CurveFit
using Test

@testset "Log Fit" begin
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
