using CurveFit
using Test
using SciMLBase

@testset "CurveFitSolution" begin
    x = 1:10
    fn(a, x) = @. 1.0 + 2.0 * x + a[1]
    y = fn([1.0], x)
    linear_sol = solve(CurveFitProblem(x, y), LinearCurveFitAlgorithm())
    nonlinear_sol = solve(NonlinearCurveFitProblem(fn, x, y, [1.0]))

    @test SciMLBase.successful_retcode(linear_sol)

    # Smoke test
    @test contains(repr(MIME"text/plain"(), linear_sol), "residuals mean:")
    @test contains(repr(MIME"text/plain"(), nonlinear_sol), "residuals mean:")
end
