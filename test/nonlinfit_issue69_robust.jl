using CurveFit
using Test
using SciMLBase

@testset "Issue #69: Robustness with larger noise" begin
    x = [1.0, 2.0, 3.0, 4.0, 5.0]

    # y ~ a/x with larger noise
    fn(a, x) = @. a[1] / x
    y_noisy = [1.0 + 0.01, 0.5 - 0.01, 1 / 3 + 0.01, 0.25, 0.2 - 0.005]
    prob = NonlinearCurveFitProblem(fn, [0.1], x, y_noisy)
    sol = solve(prob)

    @test sol.u[1] ≈ 1.0 atol = 0.1  # Looser tolerance for noisy data
    @test SciMLBase.successful_retcode(sol.retcode)
    @test !isnan(sol.u[1])
end
