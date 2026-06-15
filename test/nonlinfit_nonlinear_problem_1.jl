using CurveFit
using Test
using SciMLBase

@testset "Nonlinear Least Squares: Nonlinear Problem 1" begin
    x = 1.0:10.0
    a0 = [3.0, 2.0, 0.7]

    fn(a, x) = @. a[1] + a[2] * x^a[3]
    y = fn(a0, x)

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], x, y)
    sol = solve(prob)

    @test sol.u ≈ a0
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(a0, val)
    end
end
