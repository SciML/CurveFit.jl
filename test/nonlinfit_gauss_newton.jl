using CurveFit
using Test
using SciMLBase

@testset "Gauss-Newton curve fitting: Linear problem" begin
    U = 0.5:0.5:10
    a0 = [2.0, 1.0, 0.35]
    E = @. sqrt(a0[1] + a0[2] * U^a0[3])

    X = hcat(E, U)
    fn(a, x) = @. a[1] + a[2] * x[:, 2]^a[3] - x[:, 1]^2

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], X)
    sol = solve(prob)

    @test sol.u ≈ a0 atol = 1.0e-7
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in range(minimum(E), stop = maximum(E), length = 10)
        @test sol([val 0.0])[1] ≈ fn(a0, [val 0.0])[1] atol = 1.0e-7
    end
end
