using CurveFit
using Test
using SciMLBase

@testset "Nonlinear Least Squares: Linear Problem 2" begin
    x = 1.0:10.0
    a0 = [3.0, 2.0, 1.0]

    fn(a, x) = @. a[1] + a[2] * x[:, 1] + a[3] * x[:, 1]^2 - x[:, 2]
    P = length(x)
    X = zeros(P, 2)
    for i in 1:P
        X[i, 1] = x[i]
        X[i, 2] = fn(a0, [x[i] 0])[1]
    end

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], X)
    sol = solve(prob)

    @test sol.u ≈ a0 atol = 1.0e-7
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol([val 0.0])[1] ≈ fn(a0, [val 0.0])[1] atol = 1.0e-7
    end
end
