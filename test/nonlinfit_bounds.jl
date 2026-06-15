using CurveFit
using Test
using SciMLBase

@testset "Nonlinear Least Squares: Bounds constrain solution" begin
    # True params are [3.0, 1.0] but we constrain p[1] within [0.0, 2.0]
    fn(a, x) = @. a[1] * exp(a[2] * x)
    x = collect(range(0, 2, length = 20))
    y = fn.(Ref([3, 1]), x)

    lb = [0.0, -Inf]
    ub = [2.0, Inf]

    prob = NonlinearCurveFitProblem(fn, [1.0, 0.5], x, y; lb, ub)

    # Verify bounds are stored on the CurveFitProblem
    @test prob.lb == lb
    @test prob.ub == ub

    # Test solve() path
    sol = solve(prob)
    @test sol.u[1] <= 2.0

    # Test init+solve!() path
    cache = CurveFit.init(prob)
    sol2 = solve!(cache)
    @test sol2.u[1] <= 2.0
end
