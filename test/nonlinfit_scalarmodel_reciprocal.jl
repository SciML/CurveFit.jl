using CurveFit
using Test
using SciMLBase

@testset "ScalarModel: Reciprocal function" begin
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    a0 = [2.0]

    # Scalar function: y = a/x
    fn_scalar(a, x) = a[1] / x
    y = fn_scalar.(Ref(a0), x)

    prob = NonlinearCurveFitProblem(ScalarModel(fn_scalar), [0.5], x, y)
    sol = solve(prob)

    @test sol.u[1] ≈ a0[1] atol = 1.0e-6
    @test SciMLBase.successful_retcode(sol.retcode)
end
