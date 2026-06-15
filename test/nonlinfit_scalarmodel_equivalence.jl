using CurveFit
using Test
using SciMLBase

@testset "ScalarModel: Equivalence with vectorized @. form" begin
    # Verify that ScalarModel produces the same results as @. form

    x = 1.0:10.0
    a0 = [3.0, 2.0, 0.7]
    u0 = [0.5, 0.5, 0.5]

    # Vectorized form (standard CurveFit style)
    fn_vec(a, x) = @. a[1] + a[2] * x^a[3]
    y = fn_vec(a0, x)

    prob_vec = NonlinearCurveFitProblem(fn_vec, u0, x, y)
    sol_vec = solve(prob_vec)

    # Scalar form with ScalarModel
    fn_scalar(a, x) = a[1] + a[2] * x^a[3]

    prob_scalar = NonlinearCurveFitProblem(ScalarModel(fn_scalar), u0, x, y)
    sol_scalar = solve(prob_scalar)

    # Both should give the same result
    @test sol_vec.u ≈ sol_scalar.u
    @test SciMLBase.successful_retcode(sol_vec.retcode)
    @test SciMLBase.successful_retcode(sol_scalar.retcode)

    # Both should evaluate the same at any point
    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol_vec(val) ≈ sol_scalar(val)
    end
end
