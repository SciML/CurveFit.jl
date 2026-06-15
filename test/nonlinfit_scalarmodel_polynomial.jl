using CurveFit
using Test
using SciMLBase

# Tests for ScalarModel - Issue #46
@testset "ScalarModel: Basic polynomial fitting" begin
    x = 1.0:10.0
    a0 = [3.0, 2.0, 1.0]

    # Scalar function (no @.)
    fn_scalar(a, x) = a[1] + a[2] * x + a[3] * x^2

    # Generate y data by broadcasting the scalar function
    y = fn_scalar.(Ref(a0), x)

    # Use ScalarModel wrapper
    prob = NonlinearCurveFitProblem(ScalarModel(fn_scalar), [0.5, 0.5, 0.5], x, y)
    @test SciMLBase.isinplace(prob.nlfunc)
    sol = solve(prob)
    @test sol(x) isa Vector

    @test sol.u ≈ a0
    @test SciMLBase.successful_retcode(sol.retcode)

    # Test single-point evaluation
    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn_scalar(a0, val)
    end
end
