using CurveFit
using Test
using SciMLBase

# Regression test for https://github.com/SciML/CurveFit.jl/issues/69
@testset "Issue #69: y ~ a/x with noisy data" begin
    # Original issue: fitting y ~ a/x failed with NaN when data had noise
    x = [1.0, 2.0, 3.0, 4.0, 5.0]

    # Test case 1: Exact data (should work perfectly)
    fn(a, x) = @. a[1] / x
    y_exact = [1.0, 0.5, 1 / 3, 0.25, 0.2]
    prob1 = NonlinearCurveFitProblem(fn, [0.1], x, y_exact)
    sol1 = solve(prob1)

    @test sol1.u[1] ≈ 1.0 atol = 1.0e-6
    @test SciMLBase.successful_retcode(sol1.retcode)
    @test !isnan(sol1.u[1])

    # Test case 2: Data with tiny noise (previously failed with NaN)
    y_noisy = [1.0, 0.5, 1 / 3 + 0.00000001, 0.25, 0.2]
    prob2 = NonlinearCurveFitProblem(fn, [0.1], x, y_noisy)
    sol2 = solve(prob2)

    @test sol2.u[1] ≈ 1.0 atol = 1.0e-5
    @test SciMLBase.successful_retcode(sol2.retcode)
    @test !isnan(sol2.u[1])

    # Test case 3: y ~ ax exact data
    fn2(a, x) = @. a[1] * x
    y2_exact = [1.0, 2.0, 3.0, 4.0, 5.0]
    prob3 = NonlinearCurveFitProblem(fn2, [0.1], x, y2_exact)
    sol3 = solve(prob3)

    @test sol3.u[1] ≈ 1.0 atol = 1.0e-6
    @test SciMLBase.successful_retcode(sol3.retcode)
    @test !isnan(sol3.u[1])

    # Test case 4: y ~ ax with tiny noise (previously failed with NaN)
    y2_noisy = [1.0, 2.0, 3.000001, 4.0, 5.0]
    prob4 = NonlinearCurveFitProblem(fn2, [0.1], x, y2_noisy)
    sol4 = solve(prob4)

    @test sol4.u[1] ≈ 1.0 atol = 1.0e-5
    @test SciMLBase.successful_retcode(sol4.retcode)
    @test !isnan(sol4.u[1])
end
