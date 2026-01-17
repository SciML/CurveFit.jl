@testitem "CurveFitSolution" begin
    using SciMLBase

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

@testitem "GenericNonlinearCurveFitCache show" begin
    using CurveFit
    using SciMLBase
    using NonlinearSolve: LevenbergMarquardt, GaussNewton, TrustRegion

    x = collect(1.0:10.0)
    fn(a, x) = @. a[1] + a[2] * x
    y = fn([1.0, 2.0], x)
    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5], x, y)

    # Smoke tests with various algorithms
    cache = init(prob)
    @test_nowarn repr(MIME"text/plain"(), cache)

    cache = init(prob, LevenbergMarquardt())
    @test_nowarn repr(MIME"text/plain"(), cache)

    # Smoke test with a solved problem
    cache = init(prob)
    sol = solve!(cache)
    @test SciMLBase.successful_retcode(sol)
    @test_nowarn repr(MIME"text/plain"(), cache)
end
