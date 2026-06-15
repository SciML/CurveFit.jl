using CurveFit
using Test
using SciMLBase
using NonlinearSolveFirstOrder: LevenbergMarquardt, GaussNewton, TrustRegion

@testset "GenericNonlinearCurveFitCache show" begin
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
