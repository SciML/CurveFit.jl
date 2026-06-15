using CurveFit
using Test
using SciMLBase
using NonlinearSolveBase: NonlinearSolveBase

@testset "Nonlinear Least Squares: reinit!()" begin
    # Create an initial problem with a cache
    x = 1.0:10.0
    a0 = [3.0, 2.0, 0.7]

    fn(a, x) = @. a[1] + a[2] * x^a[3]
    y = fn(a0, x)

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], x, y)
    cache = CurveFit.init(prob)
    @test solve!(cache).u ≈ a0 atol = 1.0e-7

    # reinit!() the cache with different parameters and recheck the solve
    a0 = [4.0, 5.0, 0.2]
    x = 11.0:20.0
    y = fn(a0, x)

    CurveFit.reinit!(cache; u0 = [1.0, 1.0, 1.0], x, y)
    @test solve!(cache).u ≈ a0 atol = 1.0e-7

    # Repeat with an in-place model: NonlinearSolve wraps in-place Float64
    # problems in `AutoSpecializeCallable`, which `reinit!` must unwrap.
    fn(resid, a, x) = @. resid = a[1] + a[2] * x^a[3]
    cache = CurveFit.init(NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], x, y))
    @test cache.cache.prob.f.f isa NonlinearSolveBase.AutoSpecializeCallable
    CurveFit.reinit!(cache; u0 = [1.0, 1.0, 1.0], x, y)
    @test solve!(cache).u ≈ a0 atol = 1.0e-7
end
