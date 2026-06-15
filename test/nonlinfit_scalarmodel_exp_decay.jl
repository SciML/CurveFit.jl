using CurveFit
using Test
using SciMLBase

@testset "ScalarModel: Exponential decay (LsqFit-style migration)" begin
    # This test demonstrates migration from LsqFit.jl style
    # LsqFit: model(x, p) = p[1] * exp(-x * p[2])
    # CurveFit: model(p, x) = p[1] * exp(-x * p[2])

    x = collect(range(0, stop = 10, length = 20))
    true_params = [2.5, 0.3]

    # Scalar model function (parameter order: params first, then x)
    model(p, x) = p[1] * exp(-x * p[2])

    # Generate y data
    y = model.(Ref(true_params), x)

    # Use ScalarModel wrapper
    prob = NonlinearCurveFitProblem(ScalarModel(model), [1.0, 0.1], x, y)
    sol = solve(prob)

    @test sol.u ≈ true_params atol = 1.0e-6
    @test SciMLBase.successful_retcode(sol.retcode)

    # Verify predictions
    @test sol(5.0) ≈ model(true_params, 5.0)
end
