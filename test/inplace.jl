@testitem "In-place evaluation" begin
    using CurveFit
    using NonlinearSolve
    
    X = collect(1.0:10.0)
    θ_ref = [3.0, 2.0, 1.0]

    function fn(Y, θ, X)
    @assert size(Y) == size(X) # Fails!
    @. Y = θ[1] + θ[2] * X + θ[3] * X^2
    end
    Y = zero(X)
    fn(Y, θ_ref, X)

    nonfn = NonlinearFunction(fn)

    prob = NonlinearCurveFitProblem(nonfn, [0.5, 0.5, 0.5], X, Y)
    sol = solve(prob, LevenbergMarquardt())

    @test sol.u ≈ θ_ref

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        vec = [val]
        y1, y2 = zero(vec), zero(vec)
        @test sol(y1, vec) ≈ fn(y2, θ_ref, vec)
    end
end