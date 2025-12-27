@testitem "StatsAPI interface" begin
    using CurveFit
    using NonlinearSolve

    X = collect(1.0:10.0)
    θ_ref = [3.0, 2.0, 1.0]

    function fn(Y, θ, X)
        @. Y = θ[1] + θ[2] * exp(θ[3] * X)
    end
    Y = zero(X)
    fn(Y, θ_ref, X)

    nonfn = NonlinearFunction(fn)

    prob = NonlinearCurveFitProblem(nonfn, [0.5, 0.5, 0.5], X, Y)
    sol = solve(prob, LevenbergMarquardt())

    @test residuals(sol)
    @test coef(sol)
    @test nobs(sol)
    @test dof(sol)
    @test rss(sol)
    @test mse(sol)
    @test vcov(sol)
    @test stderror(sol)
    @test margin_of_error(sol, 0.05)
    @test confint(sol)
end    