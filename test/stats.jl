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

    @test coef(sol) ≈ θ_ref
    @test nobs(sol) == 10
    @test dof(sol) == 3
    @test dof_residual(sol) == 7
    @test rss(sol) <  1e-4
    @test mse(sol) < 1e-4
    @test vcov(sol) isa AbstractMatrix
    @test size(vcov(sol)) == (length(sol.u), length(sol.u))
    @test length(stderror(sol)) == length(coef(sol))
    @test all(stderror(sol) .>= 0) 
    @test length(margin_of_error(sol, 0.05)) == length(coef(sol))
    @test all(margin_of_error(sol, 0.05) .>= 0)
    @test length(confint(sol)) == length(coef(sol)) 
    ci = confint(sol)
    @test all(first.(ci) .<= coef(sol) .<= last.(ci))
end    
