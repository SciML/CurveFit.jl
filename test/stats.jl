@testitem "StatsAPI Integration" begin
    using StatsAPI
    using NonlinearSolve
    using LinearAlgebra

    @testset "Linear Fit" begin
        # y = 2x + 1
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = 2.0 .* x .+ 1.0
        
        prob = CurveFitProblem(x, y)
        alg = LinearCurveFitAlgorithm()
        sol = solve(prob, alg)
        
        # Check coefficients (slope, intercept)
        slope, intercept = coef(sol)
        @test slope ≈ 2.0 atol=1e-5
        @test intercept ≈ 1.0 atol=1e-5
        
        @test sum(abs2, residuals(sol)) < 1e-10
        
        @test predict(sol) ≈ y
        @test fitted(sol) ≈ y
        
        @test nobs(sol) == 5
        @test dof(sol) == 2
        @test dof_residual(sol) == 3
        
        @test rss(sol) < 1e-10
        @test mse(sol) < 1e-10
        
        # Check predict with array input (broadcast support)
        x_new = [1.0, 2.0]
        @test predict(sol, x_new) ≈ 2.0 .* x_new .+ 1.0
        
        # Exact fit -> zero variance
        @test all(isapprox.(vcov(sol), 0; atol=1e-8))
        @test all(isapprox.(stderror(sol), 0; atol=1e-8))
    end

    @testset "Nonlinear Fit" begin
        x = range(0, 2, length=10)
        a_true = 1.0
        b_true = 2.0
        c_true = 0.5
        y = @. a_true + b_true * exp(c_true * x)
        
        @. model(u, x) = u[1] + u[2] * exp(u[3] * x)
        
        u0 = [0.5, 1.0, 0.1]
        prob = NonlinearCurveFitProblem(model, u0, x, y)
        sol = solve(prob)
        
        u = coef(sol)
        @test u[1] ≈ a_true atol=1e-2
        @test u[2] ≈ b_true atol=1e-2
        @test u[3] ≈ c_true atol=1e-2
        
        @test nobs(sol) == 10
        @test dof(sol) == 3
        @test dof_residual(sol) == 7
        
        @test rss(sol) < 1e-4
        
        # Perfect fit -> near zero errors
        @test all(stderror(sol) .< 1e-2)
        
        @test size(vcov(sol)) == (3, 3)
    end
    
    @testset "Noisy Fit Statistics" begin
        # Linear fit with noise
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = 2.0 .* x .+ 1.0 .+ [0.1, -0.1, 0.2, -0.2, 0.0]
        
        prob = CurveFitProblem(x, y)
        alg = LinearCurveFitAlgorithm()
        sol = solve(prob, alg)
        
        @test mse(sol) > 0
        @test all(stderror(sol) .> 0)
        @test size(vcov(sol)) == (2, 2)
        @test isposdef(vcov(sol))
        
        # Test confidence intervals
        cis = confint(sol)
        @test length(cis) == 2
        # True params are roughly 2.0 and 1.0. CI should cover them or be close.
        # Just checking structure and non-error
        @test cis[1][1] < cis[1][2]
        
        # Test isconverged
        @test isconverged(sol)
    end
    
    @testset "Explcit API Coverage (ExpSum)" begin
        # Test ExpSumFitAlgorithm explicitly to ensure jacobian works
        # y = k + p*exp(lam*x)
        # Truth: k=1, p=2, lam=-0.5
        x_data = collect(range(0, 5, length=20))
        y_data = @. 1.0 + 2.0 * exp(-0.5 * x_data) + 0.01 * randn()
        
        prob = CurveFitProblem(x_data, y_data)
        sol = solve(prob, ExpSumFitAlgorithm(; n=1, withconst=true))
        
        @test size(vcov(sol)) == (3, 3) # k, p, lam
        @test all(stderror(sol) .> 0)
    end

    
    @testset "Polynomial Fit Array Support" begin
        x = [1.0, 2.0, 3.0]
        y = x.^2
        # y = 0 + 0x + 1x^2
        prob = CurveFitProblem(x, y)
        sol = solve(prob, PolynomialFitAlgorithm(2))
        
        @test sol.u[1] ≈ 0.0 atol=1e-5
        @test sol.u[2] ≈ 0.0 atol=1e-5
        @test sol.u[3] ≈ 1.0 atol=1e-5
        
        @test predict(sol, [2.0, 3.0]) ≈ [4.0, 9.0]
    end
end
