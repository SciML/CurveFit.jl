@testitem "StatsAPI Integration" begin
    using StatsAPI
    using NonlinearSolve
    using LinearAlgebra
    using LinearSolve

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
    
    @testset "Comprehensive API Coverage (All Types)" begin
        x_data = collect(range(0.5, 5.0, length=20)) # Avoid 0 for log/power models
        
        # 1. Log Fit: y = a + b*ln(x)
        # Truth: a=1, b=2
        y_log = @. 1.0 + 2.0 * log(x_data) + 0.01 * randn()
        prob_log = CurveFitProblem(x_data, y_log)
        sol_log = solve(prob_log, LogCurveFitAlgorithm())
        @test size(vcov(sol_log)) == (2, 2)
        @test all(stderror(sol_log) .> 0)

        # 2. Power Fit: y = b * x^a  -> ln(y) = ln(b) + a*ln(x)
        # Truth: b=2, a=0.5
        y_pow = @. 2.0 * x_data^0.5 + 0.01 * randn()
        # Power fit linearizes, so y must be positive
        prob_pow = CurveFitProblem(x_data, abs.(y_pow)) 
        sol_pow = solve(prob_pow, PowerCurveFitAlgorithm())
        @test size(vcov(sol_pow)) == (2, 2)
        @test all(stderror(sol_pow) .> 0)
        
        # 3. Exp Fit: y = b * exp(a*x) -> ln(y) = ln(b) + a*x
        # Truth: b=2, a=0.5
        y_exp = @. 2.0 * exp(0.5 * x_data) + 0.01 * randn()
        prob_exp = CurveFitProblem(x_data, abs.(y_exp))
        sol_exp = solve(prob_exp, ExpCurveFitAlgorithm())
        @test size(vcov(sol_exp)) == (2, 2)
        @test all(stderror(sol_exp) .> 0)
        
        # 4. Rational Fit (Linear): y = (p0 + p1*x)/(1 + q1*x) 
        # Truth: p0=1, p1=2, q1=0.5
        # y * (1 + 0.5x) = 1 + 2x
        sol_rat = solve(prob_log, RationalPolynomialFitAlgorithm(num_degree=1, den_degree=1))
        # Params: p0, p1, q1 (3 params)
        @test size(vcov(sol_rat)) == (3, 3)
        @test all(stderror(sol_rat) .> 0)
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

        # Modified King Fit (E^2 = A + B * U^n)
        # x corresponds to E (Voltage) in Jacobian logic, but input data order is (U, E^2)?
        # User creates CurveFitProblem(x, y). 
        # King assumes (Voltage, Velocity).
        # My implementation: x=E, y=U in Jacobian.
        # But in verify script I passed (U, E^2).
        # So I will replicate verify script logic here.
        A_true = 1.0; B_true = 2.0; n_true = 0.5
        x_k = collect(1.0:0.5:5.0)
        y_k = @. A_true + B_true * x_k^n_true
        prob_mod_king = CurveFitProblem(x_k, y_k)
        sol_mod_king = solve(prob_mod_king, ModifiedKingCurveFitAlgorithm())
        
        @test size(vcov(sol_mod_king)) == (3, 3)
        @test all(stderror(sol_mod_king) .> 0)
    end

    
    @testset "Polynomial Fit Array Support" begin
        # Use overdetermined system (5 points, 3 params) to avoid division by zero in MSE
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = x.^2
        # y = 0 + 0x + 1x^2
        prob = CurveFitProblem(x, y)
        sol = solve(prob, PolynomialFitAlgorithm(2))
        
        @test sol.u[1] ≈ 0.0 atol=1e-5
        @test sol.u[2] ≈ 0.0 atol=1e-5
        @test sol.u[3] ≈ 1.0 atol=1e-5
        
        @test predict(sol, [2.0, 3.0]) ≈ [4.0, 9.0]
        
        @test size(vcov(sol)) == (3, 3) 
        @test all(stderror(sol) .>= 0)
    end
    
    @testset "Interface & Solver Options" begin
        # 1. In-place evaluation (Issue #47)
        x_ip = [1.0, 2.0, 3.0]
        y_ip = [2.0, 4.0, 6.0]
        
        function model_inplace!(resid, u, x)
             # y = u[1]*x
             # The solver expects predictions in `resid`
             @. resid = u[1] * x
        end
        
        prob_ip = NonlinearCurveFitProblem(
            NonlinearFunction{true}(model_inplace!; resid_prototype = similar(y_ip)),
            [1.0], # u0
            x_ip,
            y_ip
        )
        sol_ip = solve(prob_ip)
        @test sol_ip.u[1] ≈ 2.0
        @test isconverged(sol_ip)
        
        # 2. LinearSolve algorithm choices (Issue #49)
        x_ls = collect(1.0:5.0) # Use Float64 to avoid InexactError in QR
        y_ls = 2.0 .* x_ls .+ 1.0
        
        # Pass LinearSolve.QRFactorization explicitly
        prob_ls = CurveFitProblem(x_ls, y_ls)
        algo_ls = PolynomialFitAlgorithm(; degree=1, linsolve_algorithm=LinearSolve.QRFactorization())
        sol_ls = solve(prob_ls, algo_ls)
        
    
        @test sol_ls.u[1] ≈ 1.0 # intercept
        @test sol_ls.u[2] ≈ 2.0 # slope
        @test size(vcov(sol_ls)) == (2, 2)

        # 3. LevenbergMarquardt with linsolve choices (Issue #49)
        # Verify passing explicit linear solver to nonlinear algorithm
        x_nl_ls = collect(1.0:10.0)
        y_nl_ls = @. 1.0 + 2.0 * exp(-0.5 * x_nl_ls)
        @. model_nl_ls(u, x) = u[1] + u[2] * exp(u[3] * x)
        prob_nl_ls = NonlinearCurveFitProblem(model_nl_ls, [0.5, 0.5, -0.1], x_nl_ls, y_nl_ls)
        
        # Test Cholesky (fast but requires PD, usually OK for LM)
        sol_lm_chol = solve(prob_nl_ls, LevenbergMarquardt(linsolve=LinearSolve.CholeskyFactorization()))
        @test sol_lm_chol.u[1] ≈ 1.0 atol=1e-3
        @test size(vcov(sol_lm_chol)) == (3, 3)
    end
    
    @testset "Extended Algorithm Test Coverage" begin
        x_data = collect(1.0:0.5:5.0)
        y_data = 2.0 .* x_data .^ 0.5 .+ 0.1
        
        # 1. King Fit (Linear alias)
        # y = A + B*sqrt(x) -> Linear fit with yfun=sqrt? No, King is xfun=abs2, yfun=sqrt? 
        # Wait, King is LinearCurveFitAlgorithm(; xfun = abs2, yfun = sqrt)
        # Check source: KingCurveFitAlgorithm() = LinearCurveFitAlgorithm(; xfun = abs2, yfun = sqrt)
        # This implies sqrt(y) = A + B*x^2. 
        # This is strictly linear in transformed space, so vcov should work.
        prob_king = CurveFitProblem(x_data, y_data)
        sol_king = solve(prob_king, KingCurveFitAlgorithm())
        @test size(vcov(sol_king)) == (2, 2)
        
        # 2. Rational Fit (Nonlinear)
        # Default defaults linear?
        # Providing u0 forces nonlinear path if alg is generic?
        # The struct default alg is QR/Linear?
        # RationalPolynomialFitAlgorithm(; alg=AbstractLinearAlgorithm, ...)
        
        # Let's force nonlinear by providing u0 or using a nonlinear solver?
        # init(...; u0) is not supported for linear path?
        # src/rationalfit.jl: if alg.alg isa AbstractLinearAlgorithm ... @assert prob.u0 === nothing
        
        # So to test nonlinear vcov for rational, we must use a nonlinear solver for `alg`.
        # RationalPolynomialFitAlgorithm(alg=LevenbergMarquardt(), ...)
        # Assuming LevenbergMarquardt is available from NonlinearSolve?
        # Or generic NewtonRaphson?
        
        prob_rat_nl = CurveFitProblem(x_data, y_data) # u0 is optional in CurveFitProblem but required for Nonlinear
        # Wait, RationalPolynomialFitAlgorithm automatically generates u0 if not provided in nonlinear path?
        # src/rationalfit.jl:111: u0 = if cache.initial_guess_cache ... else cache.prob.u0
        # And it uses a linear fit for initial guess!
        
        # So if we pass a nonlinear algorithm to RationalPolynomialFitAlgorithm, it should stay in nonlinear mode.
        sol_rat_nl = solve(prob_rat_nl, RationalPolynomialFitAlgorithm(; alg=NewtonRaphson(), num_degree=1, den_degree=1))
        
        @test size(vcov(sol_rat_nl)) == (3, 3)
        @test all(stderror(sol_rat_nl) .> 0)
    end
end
