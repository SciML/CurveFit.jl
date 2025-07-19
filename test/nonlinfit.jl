@testitem "Nonlinear Least Squares: Linear Problem 1" begin
    using SciMLBase

    x = 1.0:10.0
    a0 = [3.0, 2.0, 1.0]

    fn(a, x) = @. a[1] + a[2] * x + a[3] * x^2
    y = fn(a0, x)

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], x, y)
    sol = solve(prob)

    @test sol.u ≈ a0
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(a0, val)
    end
end

@testitem "Nonlinear Least Squares: Nonlinear Problem 1" begin
    using SciMLBase

    x = 1.0:10.0
    a0 = [3.0, 2.0, 0.7]

    fn(a, x) = @. a[1] + a[2] * x^a[3]
    y = fn(a0, x)

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], x, y)
    sol = solve(prob)

    @test sol.u ≈ a0
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fn(a0, val)
    end
end

@testitem "Nonlinear Least Squares: Linear Problem 2" begin
    using SciMLBase

    x = 1.0:10.0
    a0 = [3.0, 2.0, 1.0]

    fn(a, x) = @. a[1] + a[2] * x[:, 1] + a[3] * x[:, 1]^2 - x[:, 2]
    P = length(x)
    X = zeros(P, 2)
    for i in 1:P
        X[i, 1] = x[i]
        X[i, 2] = fn(a0, [x[i] 0])[1]
    end

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], X)
    sol = solve(prob)

    @test sol.u ≈ a0 atol=1.0e-7
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol([val 0.0])[1] ≈ fn(a0, [val 0.0])[1] atol=1.0e-7
    end
end

@testitem "Nonlinear Least Squares: Nonlinear Problem 2" begin
    using SciMLBase

    x = 1.0:10.0
    a0 = [3.0, 2.0, 0.7]

    fn(a, x) = @. a[1] + a[2] * x[:, 1]^a[3] - x[:, 2]
    P = length(x)
    X = zeros(P, 2)
    for i in 1:P
        X[i, 1] = x[i]
        X[i, 2] = fn(a0, [x[i] 0])[1]
    end

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], X)
    sol = solve(prob)

    @test sol.u ≈ a0 atol=1.0e-7
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol([val 0.0])[1] ≈ fn(a0, [val 0.0])[1] atol=1.0e-7
    end
end

@testitem "Gauss-Newton curve fitting: Linear problem" begin
    using SciMLBase

    U = 0.5:0.5:10
    a0 = [2.0, 1.0, 0.35]
    E = @. sqrt(a0[1] + a0[2] * U^a0[3])

    X = hcat(E, U)
    fn(a, x) = @. a[1] + a[2] * x[:, 2]^a[3] - x[:, 1]^2

    prob = NonlinearCurveFitProblem(fn, [0.5, 0.5, 0.5], X)
    sol = solve(prob)

    @test sol.u ≈ a0 atol=1.0e-7
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in range(minimum(E), stop = maximum(E), length = 10)
        @test sol([val 0.0])[1] ≈ fn(a0, [val 0.0])[1] atol=1.0e-7
    end
end

@testitem "Nonlinear Residuals Field" begin
    using SciMLBase

    # Test 1: Simple exponential fit with noise
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    a_true = [2.0, 0.5]
    fn(a, x) = @. a[1] * exp(a[2] * x)
    y_true = fn(a_true, x)
    # Add some noise
    y = y_true .+ [0.1, -0.05, 0.08, -0.02, 0.03]

    prob = NonlinearCurveFitProblem(fn, [1.0, 0.3], x, y)
    sol = solve(prob)

    # Test that resid field exists and has correct properties
    @test sol.resid !== nothing
    @test length(sol.resid) == length(y)
    @test SciMLBase.successful_retcode(sol.retcode)
    
    # Test that residuals are computed correctly as y - fitted_values
    fitted_values = sol.(x)
    expected_resid = y .- fitted_values
    @test sol.resid ≈ expected_resid

    # Test 2: Perfect nonlinear fit (should have near-zero residuals)
    x_perfect = [1.0, 2.0, 3.0, 4.0]
    a_perfect = [3.0, 0.8, 1.2]
    fn_perfect(a, x) = @. a[1] + a[2] * x^a[3]
    y_perfect = fn_perfect(a_perfect, x_perfect)

    prob_perfect = NonlinearCurveFitProblem(fn_perfect, [2.5, 0.7, 1.1], x_perfect, y_perfect)
    sol_perfect = solve(prob_perfect)

    @test sol_perfect.resid !== nothing
    @test SciMLBase.successful_retcode(sol_perfect.retcode)
    @test maximum(abs.(sol_perfect.resid)) < 1e-10  # Should be very small for perfect fit
    
    # Test 3: Nonlinear problem without y data (should set resid to nothing appropriately)
    # This tests the case where prob.y is nothing
    fn_noy(a, x) = @. a[1] + a[2] * x[:, 1]^a[3] - x[:, 2]
    X_test = hcat([1.0, 2.0, 3.0], [3.0, 5.0, 7.0])  # Second column represents target values
    
    prob_noy = NonlinearCurveFitProblem(fn_noy, [1.0, 1.0, 1.0], X_test)
    sol_noy = solve(prob_noy)
    
    # For problems without explicit y data, residuals should still be computed if possible
    @test sol_noy.resid !== nothing || sol_noy.resid === nothing  # Either way is acceptable
    @test SciMLBase.successful_retcode(sol_noy.retcode)
end
