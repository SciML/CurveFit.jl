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

# Regression test for https://github.com/SciML/CurveFit.jl/issues/69
@testitem "Issue #69: y ~ a/x with noisy data" begin
    using SciMLBase

    # Original issue: fitting y ~ a/x failed with NaN when data had noise
    x = [1.0, 2.0, 3.0, 4.0, 5.0]

    # Test case 1: Exact data (should work perfectly)
    fn(a, x) = @. a[1] / x
    y_exact = [1.0, 0.5, 1 / 3, 0.25, 0.2]
    prob1 = NonlinearCurveFitProblem(fn, [0.1], x, y_exact)
    sol1 = solve(prob1)

    @test sol1.u[1] ≈ 1.0 atol = 1e-6
    @test SciMLBase.successful_retcode(sol1.retcode)
    @test !isnan(sol1.u[1])

    # Test case 2: Data with tiny noise (previously failed with NaN)
    y_noisy = [1.0, 0.5, 1 / 3 + 0.00000001, 0.25, 0.2]
    prob2 = NonlinearCurveFitProblem(fn, [0.1], x, y_noisy)
    sol2 = solve(prob2)

    @test sol2.u[1] ≈ 1.0 atol = 1e-5
    @test SciMLBase.successful_retcode(sol2.retcode)
    @test !isnan(sol2.u[1])

    # Test case 3: y ~ ax exact data
    fn2(a, x) = @. a[1] * x
    y2_exact = [1.0, 2.0, 3.0, 4.0, 5.0]
    prob3 = NonlinearCurveFitProblem(fn2, [0.1], x, y2_exact)
    sol3 = solve(prob3)

    @test sol3.u[1] ≈ 1.0 atol = 1e-6
    @test SciMLBase.successful_retcode(sol3.retcode)
    @test !isnan(sol3.u[1])

    # Test case 4: y ~ ax with tiny noise (previously failed with NaN)
    y2_noisy = [1.0, 2.0, 3.000001, 4.0, 5.0]
    prob4 = NonlinearCurveFitProblem(fn2, [0.1], x, y2_noisy)
    sol4 = solve(prob4)

    @test sol4.u[1] ≈ 1.0 atol = 1e-5
    @test SciMLBase.successful_retcode(sol4.retcode)
    @test !isnan(sol4.u[1])
end

@testitem "Issue #69: Robustness with larger noise" begin
    using SciMLBase

    x = [1.0, 2.0, 3.0, 4.0, 5.0]

    # y ~ a/x with larger noise
    fn(a, x) = @. a[1] / x
    y_noisy = [1.0 + 0.01, 0.5 - 0.01, 1 / 3 + 0.01, 0.25, 0.2 - 0.005]
    prob = NonlinearCurveFitProblem(fn, [0.1], x, y_noisy)
    sol = solve(prob)

    @test sol.u[1] ≈ 1.0 atol = 0.1  # Looser tolerance for noisy data
    @test SciMLBase.successful_retcode(sol.retcode)
    @test !isnan(sol.u[1])
end
