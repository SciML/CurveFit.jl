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

# @testitem "Secant method NLS curve fitting: Linear problem" begin
#     x = 1.0:10.0
#     a0 = [3.0, 2.0, 1.0]
#     fn(x, a) = a[1] + a[2] * x + a[3] * x^2
#     y = fn.(x, Ref(a0))

#     a = CurveFit.secant_nls_fit(x, y, fn, [0.5, 0.5, 0.5], 1e-8, 30)

#     @test a[1]≈a0[1] rtol=1e-7
#     @test a[2]≈a0[2] rtol=1e-7
#     @test a[3]≈a0[3] rtol=1e-7
# end

# @testitem "Secant method NLS curve fitting: Nonlinear problem" begin
#     x = 1.0:10.0
#     fn(x, a) = a[1] + a[2] * x^a[3]
#     a0 = [3.0, 2.0, 0.7]
#     y = fn.(x, Ref(a0))
#     a = CurveFit.secant_nls_fit(x, y, fn, [0.5, 0.5, 0.5], 1e-8, 30)

#     @test a[1]≈a0[1] rtol=1e-7
#     @test a[2]≈a0[2] rtol=1e-7
#     @test a[3]≈a0[3] rtol=1e-7
# end
