@testitem "Nonlinear Least Squares: Linear Problem 1" begin
    using SciMLBase

    x = 1.0:10.0
    a0 = [3.0, 2.0, 1.0]

    fun7(a, x) = @. a[1] + a[2] * x + a[3] * x^2
    y = fun7(a0, x)

    prob = NonlinearCurveFitProblem(fun7, [0.5, 0.5, 0.5], x, y)
    sol = solve(prob)

    @test sol.coeffs ≈ a0
    @test SciMLBase.successful_retcode(sol.retcode)

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ fun7(a0, val)
    end
end

# @testitem "Nonlinear Least Squares: Nonlinear Problem 1" begin
#     x = 1.0:10.0
#     fun8(x, a) = @. a[1] + a[2] * x^a[3]
#     a0 = [3.0, 2.0, 0.7]
#     y = fun8(x, a0)

#     a = nonlinear_fit(fun8, x, [0.5, 0.5, 0.5]; target = y).u
#     @test a[1]≈a0[1] rtol=1e-7
#     @test a[2]≈a0[2] rtol=1e-7
#     @test a[3]≈a0[3] rtol=1e-7
# end

# @testitem "Nonlinear Least Squares: Linear Problem 2" begin
#     x = 1.0:10.0
#     a0 = [3.0, 2.0, 1.0]
#     fun9(x, a) = @. a[1] + a[2] * x[:, 1] + a[3] * x[:, 1]^2 - x[:, 2]
#     P = length(x)
#     X = zeros(P, 2)
#     for i in 1:P
#         X[i, 1] = x[i]
#         X[i, 2] = a0[1] + a0[2] * x[i] + a0[3] * x[i]^2
#     end

#     a = nonlinear_fit(fun9, X, [0.5, 0.5, 0.5]).u
#     @test a[1]≈a0[1] rtol=1e-7
#     @test a[2]≈a0[2] rtol=1e-7
#     @test a[3]≈a0[3] rtol=1e-7
# end

# @testitem "Nonlinear Least Squares: Nonlinear Problem 2" begin
#     x = 1.0:10.0
#     funA(x, a) = @. a[1] + a[2] * x[:, 1]^a[3] - x[:, 2]
#     a0 = [3.0, 2.0, 0.7]
#     P = length(x)
#     X = zeros(P, 2)
#     for i in 1:P
#         X[i, 1] = x[i]
#         X[i, 2] = a0[1] + a0[2] * x[i]^a0[3]
#     end

#     a = nonlinear_fit(funA, X, [0.5, 0.5, 0.5]).u
#     @test a[1]≈a0[1] rtol=1e-7
#     @test a[2]≈a0[2] rtol=1e-7
#     @test a[3]≈a0[3] rtol=1e-7
# end

# @testitem "Gauss-Newton curve fitting: Linear problem" begin
#     U = 0.5:0.5:10
#     a0 = [2.0, 1.0, 0.35]
#     E = @. sqrt(a0[1] + a0[2] * U^a0[3])

#     X = hcat(E, U)
#     funB(x, a) = @. a[1] + a[2] * x[:, 2]^a[3] - x[:, 1]^2

#     a = nonlinear_fit(funB, X, [0.5, 0.5, 0.5]).u

#     @test a[1]≈a0[1] rtol=1e-7
#     @test a[2]≈a0[2] rtol=1e-7
#     @test a[3]≈a0[3] rtol=1e-7
# end

# @testitem "Secant method NLS curve fitting: Linear problem" begin
#     x = 1.0:10.0
#     a0 = [3.0, 2.0, 1.0]
#     funC(x, a) = a[1] + a[2] * x + a[3] * x^2
#     y = funC.(x, Ref(a0))

#     a = CurveFit.secant_nls_fit(x, y, funC, [0.5, 0.5, 0.5], 1e-8, 30)

#     @test a[1]≈a0[1] rtol=1e-7
#     @test a[2]≈a0[2] rtol=1e-7
#     @test a[3]≈a0[3] rtol=1e-7
# end

# @testitem "Secant method NLS curve fitting: Nonlinear problem" begin
#     x = 1.0:10.0
#     funD(x, a) = a[1] + a[2] * x^a[3]
#     a0 = [3.0, 2.0, 0.7]
#     y = funD.(x, Ref(a0))
#     a = CurveFit.secant_nls_fit(x, y, funD, [0.5, 0.5, 0.5], 1e-8, 30)

#     @test a[1]≈a0[1] rtol=1e-7
#     @test a[2]≈a0[2] rtol=1e-7
#     @test a[3]≈a0[3] rtol=1e-7
# end
