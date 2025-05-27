@testitem "Linear Rational fit" begin
    using LinearSolve, LinearAlgebra

    x = range(1, stop = 10, length = 10)
    r = CurveFit.RationalPolynomial([1.0, 0.0, -2.0], [1.0, 2.0, 3.0])
    y = r.(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, RationalPolynomialFitAlgorithm(2, 3, QRFactorization(ColumnNorm())))

    @test sol.coeffs[1] ≈ 1.0 atol=1.0e-8
    @test sol.coeffs[2] ≈ 0.0 atol=1.0e-8
    @test sol.coeffs[3] ≈ -2.0 atol=1.0e-8
    @test sol.coeffs[4] ≈ 2.0 atol=1.0e-8
    @test sol.coeffs[5] ≈ 3.0 atol=1.0e-8
    @test sol.coeffs[6] ≈ 0.0 atol=1.0e-8

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ r(val) atol=1.0e-8
    end
end

# @testitem "Nonlinear Rational fit" begin
#     x = range(1, stop = 10, length = 10)
#     r = RationalPoly([1.0, 0.0, -2.0], [1.0, 2.0, 3.0])
#     y = r.(x)
#     f = rational_fit(x, y, 2, 3)

#     @test f[1]≈1.0 atol=1.0e-7
#     @test f[2]≈0.0 atol=1.0e-7
#     @test f[3]≈-2.0 atol=1.0e-7
#     @test f[4]≈2.0 atol=1.0e-7
#     @test f[5]≈3.0 atol=1.0e-7
#     @test f[6]≈0.0 atol=1.0e-7

#     f = curve_fit(RationalPoly, x, y, 2, 3)
#     @test f(1.5)≈r(1.5) atol=1.0e-8
#     @test f(4.5)≈r(4.5) atol=1.0e-8
# end
