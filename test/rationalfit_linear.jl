using CurveFit
using Test
using LinearSolve, LinearAlgebra

@testset "Linear Rational fit" begin
    x = range(1, stop = 10, length = 10)
    r = CurveFit.RationalPolynomial([1.0, 0.0, -2.0], [1.0, 2.0, 3.0])
    y = r.(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, RationalPolynomialFitAlgorithm(2, 3, QRFactorization(ColumnNorm())))

    @test sol.u ≈ [1.0, 0.0, -2.0, 2.0, 3.0, 0.0] atol = 1.0e-8

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val) ≈ r(val) atol = 1.0e-8
    end

    @test sol(x) ≈ y atol = 1.0e-8
end
