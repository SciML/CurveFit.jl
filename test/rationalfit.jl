@testitem "Linear Rational fit" begin
    using LinearSolve, LinearAlgebra

    x = range(1, stop = 10, length = 10)
    r = CurveFit.RationalPolynomial([1.0, 0.0, -2.0], [1.0, 2.0, 3.0])
    y = r.(x)

    prob = CurveFitProblem(x, y)
    sol = solve(prob, RationalPolynomialFitAlgorithm(2, 3, QRFactorization(ColumnNorm())))

    @test sol.u ≈ [1.0, 0.0, -2.0, 2.0, 3.0, 0.0] atol=1.0e-8

    @testset for val in (0.0, 1.5, 4.5, 10.0)
        @test sol(val)≈r(val) atol=1.0e-8
    end
end

@testitem "Nonlinear Rational fit" begin
    using NonlinearSolve

    x = range(1, stop = 10, length = 10)
    r = CurveFit.RationalPolynomial([1.0, 0.0, -2.0], [1.0, 2.0, 3.0])
    y = r.(x)

    algs = [
    # TODO: broken due to https://github.com/SciML/NonlinearSolve.jl/issues/504
    # RationalPolynomialFitAlgorithm(2, 3),
        RationalPolynomialFitAlgorithm(2, 3, LevenbergMarquardt())
    ]
    u0s = [nothing, fill(10.0, 6)]
    for alg in algs, u0 in u0s
        prob = CurveFitProblem(x, y; u0 = u0)
        sol = solve(prob, alg)

        @test sol.u ≈ [1.0, 0.0, -2.0, 2.0, 3.0, 0.0] atol=1.0e-8
        @test SciMLBase.successful_retcode(sol.retcode)

        @testset for val in (0.0, 1.5, 4.5, 10.0)
            @test sol(val)≈r(val) atol=1.0e-8
        end
    end
end
