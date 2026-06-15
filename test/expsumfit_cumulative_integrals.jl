using CurveFit
using Test

@testset "ExpSumFit: Cumulative integrals" begin
    n = 4

    x = collect(0:0.4:10)
    y = @. 1 + sin(x)
    cumints_analytic = @. [(x + 1 - cos(x)) (1 / 2 * x^2 + x - sin(x)) (
        1 / 6 * x^3 +
            x^2 / 2 +
            cos(x) - 1
    ) (
        1 /
            24 *
            x^4 +
            x^3 /
            6 +
            sin(x) -
            x
    )]

    @testset "m = 1 & n = 4" begin
        coeffs = CurveFit.__calc_integral_rules(Float64, 1:4; m = 1)

        len = length(x)
        nY, mY = 1 + (len - 1) ÷ 1, 2 * n + 1

        Y = Matrix{Float64}(undef, nY, mY)
        S = Matrix{Float64}(undef, nY, 4 - 1)

        CurveFit.__cumulative_integrals!(Y, S, coeffs, x, y, 4, 1)

        @test cumints_analytic[:, 1] ≈ Y[:, 1] rtol = 3.0e-3
        @test cumints_analytic[:, 2] ≈ Y[:, 2] rtol = 3.0e-3
        @test cumints_analytic[:, 3] ≈ Y[:, 3] rtol = 4.0e-3
        @test cumints_analytic[:, 4] ≈ Y[:, 4] rtol = 4.0e-3
    end

    @testset "m = 2 & n = 4" begin
        coeffs = CurveFit.__calc_integral_rules(Float64, 1:4; m = 2)

        len = length(x)
        nY, mY = 1 + (len - 1) ÷ 2, 2 * n + 1

        Y = Matrix{Float64}(undef, nY, mY)
        S = Matrix{Float64}(undef, nY, 4 - 1)

        CurveFit.__cumulative_integrals!(Y, S, coeffs, x, y, 4, 2)

        @test cumints_analytic[1:2:end, 1] ≈ Y[:, 1] rtol = 3.0e-5
        @test cumints_analytic[1:2:end, 2] ≈ Y[:, 2] rtol = 4.0e-5
        @test cumints_analytic[1:2:end, 3] ≈ Y[:, 3] rtol = 5.0e-5
        @test cumints_analytic[1:2:end, 4] ≈ Y[:, 4] rtol = 6.0e-5
    end

    x = collect(range(0, stop = 2, length = 13))
    y = @. exp(x)
    cumints_analytic = @. [(exp(x) - 1) (exp(x) - 1 - x) (exp(x) - 1 - x - x^2 / 2) (
        exp(x) -
            1 - x -
            x^2 /
            2 -
            x^3 /
            6
    )]

    @testset "m = 2 & n = 4" begin
        coeffs = CurveFit.__calc_integral_rules(Float64, 1:4; m = 2)

        len = length(x)
        nY, mY = 1 + (len - 1) ÷ 2, 2 * n + 1

        Y = Matrix{Float64}(undef, nY, mY)
        S = Matrix{Float64}(undef, nY, 4 - 1)

        CurveFit.__cumulative_integrals!(Y, S, coeffs, x, y, 4, 2)

        @test cumints_analytic[1:2:end, 1] ≈ Y[:, 1] rtol = 5.0e-6
        @test cumints_analytic[1:2:end, 2] ≈ Y[:, 2] rtol = 3.0e-5
        @test cumints_analytic[1:2:end, 3] ≈ Y[:, 3] rtol = 3.0e-5
        @test cumints_analytic[1:2:end, 4] ≈ Y[:, 4] rtol = 4.0e-5
    end
end
