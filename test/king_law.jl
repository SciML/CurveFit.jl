using CurveFit
using Test

@testset "King's Law" begin
    U = range(1, stop = 20, length = 20)
    A = 5.0
    B = 1.5
    n = 0.5
    E = sqrt.(A .+ B * U .^ n)

    fn(E) = ((E .^ 2 .- A) / B) .^ (1 ./ n)

    prob = CurveFitProblem(E, U)
    sol = solve(prob, KingCurveFitAlgorithm())

    @test sol.u[1] ≈ A
    @test sol.u[2] ≈ B

    @testset for val in range(minimum(E), stop = maximum(E), length = 10)
        @test sol(val) ≈ fn(val)
    end

    # Sigma not supported
    prob_sigma = CurveFitProblem(E, U; sigma = ones(length(U)))
    @test_throws AssertionError solve(prob_sigma, KingCurveFitAlgorithm())
end
