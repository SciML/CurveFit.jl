@testitem "King's Law" begin
    U = range(1, stop = 20, length = 20)
    A = 5.0
    B = 1.5
    n = 0.5
    E = sqrt.(A .+ B * U .^ n)

    fn(E) = ((E .^ 2 - A) / B) .^ (1 ./ n)

    prob = CurveFitProblem(E, U)
    sol = solve(prob, KingCurveFitAlgorithm())

    @testset for val in range(minimum(E), stop = maximum(E), length = 10)
        @test sol(val) ≈ fn(val)
    end
end

# @testitem "Modified King's Law" begin
#     U = [range(1, stop = 20, length = 20);]
#     A = 5.0
#     B = 1.5
#     n = 0.42
#     E = sqrt.(A .+ B * U .^ n)
#     fun6(E) = ((E^2 - A) / B)^(1 / n)

#     f = king_fit(E, U)
#     @test f[1]≈A atol=1.0e-7
#     @test f[2]≈B atol=1.0e-7
#     @test f[3]≈n atol=1.0e-7
#     f = curve_fit(KingFit, E, U)
#     @test f(3.0)≈fun6(3.0) atol=1.0e-5
# end
