@testitem "King's Law" begin
    U = [range(1, stop = 20, length = 20);]
    A = 5.0
    B = 1.5
    n = 0.5
    E = sqrt.(A .+ B * U .^ n)
    fun5(E) = ((E .^ 2 - A) / B) .^ (1 ./ n)

    f = linear_king_fit(E, U)
    @test f[1]≈A atol=1.0e-7
    @test f[2]≈B atol=1.0e-7
    f = curve_fit(LinearKingFit, E, U)
    @test f(3.0)≈fun5(3.0) atol=1.0e-7
end

@testitem "Modified King's Law" begin
    U = [range(1, stop = 20, length = 20);]
    A = 5.0
    B = 1.5
    n = 0.42
    E = sqrt.(A .+ B * U .^ n)
    fun6(E) = ((E^2 - A) / B)^(1 / n)

    f = king_fit(E, U)
    @test f[1]≈A atol=1.0e-7
    @test f[2]≈B atol=1.0e-7
    @test f[3]≈n atol=1.0e-7
    f = curve_fit(KingFit, E, U)
    @test f(3.0)≈fun6(3.0) atol=1.0e-5
end
