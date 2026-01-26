@testitem "King's Law" begin
    U = range(1, stop = 20, length = 20)
    A = 5.0
    B = 1.5
    n = 0.5
    E = sqrt.(A .+ B * U .^ n)

    fn(E) = ((E .^ 2 .- A) / B) .^ (1 ./ n)

    prob = CurveFitProblem(E, U)
    sol = solve(prob, KingCurveFitAlgorithm())

    @testset for val in range(minimum(E), stop = maximum(E), length = 10)
        @test sol(val) ≈ fn(val)
    end

    # Sigma not supported
    prob_sigma = CurveFitProblem(E, U; sigma = ones(length(U)))
    @test_throws AssertionError solve(prob_sigma, KingCurveFitAlgorithm())
end

@testitem "Modified King's Law" begin
    using NonlinearSolve

    U = range(1, stop = 20, length = 20)
    A = 5.0
    B = 1.5
    n = 0.42
    E = sqrt.(A .+ B * U .^ n)

    fn(E) = ((E .^ 2 .- A) ./ B) .^ (1 ./ n)

    algs = [
        # TODO: broken due to https://github.com/SciML/NonlinearSolve.jl/issues/504
        # nothing,
        LevenbergMarquardt(),
    ]
    u0s = [nothing, fill(1.0, 3)]
    for alg in algs, u0 in u0s
        prob = CurveFitProblem(E, U; u0 = u0)
        sol = solve(prob, ModifiedKingCurveFitAlgorithm(alg))

        @test sol.u ≈ [A, B, n] atol = 1.0e-8
        @test SciMLBase.successful_retcode(sol.retcode)

        @testset for val in range(minimum(E), stop = maximum(E), length = 10)
            @test sol(val) ≈ fn(val) atol = 1.0e-8
        end
    end
end
