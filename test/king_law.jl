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

    @test sol(E) ≈ U

    # Sigma not supported
    prob_sigma = CurveFitProblem(E, U; sigma = ones(length(U)))
    @test_throws AssertionError solve(prob_sigma, KingCurveFitAlgorithm())
end

@testset "King vcov (delta method)" begin
    # King fits sqrt(U) = b + a·E² by OLS, then A = -b/a, B = 1/a. vcov must
    # describe that transformed-space estimator, not the U-space Jacobian.
    E = collect(1.0:0.25:3.0)
    A_true, B_true = 1.5, 0.8
    U = ((E .^ 2 .- A_true) ./ B_true) .^ 2 .+ 1.0e-3 .* sin.(1:length(E))

    sol = solve(CurveFitProblem(E, U), KingCurveFitAlgorithm())

    # Reference: plain linear OLS on (E², sqrt(U)) gives the (slope a,
    # intercept b) covariance; delta-method it through (A, B) = (-b/a, 1/a).
    sol_lin = solve(CurveFitProblem(E .^ 2, sqrt.(U)), LinearCurveFitAlgorithm())
    a, b = sol_lin.u            # (slope, intercept)
    Vlin = vcov(sol_lin)        # ordered (a, b)
    Σ = [Vlin[2, 2] Vlin[2, 1]; Vlin[1, 2] Vlin[1, 1]]  # reordered to (b, a)
    G = [-1 / a b / a^2; 0.0 -1 / a^2]                        # ∂(A, B)/∂(b, a)
    Vref = G * Σ * G'

    @test vcov(sol) ≈ Vref
end
