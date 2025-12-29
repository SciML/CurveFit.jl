@testitem "Exposing LevenbergMarquardt linear solve choices" begin
    using CurveFit
    using NonlinearSolve
    
    X = collect(1.0:10.0)
    θ_ref = [-1.0, 3.0, 0.5]

    function f(θ, X)
        return @. θ[1] + θ[2] * X + θ[3] * X^2
    end
    Y = f(θ_ref, X)
    
    nf = NonlinearFunction(f)

    prob = NonlinearCurveFitProblem(nf, [0.0, 0.0, 0.0], X, Y)

    @testset "LM_QR" begin
        sol_qr = solve(prob, LM_QR())
        @test sol_qr.alg.alg.trustregion isa NonlinearSolveFirstOrder.LevenbergMarquardtTrustRegion
        @test sol_qr.alg.alg.descent.descent.linsolve isa LinearSolve.QRFactorization
        @test sol_qr.u ≈ θ_ref
    end

    @testset "LM_CH" begin
        sol_ch = solve(prob, LM_CH())
        @test sol_ch.alg.alg.trustregion isa NonlinearSolveFirstOrder.LevenbergMarquardtTrustRegion
        @test sol_ch.alg.alg.descent.descent.linsolve isa LinearSolve.CholeskyFactorization
        @test sol_ch.u ≈ θ_ref
    end
end