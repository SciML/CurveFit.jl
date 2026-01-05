@testitem "Linear solver choice testing" begin
    @testset "LUFactorization" begin
        using CurveFit
        using NonlinearSolve

        X = collect(1.0:10.0)
        θ_true = [3.0, 2.0]

        function f(θ, X)
            return @. θ[1] + X^(θ[2])
        end
        Y = f(θ_true, X)

        nonf = NonlinearFunction(f)
        alg = LevenbergMarquardt(linsolve = LUFactorization())

        prob = NonlinearCurveFitProblem(nonf, [0.5, 0.1], X, Y)
        sol = solve(prob, alg)

        @test sol.u ≈ θ_true
        @test alg.descent.descent.linsolve isa LUFactorization
    end

    @testset "QRFactorization" begin
        using CurveFit
        using NonlinearSolve

        X = collect(1.0:10.0)
        θ_true = [3.0, 2.0]

        function f(θ, X)
            return @. θ[1] + X^(θ[2])
        end
        Y = f(θ_true, X)

        nonf = NonlinearFunction(f)
        alg = LevenbergMarquardt(linsolve = QRFactorization())

        prob = NonlinearCurveFitProblem(nonf, [0.5, 0.1], X, Y)
        sol = solve(prob, alg)

        @test sol.u ≈ θ_true
        @test alg.descent.descent.linsolve isa QRFactorization
    end

    @testset "CholeskyFactorization" begin
        using CurveFit
        using NonlinearSolve

        X = collect(1.0:10.0)
        θ_true = [3.0, 2.0]

        function f(θ, X)
            return @. θ[1] + X^(θ[2])
        end
        Y = f(θ_true, X)

        nonf = NonlinearFunction(f)
        alg = LevenbergMarquardt(linsolve = CholeskyFactorization())

        prob = NonlinearCurveFitProblem(nonf, [0.5, 0.1], X, Y)
        sol = solve(prob, alg)

        @test sol.u ≈ θ_true
        @test alg.descent.descent.linsolve isa CholeskyFactorization
    end
end
