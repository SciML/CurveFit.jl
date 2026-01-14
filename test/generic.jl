@testitem "CurveFitSolution" begin
    using SciMLBase

    x = 1:10
    fn(x) = 1.0 + 2.0 * x
    y = fn.(x)
    prob = CurveFitProblem(x, y)
    sol = solve(prob, LinearCurveFitAlgorithm())

    @test SciMLBase.successful_retcode(sol)
end
