using CurveFit
using Test

@testset "ExpSumFit: Integration rules" begin
    # Trapezoidal rule
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 1:1; m = 1) == [1 // 2 1 // 2]

    # Simpson's first (1/3) rule
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 1:1; m = 2) ==
        [1 // 3 4 // 3 1 // 3]
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 2:2; m = 2) ==
        [2 // 3 4 // 3 0 // 1]
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 3:3; m = 2) ==
        [3 // 5 4 // 5 -1 // 15]

    # Simpson's second (3/8) rule
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 1:1; m = 3) ==
        [3 // 8 9 // 8 9 // 8 3 // 8]
    @test CurveFit.__calc_integral_rules(Rational{Int64}, 2:2; m = 3) ==
        [39 // 40 27 // 10 27 // 40 3 // 20]
end
