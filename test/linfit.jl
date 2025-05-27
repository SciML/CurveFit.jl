@testitem "Linear Fit" begin
    x = range(1, stop = 10, length = 10)

    fun0(x) = 1.0 + 2.0 * x
    y = fun0.(x)
    f = linear_fit(x, y)

    @test f[1]≈1.0 atol=1.0e-7
    @test f[2]≈2.0 atol=1.0e-7

    f = curve_fit(LinearFit, x, y)
    @test f(1.5)≈fun0(1.5) atol=1.0e-7
end

@testitem "Log Fit" begin
    x = range(1, stop = 10, length = 10)

    fun1(x) = 1.0 + 2.0 * log(x)
    y = fun1.(x)
    f = log_fit(x, y)
    @test f[1]≈1.0 atol=1.0e-7
    @test f[2]≈2.0 atol=1.0e-7
    f = curve_fit(LogFit, x, y)
    @test f(1.5)≈fun1(1.5) atol=1.0e-7
end

@testitem "Power Fit" begin
    x = range(1, stop = 10, length = 10)

    fun2(x) = 2.0 * x .^ 0.8
    y = fun2(x)
    f = power_fit(x, y)
    @test f[1]≈2.0 atol=1.0e-7
    @test f[2]≈0.8 atol=1.0e-7
    f = curve_fit(PowerFit, x, y)
    @test f(1.5)≈fun2(1.5) atol=1.0e-7
end

@testitem "Exp Fit" begin
    x = range(1, stop = 10, length = 10)

    fun3(x) = 2.0 * exp.(0.8 * x)
    y = fun3(x)
    f = exp_fit(x, y)
    @test f[1]≈2.0 atol=1.0e-7
    @test f[2]≈0.8 atol=1.0e-7
    f = curve_fit(ExpFit, x, y)
    @test f(1.5)≈fun3(1.5) atol=1.0e-7
end

@testitem "Polynomial" begin
    using Polynomials

    x = range(1, stop = 10, length = 10)

    fun4(x) = 1.0 + 2.0 * x + 3.0 * x^2 + 0.5 * x^3
    y = fun4.(x)
    f = poly_fit(x, y, 4)
    @test f[1]≈1.0 atol=1.0e-7
    @test f[2]≈2.0 atol=1.0e-7
    @test f[3]≈3.0 atol=1.0e-7
    @test f[4]≈0.5 atol=1.0e-7
    @test f[5]≈0.0 atol=1.0e-7

    f = curve_fit(Polynomial, x, y, 4)
    @test f(1.5)≈fun4(1.5) atol=1.0e-7

    # Polynomials with large numbers
    coefs = [80.0, -5e-18, -7e-20, -1e-36]
    P = Polynomial(coefs)
    x1 = 1e10 * (0:0.1:5)
    y1 = P.(x1)
    P2 = curve_fit(Polynomial, x1, y1, 3)
    c = coeffs(P2)
    @test coefs[1]≈c[1] rtol=1e-5
    @test coefs[2]≈c[2] rtol=1e-5
    @test coefs[3]≈c[3] rtol=1e-5
    @test coefs[4]≈c[4] rtol=1e-5
end
