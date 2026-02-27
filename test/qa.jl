@testitem "Aqua tests" tags = [:qa, :nopre] begin
    using Aqua
    Aqua.test_all(CurveFit)
end

@testitem "Explicit Imports" tags = [:qa] begin
    using ExplicitImports
    using StatsAPI

    @test check_no_implicit_imports(CurveFit) === nothing
    @test check_no_stale_explicit_imports(CurveFit) === nothing
    @test check_no_self_qualified_accesses(CurveFit) === nothing
    @test check_all_qualified_accesses_via_owners(CurveFit) === nothing
end
