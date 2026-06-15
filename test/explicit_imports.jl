using CurveFit
using Test
using ExplicitImports
using StatsAPI

@testset "Explicit Imports" begin
    @test check_no_implicit_imports(CurveFit) === nothing
    @test check_no_stale_explicit_imports(CurveFit) === nothing
    @test check_no_self_qualified_accesses(CurveFit) === nothing
    @test check_all_qualified_accesses_via_owners(CurveFit) === nothing
end
