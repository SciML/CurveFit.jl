using CurveFit, BenchmarkTools
using CommonSolve: solve
using StableRNGs

const SUITE = BenchmarkGroup()
const rng = StableRNG(123)

x = collect(range(0.0, 5.0, length = 200))
y_lin = 2.0 .+ 3.0 .* x
y_pow = 2.0 .* x .^ 1.5
y_log = 1.0 .+ 2.0 .* log.(x .+ 1.0)
y_exp = exp.(-0.1 .* x) .+ 0.5 .* exp.(-1.0 .* x)

# =============================================================================
# Curve fits
# =============================================================================

SUITE["fit"] = BenchmarkGroup()

SUITE["fit"]["linear"] = @benchmarkable solve(
    $(CurveFitProblem(x, y_lin)), LinearCurveFitAlgorithm()
)
SUITE["fit"]["log"] = @benchmarkable solve(
    $(CurveFitProblem(x, y_log)), LogCurveFitAlgorithm()
)
SUITE["fit"]["power"] = @benchmarkable solve(
    $(CurveFitProblem(x, y_pow)), PowerCurveFitAlgorithm()
)
SUITE["fit"]["expsum"] = @benchmarkable solve(
    $(CurveFitProblem(x, y_exp)), ExpSumFitAlgorithm(; n = 2, m = 1)
)
SUITE["fit"]["king"] = @benchmarkable solve(
    $(CurveFitProblem(x, y_lin)), KingCurveFitAlgorithm()
)

# =============================================================================
# Construction
# =============================================================================

SUITE["construct"] = BenchmarkGroup()
SUITE["construct"]["problem"] = @benchmarkable CurveFitProblem($x, $y_lin)
