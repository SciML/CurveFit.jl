```@meta
CurrentModule = CurveFit
```

# Changelog

This documents notable changes in CurveFit.jl. The format is based on [Keep a
Changelog](https://keepachangelog.com).

## [v1.10.0] - 2026-06-28

### Added
- Added support for passing arrays to solutions of
  [`KingCurveFitAlgorithm`](@ref), [`ModifiedKingCurveFitAlgorithm`](@ref), and
  [`RationalPolynomialFitAlgorithm`](@ref) ([#115]).

### Fixed
- Fixed the statistics functions for [`ExpSumFitAlgorithm`](@ref) to handle
  `withconst=true` correctly ([#115]).
- Fixed `sol(x::Number)` of [`ExpSumFitAlgorithm`](@ref) to return scalars for
  consistency with the other solutions ([#115]).
- Previously the original [`CurveFitProblem`](@ref) from a nonlinear fit was
  always copied into the solution, even after calling `reinit!(cache,
  ...)`. This meant that the statistics functions like [`margin_error()`](@ref)
  etc would incorrectly return values for the original problem rather than the
  one actually solved. Now `sol.prob` is reconstructed using the correct inputs
  ([#115]).

## [v1.9.4] - 2026-06-28

### Changed
- Previously nonlinear fits would compute the residuals as `ŷ − y`, they are now
  computed as `y − ŷ` to be consistent with the linear fits ([#114]).

### Fixed
- Corrected [`margin_error()`](@ref) to use the residual degrees of freedom
  rather than the degrees of freedom of the model ([#114]).
- Fixed the covariance calculation in [`vcov()`](@ref) to correctly handle the
  uncertainties produced by linear fits of a transformed nonlinear model
  (e.g. from [`PowerCurveFitAlgorithm`](@ref), [`ExpCurveFitAlgorithm`](@ref),
  and [`KingCurveFitAlgorithm`](@ref)) by using the delta method ([#114]).

## [v1.9.3] - 2026-06-26

### Fixed
- [`LinearCurveFitAlgorithm`](@ref) will now automatically invert the intercept
  when `yfun` is given to ensure that the returned parameters match the values
  being fitted. Also affects [`ExpCurveFitAlgorithm`](@ref) and
  [`PowerCurveFitAlgorithm`](@ref). *This is considered a bugfix rather than a
  breaking change.*
- Fixed the parameter handling and Jacobian of [`KingCurveFitAlgorithm`](@ref)
  ([#112]).
- Fixed the statistics methods for [`LinearCurveFitAlgorithm`](@ref) when a
  transform is applied with `yfun` ([#112]). Previously the Jacobian for a
  linear function would be computed on the residuals stored in the original
  y-space.

## [v1.9.2] - 2026-06-24

### Changed
- Various improvements to CI.

## [v1.9.1] - 2026-04-25

### Fixed
- Fixed support for `reinit!()`'ing nonlinear fit caches when using the new
  `AutoSpecializeCallable` wrapper from NonlinearSolve ([#98]).

## [v1.9.0] - 2026-04-24

### Added
- Added an `absolute_sigma` argument to [`vcov()`](@ref) and related functions
  to control whether the covariance matrix is rescaled by reduced χ² ([#97]).
- Added a `weighted` argument to [`residuals()`](@ref) and related functions to
  control whether the returned residuals are scaled by the problem weights, if
  any ([#97]).

### Changed
- The `resid` field of [`CurveFitSolution`](@ref) now consistently stores the
  *weighted* residuals for both linear and nonlinear fits ([#97]). Previously
  the unweighted residuals were stored for linear fits.

### Fixed
- [`vcov()`](@ref) and related functions previously ignored the problem weights,
  they are now taken into account when present ([#97]).

## [v1.8.1] - 2026-04-13

### Changed
- Added support for SciMLBase v3 ([#95]).

## [v1.8.0] - 2026-04-07

### Changed
- Added support for RecursiveArrayTools 4.0 ([#94]).

## [v1.7.0] - 2026-03-02

### Added
- Added a precompilation workload to reduce TTFX ([#90]).

## [v1.6.0] - 2026-02-28

### Added
- Implemented support for bounds for some algorithms ([#87]).

## [v1.5.1] - 2026-02-16

### Changed
- Fixed compatibility with NonlinearSolveFirstOrder.jl v2 ([#86]).

## [v1.5.0] - 2026-02-14

### Changed
- CurveFit now depends only on NonlinearSolveFirstOrder.jl to reduce
  dependencies ([#85]). The default algorithm remains the same.

## [v1.4.0] - 2026-01-31

### Added
- Implemented [`margin_error()`](@ref) ([#81]).
- Added support for standard deviation weights for linear fits ([#80]).

### Changed
- [`ScalarModel()`](@ref)'s will now operate in-place for improved performance
  ([#82]).

## [v1.3.0] - 2026-01-26

### Added
- Added support for standard deviation weights for nonlinear fits ([#79]).

### Changed
- **Breaking**: `reinit!(::GenericNonlinearCurveFitCache)` now takes in `u0` as
  a keyword argument rather than a positional argument for consistency with
  NonlinearSolve.jl ([#79]).

### Fixed
- Fixed `reinit!(::GenericNonlinearCurveFitCache)` to allow passing a new
  `x`/`y` as well as `u0` ([#79]).

## [v1.2.0] - 2026-01-21

### Added
- Implemented [`ScalarModel`](@ref) to allow using scalar functions as models
  ([#75]).
- Implemented `SciMLBase.successful_retcode()` for [`CurveFitSolution`](@ref)
  ([#78]).
