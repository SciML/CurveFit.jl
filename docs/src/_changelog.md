```@meta
CurrentModule = CurveFit
```

# Changelog

This documents notable changes in CurveFit.jl. The format is based on [Keep a
Changelog](https://keepachangelog.com).

## v[1.9.1] - 2026-04-25

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
