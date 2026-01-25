```@meta
CurrentModule = CurveFit
```

# Changelog

This documents notable changes in CurveFit.jl. The format is based on [Keep a
Changelog](https://keepachangelog.com).

## Unreleased

### Added
- Added support for standard deviation weights for linear and nonlinear fits
  ([#79], [#80]).

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
