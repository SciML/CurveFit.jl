# Add Documenter.jl Documentation for CurveFit.jl

## Summary
This PR adds comprehensive documentation for CurveFit.jl using Documenter.jl, following the documentation structure and patterns used by NonlinearSolve.jl and other SciML packages.

## Changes
- Set up Documenter.jl infrastructure with `docs/make.jl` and `docs/pages.jl`
- Added comprehensive documentation structure including:
  - Getting started tutorial
  - Advanced nonlinear fitting tutorial  
  - API reference documentation
  - Basics section covering problems, solutions, and FAQ
  - References with citations
- Created sample source files with docstrings (`types.jl`, `linear.jl`)
- Added GitHub Actions workflow for automatic documentation deployment
- Included bibliography with key references for curve fitting algorithms

## Documentation Structure
```
docs/
├── make.jl                 # Main documentation build script
├── pages.jl                # Page organization
├── Project.toml            # Documentation dependencies
└── src/
    ├── index.md            # Home page
    ├── references.md       # Bibliography
    ├── refs.bib           # BibTeX references
    ├── api/
    │   └── api.md         # API reference
    ├── basics/
    │   ├── curve_fitting_problems.md
    │   ├── fitting_functions.md
    │   ├── solutions.md
    │   └── faq.md
    └── tutorials/
        ├── getting_started.md
        └── nonlinear_fitting.md
```

## Key Features
- **Comprehensive tutorials**: Step-by-step guides for linear, polynomial, exponential, and custom nonlinear fitting
- **Detailed API documentation**: Full reference for all exported functions and types
- **Mathematical background**: Explains least squares formulation and fitting theory
- **Practical examples**: Real-world use cases including robust fitting, constrained optimization, and multi-response fitting
- **Troubleshooting guide**: FAQ section with solutions to common problems
- **Citation support**: Uses DocumenterCitations.jl for proper academic references

## Documentation Preview
The documentation covers:
1. Basic usage and getting started
2. Different types of curve fitting problems
3. Advanced techniques for nonlinear fitting
4. Model selection and validation
5. Handling errors and uncertainties
6. Performance optimization tips

## Next Steps
After merging, the documentation will be automatically built and deployed via GitHub Actions when:
- Pushing to main/master branch
- Creating new tags
- Opening pull requests (for preview)

## Testing
To build the documentation locally:
```bash
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```

The built documentation will be in `docs/build/`.

---

🤖 Generated with [Claude Code](https://claude.ai/code)

Co-Authored-By: Claude <noreply@anthropic.com>