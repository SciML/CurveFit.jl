import Changelog
using Documenter
using CurveFit
using CommonSolve
using NonlinearSolve: NonlinearSolve

include("pages.jl")

# Revise if possible to pick up any docstring changes when running under
# LiveServer.jl.
if isdefined(Main, :Revise)
    Revise.revise()
end

# Build the changelog. The actual contents are in `_changelog.md` and the built
# file is in `changelog.md`.
Changelog.generate(
    Changelog.Documenter(),
    joinpath(@__DIR__, "src/_changelog.md"),
    joinpath(@__DIR__, "src/changelog.md"),
    repo = "SciML/CurveFit.jl"
)

makedocs(;
    sitename = "CurveFit.jl",
    authors = "CurveFit Contributors",
    modules = [CurveFit, CommonSolve],
    clean = true,
    doctest = false,
    warnonly = [:missing_docs, :docs_block, :cross_references],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://docs.sciml.ai/CurveFit/stable/",
        assets = String[],
    ),
    pages = pages
)

deploydocs(;
    repo = "github.com/SciML/CurveFit.jl.git",
    devbranch = "master",
    push_preview = true,
)
