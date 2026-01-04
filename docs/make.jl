using Documenter
using CurveFit
using CommonSolve

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
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Tutorial" => "tutorial.md",
        "Manual" => "manual.md"
    ],
)

deploydocs(;
    repo = "github.com/SciML/CurveFit.jl.git",
    devbranch = "master",
    push_preview = true,
)
