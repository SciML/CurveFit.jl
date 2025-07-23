using Documenter, DocumenterCitations
using CurveFit

cp(joinpath(@__DIR__, "..", "Project.toml"),
    joinpath(@__DIR__, "src", "assets", "Project.toml"); force = true)

include("pages.jl")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    sitename = "CurveFit.jl",
    authors = "CurveFit Contributors",
    modules = [CurveFit],
    clean = true,
    doctest = false,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://docs.sciml.ai/CurveFit/stable/",
        assets = ["assets/favicon.ico"],
        mathengine = MathJax3(),
    ),
    plugins = [bib],
    pages = pages,
)

deploydocs(;
    repo = "github.com/YourOrg/CurveFit.jl.git",
    devbranch = "main",
    push_preview = true,
)