using Documenter
using CurveFit

makedocs(;
    sitename = "CurveFit.jl",
    authors = "CurveFit Contributors",
    modules = [CurveFit],
    clean = true,
    doctest = false,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://docs.sciml.ai/CurveFit/stable/",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/SciML/CurveFit.jl.git",
    devbranch = "master",
    push_preview = true,
)