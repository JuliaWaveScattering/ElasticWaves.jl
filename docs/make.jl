# The examples in the manual are run when building the docs. The line below lets Plots save figures without a display.
ENV["GKSwstype"] = "100"

using Documenter, ElasticWaves

makedocs(
    sitename = "ElasticWaves.jl",
    authors = "Art Gower, Jess Kent, and Matheus Loures",
    modules = [ElasticWaves],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true", # changes urls if built locally
        canonical = "https://JuliaWaveScattering.github.io/ElasticWaves.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "manual/bearing.md",
            "manual/scattering.md"
        ],
        "Library" => "library/library.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaWaveScattering/ElasticWaves.jl.git",
    devbranch = "main"
)
