using Documenter

makedocs(;
    doctest  = false,
    format   = Documenter.HTML(;
        prettyurls=get(ENV,"CI","false") == "true",
        canonical="https://natgeo-wong.github.io/ExploreWTGSpace",
        assets=String[],
    ),
    authors  = "Nathanael Wong <natgeo.wong@outlook.com>",
    sitename = "ExploreWTGSpace",
    pages    = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/natgeo-wong/ExploreWTGSpace.git",
    devbranch = "main"
)
