using Documenter
using DocumenterCitations
using Literate

bib_filepath = joinpath(dirname(@__FILE__), "references.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

makedocs(;
    doctest  = false,
    format   = Documenter.HTML(collapselevel=1,prettyurls=false),
    authors  = "Nathanael Wong <natgeo.wong@outlook.com>",
    sitename = "ExploreWTGSpace",
    pages    = [
        "Home"                        => "index.md",
    ]
)

deploydocs(
    repo = "github.com/natgeo-wong/ExploreWTGSpace.git",
    devbranch = "main"
)
