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
        "Introduction"        => "index.md",
        "WTG Schemes"         => "schemes.md",
        "Experimental Setups" => "setups.md",
        "Comparing WTG Schemes" => [
            "Results"                      => "comparison/results.md",
            "Vertical Mode Decomposiition" => "comparison/verticalmodes.md",
            "Gross Moist Stability"        => "comparison/grossmoist.md",
        ],
        "Idealized Radiation vs RRTM" => [
            "Results"                           => "radiation/results.md",
            "Implications for Self-Aggregation" => "radiation/selfaggregation.md",
        ],
        "Convectively Coupled Waves?"=> [
            "Results"              => "ccw/results.md",
            "Diurnal vs Perpetual" => "ccw/diurnalcycle.md",
        ],
    ]
)

deploydocs(
    repo = "github.com/natgeo-wong/ExploreWTGSpace.git",
    devbranch = "main"
)
