using Documenter, DocumenterCitations, VibratoMonteCarlo, FinancialMonteCarlo
# Documenter takes the logo and the favicon from the assets: etc/ holds the only copies, and git keeps no empty folder.
assets = mkpath(joinpath(@__DIR__, "src", "assets"))
for asset in ("logo.png", "favicon.ico")
    cp(joinpath(@__DIR__, "..", "etc", asset), joinpath(assets, asset), force = true)
end
# The figures are computed with the package, in a light and a dark variant (see figures.jl and assets/figures.css), and
# the derived formulas are generated with Symbolics (formulas.jl, called by the @eval blocks of the pages): each in a
# module of its own, so that Plots and Symbolics do not clash with the names of the packages in Main.
module DocumentationFigures
include(joinpath(@__DIR__, "figures.jl"))
end
module DocumentationFormulas
include(joinpath(@__DIR__, "formulas.jl"))
end
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)
makedocs(
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/favicon.ico", "assets/figures.css", "assets/style.css", "assets/citations.css"], edit_link = "master", size_threshold = 400 * 2^10),
    repo = Remotes.GitHub("rcalxrc08", "VibratoMonteCarlo.jl"),
    sitename = "VibratoMonteCarlo.jl",
    modules = [VibratoMonteCarlo],
    checkdocs = :exports,
    plugins = [bib],
    pages = ["index.md", "Theory" => ["lrm.md", "vibrato.md", "refinements.md", "american.md"], "methods.md", "ad.md", "api.md", "references.md"],
)
