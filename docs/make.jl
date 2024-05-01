using Documenter, VibratoMonteCarlo, Literate
# EXAMPLE = joinpath(@__DIR__, "..", "examples", "getting_started.jl")
# OUTPUT = joinpath(@__DIR__, "src")
# Literate.markdown(EXAMPLE, OUTPUT; documenter = true)
# makedocs(format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/favicon.ico"]), sitename = "VibratoMonteCarlo", modules = [VibratoMonteCarlo], pages = ["index.md", "getting_started.md", "types.md", "stochproc.md", "parallel_vr.md", "payoffs.md", "metrics.md", "multivariate.md", "intdiffeq.md", "extends.md"])
makedocs(format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/favicon.ico"], repolink = "https://gitlab.com/rcalxrc08/VibratoMonteCarlo.jl"), repo = "https://gitlab.com/rcalxrc08/VibratoMonteCarlo.jl.git", sitename = "VibratoMonteCarlo", modules = [VibratoMonteCarlo], pages = ["index.md"])
get(ENV, "CI", nothing) == "true" ? deploydocs(repo = "https://gitlab.com/rcalxrc08/VibratoMonteCarlo.jl.git") : nothing