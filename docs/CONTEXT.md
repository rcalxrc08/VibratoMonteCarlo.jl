# Documentation

## Inputs

| Load | Why |
|:--|:--|
| `make.jl`, `figures.jl`, `formulas.jl` | build, figures computed at build time, formulas generated with Symbolics |
| `src/*.md` | the pages (Theory, Guides, Reference) |
| `references/docs_design.md` | the design shared with the other packages, pages, formulas, figures |

## Process

Formulas are generated, never typed by hand; figures come in a light and a dark variant; the shared CSS stays identical
across repositories; the original-research attributions stay.

## Outputs

`julia --project=docs docs/make.jl` builds `docs/build`; `DeployPage.yml` publishes it.
