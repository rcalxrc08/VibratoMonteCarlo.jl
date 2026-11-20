# Examples, README and logo

The README examples and the artwork, both copied elsewhere by a script.

## Inputs

| Load | Why |
|:--|:--|
| `readme.jl` | copies each `examples/<name>.jl` into `README.md` between `<!-- BEGIN examples/<name>.jl -->` and `<!-- END ... -->` (GitHub cannot include a file) |
| `examples/hello_world.jl` | the README example: a European call and its delta with DualNumbers, against `blsprice` |
| `logo.jl` | draws `logo.png` (Luxor and the Georgia font, in no Project.toml); `favicon.ico` is the logo padded to a square |
| `references/docs_design.md` | `docs/make.jl` copies the logo and the favicon into the docs assets |
| `references/testing.md` | `test/examples.jl` checks the README and runs the example |

## Process

- After changing an example, run `julia etc/readme.jl` from the repository root, then `test/examples.jl`.
- A new example is `examples/<name>.jl` plus the two markers in `README.md`; examples use the public API only.
- The logo and favicon live in `etc/` only (the docs copies are gitignored); redraw them with
  `julia etc/logo.jl` in an environment with Luxor.

## Outputs

`README.md` in sync with `examples/`, `logo.png`, `favicon.ico`. Checked by `test/examples.jl`.
