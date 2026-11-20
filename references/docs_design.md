# Documentation design

_The Documenter site: shared design, pages, formulas, figures._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

- `DeployPage.yml` builds the Documenter site in `docs/` (`docs/make.jl`, Julia 1, the private registry added before
  `instantiate`) and publishes it with GitHub Pages from Actions (`upload-pages-artifact` + `deploy-pages`, Pages source
  "GitHub Actions"): no `gh-pages` branch and no `DOCUMENTER_KEY`.
- The logo is drawn by `etc/logo.jl` (Luxor and the Georgia font, not in any Project.toml); `etc/favicon.ico` is the
  logo padded to a square. `docs/make.jl` copies both to `docs/src/assets` (gitignored): they live in `etc/` only.
- The design is shared by all the rcalxrc08 docs: `docs/src/assets/style.css` (Dark Quant), `dark.js`, `figures.css`,
  `site.css` (page headers with kicker, lede and section links, numbered sections, cards) and `home.css` (landing page)
  are identical copies in every repository: keep them in sync.
- The home page is that landing page; its hero is `hero_figure` at the end of `figures.jl`: the gamma of a binary call
  on the same 10^4 paths by finite differences of plain Monte Carlo (noise), by Vibrato with hyper-dual numbers, and the
  formula.
- Vibrato Saltando and the Spectral Theorem are original work of the author (Nicola Scaramuzzino): the home page
  ("Original research"), the page header of `spectral.md`, its History admonition and the "Original contribution"
  admonition of `vibrato.md` say so; keep it.
- The sidebar groups Theory, Guides and Reference. Pages: `index.md` (hero, quick start), the theory in `lrm.md`,
  `spectral.md` (the Spectral Theorem: primal in price space, its dual in frequency space with the Bell multipliers of
  `L = log E[e^{w X_T}]`, the link with FinancialFFT; and the general dual theorem on the complex plane, for every
  parameter of Merton, Kou, VG and NIG, from `fourier_greeks_method.md`; formulas `hermite_r_weights`, `dual_symbols`,
  `dual_multipliers`, `stein_bridge`, `frequency_scores`, the last one asserting its exponents against
  `FinancialFFT.characteristic_exponent_i`; the identities and the general theorem were checked numerically against
  Taylor1 through the Fourier price and `LewisMethod`, ≤ 1e-9 up to the third order; implementing the multipliers in
  FinancialFFT is the next step), `vibrato.md`, `refinements.md`, `american.md`, then `methods.md` (usage, inner
  integrations, limits), `ad.md` (backends, GPU), `api.md` (`@docs`), `references.md`.
- Derived formulas are never typed by hand: `docs/formulas.jl` generates them with Symbolics at build time from
  `VibratoMonteCarlo.log_density` (`latex_weight`, `generate_faa_di_bruno_exp`), and each is shown by an `@eval` block
  calling `Main.DocumentationFormulas.<name>()`, which returns Markdown with display math. Simplify with `δ = √h` as a
  plain symbol, since Symbolics does not cancel `sqrt(T)`. Definitions, model inputs and general identities stay in the
  text. Check new formulas numerically before publishing them.
- Figures: `docs/figures.jl` computes them with the package at build time, in a light and a dark variant (`-light.svg`,
  `-dark.svg`, shown by theme through `figures.css`), with the palette and markers of FinancialFFT's figures. They go to
  `docs/src/assets/figures` (gitignored).
- Math is ``` ``x`` ``` inline and ```` ```math ```` blocks. Citations are `[key](@cite)` (DocumenterCitations,
  author-year) with the entries in `docs/src/refs.bib`.
- To build locally: `julia --project=docs docs/make.jl` (the `[sources]` of `docs/Project.toml` points at the package;
  several minutes, most of it Symbolics and the figures), then serve `docs/build` with a static server.
