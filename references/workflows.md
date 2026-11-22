# Workflows: commands, CI, sub-environments, untracked files

_How to run things and where they run._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

## Commands

```bash
julia --project -e 'import Pkg; Pkg.instantiate()'
julia --project -e 'import Pkg; Pkg.test()'                     # runs test/runtests.jl (~12 min on Julia 1.11)
VIBRATO_TEST_ACCURACY=true julia --project -e 'import Pkg; Pkg.test()'   # also runs test/check/accuracy.jl
julia --project=bench -e 'import Pkg; Pkg.develop(path=pwd()); Pkg.instantiate(); include("bench/bench.jl")'
julia -e 'using JuliaFormatter; format(".")'                      # uses .JuliaFormatter.toml
julia --project=test/cuda -e 'import Pkg; Pkg.instantiate(); include("test/cuda/runtests.jl")'   # GPU tests, needs an NVIDIA GPU (~4 min on an MX150)
julia --project=test/cuda test/cuda/timings.jl                    # GPU timings of the parts moved by #51 (BENCH_PART=sub|jumps|fft|am)
julia --project=docs -e 'import Pkg; Pkg.instantiate()'; julia --project=docs docs/make.jl   # documentation in docs/build
julia etc/readme.jl                                               # copies examples/*.jl into README.md
```

`Manifest.toml` files are gitignored, so the first run resolves dependencies fresh. The private registry must be added
first (`references/dependencies.md`).

## CI

- **GitHub Actions** is the CI:
  - `CI.yml` runs the test suite on every push (Julia 1.10 and latest), with the private registry added through
    `julia-buildpkg`'s `localregistry` input, and uploads the coverage to Codecov (needs the `CODECOV_TOKEN` secret).
  - `CompatHelper.yml` runs nightly and opens the compat-bump PRs that make up most recent commits.
  - Every action runs on Node.js 24 (`checkout@v7`, `setup-julia@v3`, `cache@v3`, `codecov-action@v7`, the Pages
    actions; the `julia-buildpkg`, `julia-runtest` and `julia-processcoverage` v1 are composite). CompatHelper installs
    Julia with `setup-julia` instead of relying on the runner image. Keep new actions on Node 24 versions.
  - `DeployPage.yml` builds the Documenter site (`references/docs_design.md`).
- GitLab is no longer used: the `origin2` remote is legacy. Push only to `origin` (GitHub).
- The GPU tests do not run on CI (`references/testing.md`).

## Sub-environments

- `test/`, `test/cuda/`, `docs/`, `bench/` have their own `Project.toml`.
- **Compat of the package itself in the sub-environments** (`docs/` and `test/cuda/`, which load it through `[sources]`
  with a path): each has a `[compat]` entry for the package, so CompatHelper does not open PRs to add one (it has no
  option to skip a dependency). Bump these entries in the same commit as a breaking version bump, otherwise the local
  package no longer satisfies them and the environment fails to resolve.
- A FinancialMonteCarlo compat bump goes in the package, tests, `test/cuda`, docs and `bench` together
  (`references/dependencies.md`).

## Releases

Bump `version`, push to `master`, wait for CI, then register the commit on the private registry `localregistry`
(`rcalxrc08/LocalRegistry`, not General) with `LocalRegistry.register`, from Julia 1.12, where LocalRegistry.jl is
installed (`julia +1.12`), and push the registry. The private registry has no TagBot: tag `vX.Y.Z` on the registered
commit and create the GitHub release by hand, on `origin` (GitHub; never the GitLab mirror `origin2`). FinancialFFT is
in `localregistry` too; FinancialMonteCarlo and FinancialToolbox are in General. When FinancialFFT makes a breaking
release, widen its compat here and in `docs/`, `test/` and `test/cuda/`.

## Untracked and experimental files: not package code

`src/VibratoMonteCarlo.jl` does not include `src/deps/zygote_dependencies.jl` (untracked, marked `#wrong`), nor the
Zygote reverse-mode experiments with hand-written `rrule`s (`vibrato_zygote*.jl`, `simulate_zygote*.jl`,
`vibrato_path.jl`, untracked when present). `examples/plot_implied_vol.jl`, `fourier_greeks_method.md` and `pluto/` are
untracked notes. Don't wire them into the module, and don't delete them, unless asked.
