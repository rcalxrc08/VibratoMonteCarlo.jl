# Status

_What is in flight. Layer 4: volatile, rewritten freely; GitHub issues are the source of truth, this file is the index
an agent reads first._ Updated: 2026-11-21

## Current focus

Separate the device-only code from the host code: files, names, contracts (#75), now that #74 (Standard inner
integration allocations) is closed with the weights fix 698ddbd.

## Active work

| Item | Branch | Issue | Next step |
|:--|:--|:--|:--|
| Device-only code split | `master` | #75 | list the device-only functions per stage, then move them with their contracts |

## Known issues

| Issue | Symptom | Workaround |
|:--|:--|:--|
| #70 | `log_gamma` promotes Float32 to Float64 | — |
| #67 | VibratoSaltando on VG/NIG and American options price Float32 models in Float64 | — |
| #68 | Vibrato/PureLRM drift `r - d` assumes scalar rates | scalar rates only; unblocked: FinancialMonteCarlo#81 (lazy `r - d`, no arithmetic between curves) shipped in 0.9.0, so fix with `integral(r, T) - integral(d, T)` and a curve test |
| #54 | FFT density path: binaries converge only O(dx) | finer grid |
| #34 | Kou `vibrato_saltando` misses the 1e-3 tolerance in `check/accuracy.jl` | — |
| #28 | FFT vibrato cannot differentiate κ of VG and NIG | — |

## Technical debt and enhancements

- #71 ChainRulesCore is a stale hard dependency.
- #62 remaining multiply-adds without `muladd`; #61 Heston antithetic Euler step (reuse FinancialMonteCarlo kernels).
- #57 American bisection → `FinancialFFT.brent_root`; #44 reuse FinancialMonteCarlo's jumps from uniforms.
- #53 GPU: antithetic G_T of VG/NIG drawn on the host; #59 exact Heston last step; #60 Brownian bridge for barriers.
- #33 generic `ItoProcess` methods; #32 λ and T derivatives for Kou/Merton; #38 docs: dual spectral theorem.

## Working notes

Notes of an active initiative go next to this file as `state/YYYY-MM-DD-<topic>.md`; when it ends, what stays true
moves to `references/` or `references/decisions.md` and the note is deleted. Untracked today (none is included by the module):
`fourier_greeks_method.md`, `examples/plot_implied_vol.jl`, `pluto/`, `src/deps/`.
