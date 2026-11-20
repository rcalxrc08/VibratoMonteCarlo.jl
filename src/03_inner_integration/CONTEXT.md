# Stage 03: inner integration of the last step

The vibrato step: the payoff integrated against the density of the last step, derivatives only in the density weight.

## Inputs

| Load | Why |
|:--|:--|
| `lrm.jl` | `log_density`, `v_value`, `v_mod`, `eps_value`, `carries_derivatives`, the integrand with its baseline, `lrm_interface_vec!` per config × inner integration, `singular_points`, `lrm_paths` and the device matrices |
| `lrm_aug.jl` | the FFT density of the last step: `FFTStepLaw`, `fft_step_law`, `value_law`, `lrm_interface_aug!`, `aug_paths` (blocks of `AUG_BLOCK_NODES` on the device), `floor_density` |
| `references/inner_integration.md` | the core idea, the rules, smiles, inference |
| `references/autodiff.md` | `v_value` must strip every level; Mooncake overlays `v_mod` and `carries_derivatives` |
| `references/gpu.md` | the paths × nodes matrices |
| `references/method_limits.md` | bounds and node counts known to fail |

## Process

- Sample points from values only (`v_value(mu) + v_value(sigma) * z`); every derivative flows through `v_mod`.
- Keep the baseline: the value is unchanged, the derivative noise is cut 10-80×.
- `LRMInterfaceArguments`/`IntegrandArguments` keep `@.` broadcasting over the paths only.
- The host loop and the device matrix of a rule must agree (tested to ≤ 4e-11): change both together.
- A smile is one integral per option on shared nodes (fixed-node rules), one column per option on the device.

## Outputs

`lrm_paths(mu, sigma, args, mcBaseData, vb_mc)` and `aug_paths(law, option, vb_mc, S_end, mcBaseData)`: one undiscounted
value per path (a matrix for a smile). Checked by `test/gauss_legendre.jl`, `test/taylor.jl`, `test/regressions.jl`,
`test/quality.jl` (inference), `test/cuda/runtests.jl`.
