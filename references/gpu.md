# GPU: CudaMode

_What runs on the device, what stays on the host, timings._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

`FinancialMonteCarlo.CudaMode`, with `using CUDA`: no extension is needed, `AbstractCudaMode` exists without CUDA, the
arrays move with `on_device` (`FinancialMonteCarlo.to_mode`) and the normals come from `normals!`
(`FinancialMonteCarlo.draw_normals!`: CURAND for a `CuArray`), and every transfer is a copy of a whole array
(`to_host`, `on_device`), so no GPU array is indexed by scalars. The helpers and `CudaConfig` are in
`src/01_methods/utils.jl`.

## Inner integration

The metrics call `path_mean(lrm_paths(mu, sigma, args, mcBaseData, vb_mc))`. `lrm_paths` broadcasts `lrm_interface_vec!`
per path; for a `CudaConfig` it evaluates the integrand on a matrix of paths × nodes in one broadcast and reduces along
the nodes (`device_samples`, `device_option_paths`), since a kernel can run neither the loops of `lrm_interface_vec!` nor
take the configuration (CPU rng). A smile is a matrix with one column per option (`path_mean`, `subordinator_mean`).
`VibratoMonteCarloStandard` draws its inner samples with CURAND.

## What runs where

- On the device: the paths of FinancialMonteCarlo, `simulate_with_vol` (Heston, same code as the CPU with `normals!`),
  the integration, `G_T` of plain MC (FinancialMonteCarlo's `draw_subordinator!`: a kernel for the Gamma of VG and the
  inverse Gaussian of NIG, seeded from the host generator, so `set_seed!` reproduces it) and its weight on `κ` and `T`
  (one broadcast, the process is isbits). The last jumps of Kou and Merton: the `JumpList` is drawn on the host by
  FinancialMonteCarlo (as its own GPU paths do), the paths stay on the device and the diffusion at the last jump is
  gathered there.
- On the host, then uploaded: the antithetic `G_T` (quantiles: the device sampler rejects, it cannot be paired, #53;
  `SobolMC` is refused in `CudaMode` by FinancialMonteCarlo).
- The FFT density of the last step (`Vibrato` and `PureLRM` on non-Ito processes): `fft_step_law` computes it with
  `FFTDensityInverter` in the mode of the configuration (it stays on the device), and `aug_paths` evaluates it with
  `FinancialFFT.density_values` on a matrix of paths × nodes, reduced along the nodes as `lrm_paths` does (a smile is
  one column per option). `density_values` builds about ten arrays of the size of that matrix, so the paths go in blocks
  of `AUG_BLOCK_NODES` (2^20 nodes, 32 MB per array with hyper duals). `PureLRM` passes `S0` as a one element device
  vector. The two FFTs differ by about 2e-15, which moves the grid by 1e-14: device and host prices agree to ~1e-10, and
  to 1e-12 on the same law. The floor of the density at eps is `floor_density`, an `ifelse` in the type of the density:
  `max(::Hyper, ::Float64)` returns either argument unconverted (`Hyper` is not `<: Real`), which is not inferred and
  fails a GPU broadcast (#52).
- American `Vibrato`: the Longstaff–Schwartz regressions and the exercise regions run on the host, on values (the
  `v_value`s of the paths, copied once); the rest on the device: the smoothed payoffs of the frozen method
  (`step_payoffs` of `FrozenFFT` / `FrozenHeston`, for a set of paths each at its column, gathered by linear index with
  `gather_at`), values only for the regressions and with derivatives for the cash flows, and the bridge vibrato of
  Black–Scholes (`bridge_step`, one exercise date for all the paths: each interval of the region, the same for every
  path, clamped to ±`BRIDGE_Z` standard deviations, and the gaps between them, so every path has the same pieces, an
  empty one of zero width; a path with no exercise piece keeps `W`, as `bridge_vibrato`).
- `PureLRM` on Ito processes: `mu` as a one element device vector through `lrm_paths`.
- So the host only draws the jump times and sizes and the antithetic `G_T`, and runs the regressions.
- Mind `@.` around `size(M, 1)`: it broadcasts `size` over `M` (use `$(size(M, 1))`).
- AD numbers must be isbits on the device: hyper duals and ForwardDiff (nested too) are tested, DualNumbers is isbits
  but untested, TaylorSeries and ReverseDiff cannot run there.

## Timings

On the MX150 (`test/cuda/timings.jl`, hyper dual `S0`, best of 3), 0.5.0 (host parts) → #51:

| part | 0.5.0 | #51 |
|:--|:--|:--|
| `G_T` of VG / NIG, 1M draws | 28 / 50 ms | 9.9 / 8.2 ms |
| last jumps, Kou plain / antithetic, 100k × 10 | 53 / 61 ms | 32 / 36 ms |
| last jumps, Merton plain / antithetic, 100k × 10 | 52 / 87 ms | 55 / 44 ms |
| `VibratoSaltando` Kou / Merton antithetic, 100k, Gauss-Legendre | 225 / 215 ms | 157 / 162 ms |
| FFT `Vibrato` Kou / NIG, 20k × 10 (the CPU: 0.85 / 0.83 s) | 1.34 / 0.88 s | 0.22 / 0.20 s |
| American put, BS / Heston, 20k × 10, Gauss-Legendre | 1.10 / 0.64 s | 0.30 / 0.13 s |
| American put, Kou, 20k × 10, FFT (the CPU: 2.4 s) | 2.24 s | 0.68 s |
| `PureLRM` Kou / NIG / BS, one point | 3.1 / 2.4 / 0.9 ms | 4.2 / 3.5 / 2.1 ms |

`VibratoSaltando` on VG and NIG stays at 0.12 s (`G_T` is a small part of it), `PureLRM` pays the kernel launches of one
point. The jumps are still drawn on the host (`add_jumps!`), most of what is left of their time.
