# CONTEXT.md: task routing (Layer 1)

Find the task, open the contract it points to, and load only the inputs listed there.

| Task | Contract |
|:--|:--|
| A method type, an inner integration (nodes, bounds), input validation, execution modes | `src/01_methods/CONTEXT.md` |
| Paths with extra state: last jumps of Kou/Merton, Heston variance, the subordinator of VG/NIG | `src/02_simulation/CONTEXT.md` |
| The vibrato step: `v_value`/`v_mod`, the integrand and baseline, a rule per configuration, the FFT density of the last step | `src/03_inner_integration/CONTEXT.md` |
| A European pricer: `Vibrato`, `VibratoSaltando`, `PureLRM`, `Conditional*`, smiles | `src/04_european/CONTEXT.md` |
| American options: Longstaff–Schwartz, bridge vibrato, frozen decisions | `src/05_american/CONTEXT.md` |
| A new AD number type, Mooncake rules | `ext/CONTEXT.md` |
| Adding or running tests, seeds, tolerances, the GPU tests | `test/CONTEXT.md` |
| Documentation pages, formulas, figures | `docs/CONTEXT.md` |
| Commands, CI, compat bumps, sub-environments, untracked files | `references/workflows.md` |

Shared references (Layer 3), each loaded only when a contract lists it: `public_api.md`, `conventions.md`,
`autodiff.md`, `gpu.md`, `dependencies.md`, `workflows.md`, `testing.md`, `method_limits.md`, and one per stage:
`inner_integration.md` (01, 03), `simulation.md`, `european.md`, `american.md`; `docs_design.md` for the site.

The stages communicate through Julia values: 01 builds the method (`mc`, `vb_mc`); 02 returns the paths with the state to
condition on (last jump, variance, `G_T` and its weight); 03 returns one value per path from `lrm_paths`/`aug_paths`;
04 and 05 average, discount and return the price through `pricer`.
