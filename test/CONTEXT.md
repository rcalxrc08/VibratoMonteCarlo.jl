# Tests

## Inputs

| Load | Why |
|:--|:--|
| `runtests.jl` | the explicit list of included files |
| `utils.jl` | backends loaded, model constructors, `test_value_and_derivatives` |
| `context.jl` | the context layers follow the code: run it after moving, adding or renaming a file (seconds) |
| `references/testing.md` | what each file covers, seeds and tolerances, durations |
| `references/method_limits.md` | settings known to miss a tolerance |
| `references/workflows.md` | commands, CI, the GPU project |

## Process

Add a test to the file of the stage it checks; a new file goes into `runtests.jl` and guards on `TEST_UTILS_LOADED`.
Calibrate a tolerance on the spread over several seeds, not on one run. GPU tests go to `cuda/runtests.jl` (opt-in, own
project).

## Outputs

`julia --project -e 'import Pkg; Pkg.test()'`: about 12 minutes on Julia 1.11; `cuda/runtests.jl` about 4 minutes on an
MX150.
