# PowerFlowBenchmarkPES27

Benchmarking suite for AC Optimal Power Flow (OPF) flexibility analysis. Evaluates how different solver initializations (ACOPF, DCOPF, SOC, and QC relaxations) affect the minimum number of flexible generators needed to achieve feasibility under load shedding conditions, across a large set of pglib-opf test cases.

## Requirements

- Julia 1.11
- NVIDIA GPU with CUDA support

## Installation

```bash
git clone https://github.com/mit-shin-group/power-flow-benchmark
cd PowerFlowBenchmarkPES27
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

The benchmark data (pglib-opf `.m` files) is included in the `data/` directory.

## Running the Benchmarks

Open a Julia REPL from the project root and run:

```julia
include("benchmark/benchmark_cases.jl")
benchmark_ac_vs_dc(full_cases)
```

Results are saved to `saved_data/` as CSV files and PNG plots per test case.

## Generating Summary Plots

After benchmarks complete, call `make_plots` in the same session:

```julia
make_plots("flex_cdf", "flex_vs_pshed", "cost_cdf", "flex_fraction", "cost_cdf_per_method"; cases=full_cases)
```

Output is saved to `saved_data/formatted_data/`.

## Benchmark Cases

Case lists are defined in `benchmark/benchmark_cases.jl`. The default run uses all cases. 

| Variable | Description |
|---|---|
| `short_cases` | 5 small cases for quick testing |
| `full_cases` | All 63 pglib-opf cases |
| `half_cases` | First 50 cases |
| `back_half_cases` | 8 largest cases |

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
