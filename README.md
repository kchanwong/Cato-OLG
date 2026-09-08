# Social Security OLG Model (Julia)

This is an overlapping generations (OLG) model of the United States economy. The
households in the model are married couples. Use the model to score Social
Security reforms. To score a reform is to measure its effect on Social Security
finances. 

## The two parts of the model

The model has two parts. 

**1. The steady state.** The function `solve_steady_state` calculates a general
equilibrium. Households choose how much to consume, save, and work. Wages and
interest rates adjust until the capital and labor markets are in balance.
Changes in household behavior occur in this part of the model. This part is the
slower of the two. It takes about 1 minute with 6 threads.

**2. The projection.** The function `project_economy` takes one or two steady
states. It then calculates a 75 year path of Social Security finances in nominal
dollars. The path includes taxable payroll, outlays, the trust fund, and the
actuarial balance. Productivity, inflation, population, and the dependency ratio
are exogenous inputs to this part. The score comes from this part.

To score a reform, solve a second steady state with the reform policy. Then give
both steady states to `project_economy`. The projection moves from baseline
behavior to reform behavior after the year `reform_year`.

## Setup

Install the packages:

```julia
julia> import Pkg
julia> Pkg.add(["Distributions", "Optim", "JLD2", "DataFrames", "XLSX", "CSV", "Plots"])
```

Start Julia with more than one thread. The steady state solver calculates the 16
couple types in parallel.

```
$ julia --threads=auto
```

Load the model:

```julia
julia> include("source/setup.jl")
```

The file `setup.jl` sets the variable `REPO` to the path of the repository root.
Use `REPO` to make output paths. The file also loads two objects from
`cached_objects.jld2`:

| Name | Description |
|---|---|
| `ss` | The baseline steady state under current law. A new calculation takes about 1 minute. It is in the cache so that every analysis starts from the same baseline. |
| `dep_fit` | A path of the ratio of retirees to workers, with one value per projection year. The path is fitted so that the baseline projection agrees with the SSA Trustees Report. |

Always use the cached `ss` and `dep_fit` for the baseline. Then each analysis
starts from the same numbers. To calculate them again is a separate task. Refer
to [docs/EXTENDING.md](docs/EXTENDING.md).

## How to run an analysis

Copy the file [`analysis/TEMPLATE.jl`](analysis/TEMPLATE.jl) into
`analysis/new_policies/`. Then change the copy. Each analysis has the same four
steps:

```julia
include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))

# Step 1. Define the reform as a change to create_params().
par_reform = create_params()
par_reform[:J_retire] = 70                    # Raise the full retirement age.

# Step 2. Solve the steady state of the reform. This takes about 1 minute.
ss_reform = solve_steady_state(par_reform)

# Step 3. Calculate the baseline projection.
proj_base = project_economy(ss,
    n_years = 75, start_year = 2025,
    g_A = 0.0114, g_pop = 0.005, inflation = 0.024,
    dep_path = dep_fit, ss_cola = "wage",
    trust_fund_init = 2.8e12, trust_fund_rate = 0.047,
    gdp_anchor = 28e12, label = "Current Law")

# Step 4. Calculate the reform projection. Use the same assumptions. Add
# ss_reform, reform_year, and baseline_proj. Then the score table prints.
proj_reform = project_economy(ss,
    ss_reform = ss_reform, reform_year = 2034,
    n_years = 75, start_year = 2025,
    g_A = 0.0114, g_pop = 0.005, inflation = 0.024,
    dep_path = dep_fit, ss_cola = "wage",
    trust_fund_init = 2.8e12, trust_fund_rate = 0.047,
    gdp_anchor = 28e12, baseline_proj = proj_base, label = "FRA 70")
```

**Caution:** The baseline projection and the reform projection must have the same
assumptions. The score is a difference between the two projections. If `g_A` or
`trust_fund_init` is different, the difference looks like an effect of the
policy.

Both functions return a NamedTuple. It contains 75 element vectors, such as
`ss_cash_flow_nom`, `trust_fund_nom`, and `ss_cost_rate`. It also contains
single values, such as `actuarial_balance_pct`. For the
full list, refer to
[docs/PARAMETERS.md](docs/PARAMETERS.md#3-projection-outputs).

The scripts write spreadsheets to the `output/` directory. The `setup.jl` file
makes this directory. Git does not track it.

## Static scores and dynamic scores

The function `static_score_economy` holds the macroeconomic aggregates at the
baseline steady state. The policy changes, but household behavior does not
change. Run both functions to find the size of the general equilibrium effect.

Almost all keywords are the same in the two functions. There are five
exceptions. Refer to
[docs/PARAMETERS.md](docs/PARAMETERS.md#2-projection-keywords).

The function `compare_projections(ss, ss_reform, reform_year)` calculates the
baseline projection, the static projection, and the dynamic projection. It then
prints the three results together.

## Run in a container

The file `run_job.jl` scores one reform without a script. It takes one JSON
object. The object can set any model parameter and any projection keyword.

```
julia --threads=auto run_job.jl '<json>'     # A JSON string,
julia --threads=auto run_job.jl job.json     # or a file.
```

The JSON object has these fields. All of them are optional.

| Field | Default | Description |
|---|---|---|
| `label` | `"Reform"` | The name of the run in the printed output. |
| `reform_year` | | The year in which the reform starts. A reform needs it. |
| `params` | `{}` | The keys of `create_params()`. Refer to [docs/PARAMETERS.md](docs/PARAMETERS.md#1-model-parameters). The job stops if a key is not a parameter, so a spelling mistake is not silent. Give `"Inf"` for a parameter with no limit. |
| `benefit_fn` | `"benefits_current_law"` | A benefit rule, by name. |
| `aime_modifier` | `null` | `"caregiver_modifier"`, or nothing. |
| `projection` | `{}` | The keywords of `project_economy` and `static_score_economy`. Refer to [docs/PARAMETERS.md](docs/PARAMETERS.md#2-projection-keywords). `dep_path` defaults to the cached `dep_fit`. The job gives each keyword to the function that accepts it. |
| `static` | `true` | Add a static score. |
| `output_prefix` | `"score"` | The name of the two output files. |

**Caution:** a reform that changes the number of retirees also needs a different
`dep_path`. Give the new path as a vector of numbers in `projection`. The file
[`analysis/new_policies/raise_retirement_age.jl`](analysis/new_policies/raise_retirement_age.jl)
shows how to scale the cached path.

The job always calculates the baseline projection. With no `params`, no
`benefit_fn` and no `aime_modifier`, it calculates the baseline only. In every
other case it solves the reform steady state, then it calculates the dynamic
score and the static score against the same baseline.

The job writes `<output_prefix>.xlsx`, with one sheet for each projection, and
`<output_prefix>.json`, with the headline results. Set `OUTPUT_DIR` to choose
the local directory; it defaults to `/data`.

### Build the image

```
docker build -t olg-model .
```

The image holds `cached_objects.jld2`, so each run starts from the same
baseline. The build runs `run_job.jl --selftest`, which checks that the JSON
becomes parameters of the correct type.

### Run it on this machine

The file `docker-compose.yml` mounts `./output` on `/data`, so the results of a
job appear in the repository.

```
docker compose run --rm model '{"label":"FRA 70","reform_year":2034,
    "output_prefix":"fra70","params":{"J_retire":70}}'
```

With no argument, the service runs the baseline. The container writes as root,
so the files in `output/` belong to root.

## Documentation

- [docs/PARAMETERS.md](docs/PARAMETERS.md) gives the meaning and the effect of
  each parameter and each projection keyword. It also lists the parameters that
  do nothing.
- [docs/EXTENDING.md](docs/EXTENDING.md) tells you when a reform is only a
  parameter change, when it needs a new benefit function, and when it needs a
  change to the model code.

## Files

```
source/setup.jl                     Load this file from each analysis script.
source/functions_pwbm_w_spouse.jl   All of the model code.
cached_objects.jld2                 The cached baseline steady state and dep_fit.
analysis/TEMPLATE.jl                The start point for a new analysis.
analysis/new_policies/              Scores of reforms.
analysis/counterfactuals/           Economic scenarios with no policy change.
output/                             Spreadsheets that the scripts write.
run_job.jl                          One reform from one JSON object.
Dockerfile                          The container image.
```
