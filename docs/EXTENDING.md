# How to add a reform to the model

A reform can be at one of four levels. This document gives the levels in order
of cost, from the lowest cost to the highest. Find the level of your reform
before you start work. The amount of work at each level is very different.

## Level 1. A projection keyword

At this level you do not solve a steady state. The run takes a few seconds.

Use this level for two types of work. The first type is a reform that changes
the accounting of Social Security but does not change how much people work or
save. The second type is a change to an economic assumption, not to a policy.

These are some examples:

- A different measure of inflation for the COLA (`cola_inflation`).
- The COLA cap (`cola_cap`).
- Price indexing of the PIA factors inside the projection
  (`pia_factor_indexing`).
- A different interest rate on the trust fund (`trust_fund_rate`).
- A scenario with high inflation and low growth (`scenario_periods`).

Give the baseline projection and the reform projection the same macroeconomic
assumptions. Then give the reform run `baseline_proj = proj_base`, so that the
score table prints. For an example, refer to
`analysis/counterfactuals/inflation.jl`.

## Level 2. A parameter change

At this level you solve one new steady state. This takes about 1 minute with 6
threads. Use this level if the reform changes a number that is already in the
model.

```julia
par_reform = create_params()
par_reform[:payroll_rate] = 0.124        # Increase the payroll tax rate.
par_reform[:payroll_cap]  = Inf          # Remove the taxable maximum.
par_reform[:J_retire]     = 70           # Increase the full retirement age.
par_reform[:n_years]      = 40           # Make the AIME period longer.
ss_reform = solve_steady_state(par_reform)
```

**Caution:** `par` is a plain dictionary. If you spell a key incorrectly, the
model ignores it and the reform has no effect. Check each key in
[PARAMETERS.md](PARAMETERS.md). Then confirm the result:

```julia
extract_ss_ratios(ss).ss_ben_per_retiree, extract_ss_ratios(ss_reform).ss_ben_per_retiree
```

Two reforms look like a parameter change, but they need more work:

**An increase in the retirement age also changes the number of retirees.** The
key `par[:J_retire] = 70` gives the effect on households. They work longer and
they get a higher AIME. But the projection counts retirees from `dep_path`. You
must scale that path yourself. For an example, refer to the variable
`dep_path_le` in `analysis/new_policies/raise_retirement_age.jl`.

**Removal of the taxable maximum needs an extra keyword in a static score.**
Give `static_score_economy` the keyword `eliminates_cap = true` in addition to
`par[:payroll_cap] = Inf`. A static score has no change in behavior. Therefore
it cannot calculate the new share of earnings above the maximum.

## Level 3. A new benefit rule

At this level you solve one new steady state, but you do not change the model.
Use this level if the reform changes how the model calculates benefits from a
history of earnings. This level is the primary extension point of the model.
Most reforms are at this level.

A benefit function takes the two AIME values of a couple. It returns the two
benefits. It must give a result for each of the three marital states.

```julia
function benefits_my_reform(AIME1::Float64, AIME2::Float64,
                            par::Dict, status::Symbol)
    PIA1 = compute_PIA(AIME1, par)
    PIA2 = compute_PIA(AIME2, par)
    if status == :couple
        return (PIA1, PIA2)                       # No spousal supplement.
    elseif status == :survivor1                   # Spouse 1 is alive.
        return (max(PIA1, PIA2), 0.0)
    elseif status == :survivor2                   # Spouse 2 is alive.
        return (0.0, max(PIA2, PIA1))
    end
    (PIA1, PIA2)
end

par = create_params()
ss_reform = solve_steady_state(par; benefit_fn = benefits_my_reform)
```

**Caution:** give the rule to `solve_steady_state`. If you give it only to
`create_params`, the model ignores it and scores current law without an error.
Refer to the known problems below.

The function must obey four rules:

- Return a tuple of two values.
- Put spouse 1 first and spouse 2 second, in all states.
- Return `0.0` for the spouse who is dead in a survivor state.
- Return values in model units, not in dollars.

The model calls this function many times. The household solver calls it once
for each combination of the two productivity states and the age. The function
`compute_aggregates_couples` calls it inside its innermost loop, which is the
loop over assets. Therefore keep the function fast, and do not allocate memory
in it.

The file `source/functions_pwbm_w_spouse.jl` contains five rules that you can
use as examples: `benefits_current_law`, `benefits_capped_spousal`,
`benefits_earnings_sharing`, `benefits_caregiver_credits`, and
`benefits_price_index_floor`. The analysis scripts contain two more rules: a
phase out of the spousal benefit for the top quintile, and flat benefits.

Some reforms change creditable earnings and not the benefit formula. For these
reforms, use `aime_modifier` instead. This function changes the two histories of
earnings before the model averages them into AIME.

```julia
par = create_params()
ss_reform = solve_steady_state(par;
                               benefit_fn    = benefits_caregiver_credits,
                               aime_modifier = caregiver_modifier)
```

The modifier receives five arguments: `(e1, e2, n1_path, n2_path, par)`. The
first two arguments are the indexed earnings of each spouse. The next two are
the hours of work of each spouse. Each of the four is a vector over the working
life. The modifier returns a new `(e1, e2)`.

Credits for care of a family member need both functions. The modifier gives the
credit. The benefit rule then removes the spousal supplement that the credit
replaces.

Give your rule a name that is not the same as any other name. The function
`solve_steady_state` prints the name of the rule that it uses. Read that line
before you use a result.

## Level 4. A change to the model code

All other reforms are at this level. This work takes days. You must change
several functions that depend on each other.

You are at this level if the reform changes the state space or the choices of
the household.

| Reform | Reason |
|---|---|
| A choice of when to claim benefits, with an actuarial adjustment | The age at which a household claims benefits is a parameter, not a choice. The model needs a new choice and a new state. |
| Private accounts, with a portfolio choice | The file `analysis/new_policies/privatisation.jl` moves FICA revenue inside the projection only. Households do not hold the accounts and do not see them. |
| Divorce, remarriage, or households that never marry | A change of marital status is one way only, from `couple` to `survivor`. There is no type for a household that never marries. |
| Fertility or immigration as a choice | Population is exogenous. The keyword `g_pop` and the vector `dep_path` set it. |
| A benefit test on income or on assets | A benefit function receives AIME only. It does not receive assets. You must pass the asset state to it. |
| An open economy | The parameters for the foreign sector exist, but no function reads them. |
| Different mortality | The arrays `SSA_Q_MALE` and `SSA_Q_FEMALE` hold the life tables. To replace the arrays is easy. To make mortality depend on income is not easy. |
| Disability insurance | The model covers OASI only. |

At this level, four functions do most of the work:

- `run_egm_couple` and `run_egm_survivor` solve the household problem with the
  endogenous grid method.
- `compute_distribution_couples` calculates the population forward in time.
- `compute_aggregates_couples` adds up the results.

The state space of a couple is `(assets, z1, z2, age)`. The state space of a
survivor is `(assets, z, age)`. Each of the functions loops over one of these
two state spaces, or over both. Therefore, if you add a state, you must change
all four functions in the same way. You must also change
`payroll_stats_couples`.

---

## How to build the cache again

The file `cached_objects.jld2` holds two objects. The first is `ss`, the
baseline steady state under current law. The second is `dep_fit`, the fitted
path of the dependency ratio.

Build the cache again only if the baseline changes. There are four reasons to do
this:

- A change to the calibration, such as `target_avg_earn`, `beta`, or `gamma`.
- A change to a grid.
- A correction to current law in the model.
- A new Trustees Report to fit against.

A reform is never a reason to build the cache again. Each reform gets its own
steady state.

```julia
include(joinpath("source", "setup.jl"))     # This gives you the model code.

ss = solve_steady_state(create_params())    # About 1 minute with --threads=auto.

# Fit the dependency path, so that the baseline projection agrees with the
# Trustees Report. The dictionary `target` maps a calendar year to the target
# cash flow balance, as a percentage of taxable payroll.
fit = fit_dep_path_vector(ss, dep_ratio_guess, target)
dep_fit = fit.dep_path_fit

using JLD2
@save joinpath(REPO, "cached_objects.jld2") ss dep_fit
```

The function `fit_dep_path_vector` does one Brent search for each year of the
projection. Each step of each search calculates a full projection. The function
is therefore slow. But it uses `project_economy` only. It does not solve a
steady state.

Keep the old cache file until you confirm that the new baseline projection still
agrees with the Trustees Report.

  A one line change to the defaults of `solve_steady_state` corrects this, at
  the cost of new numbers for the two `primus_reform.jl` scenarios.
