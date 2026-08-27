# Parameter reference

The model has three sets of controls. This document describes them in the order
in which you meet them.

1. [Model parameters](#1-model-parameters). These are the keys of the `par`
   dictionary that `create_params()` returns. They change household behavior and
   the general equilibrium. If you change one, you must solve a new steady state.
2. [Projection keywords](#2-projection-keywords). These are the keyword
   arguments of `project_economy` and `static_score_economy`. They set the
   exogenous macroeconomic paths and the accounting assumptions. They are cheap
   to change.
3. [Projection outputs](#3-projection-outputs). These are the results.

## Units

The model calculates in model units. One model unit is `par[:mu_dollar]`
dollars.

The solver calibrates `mu_dollar` while it runs. It sets `mu_dollar` so that
average earnings agree with `target_avg_earn`. Therefore, do not store a
parameter as a value in model units. Store it in dollars and divide it by
`mu_dollar` at the point of use. The solver recalculates two parameters at each
equilibrium iteration for this reason: `spousal_cap` and `min_benefit_floor`.

---

## 1. Model parameters

The function `create_params()` sets these parameters. To change one, assign a new
value before you solve the steady state:

```julia
par = create_params()
par[:payroll_rate] = 0.124
ss_reform = solve_steady_state(par)
```

**Caution:** `par` is a `Dict{Symbol,Any}`. If you spell a key incorrectly, Julia
accepts the new key and no function reads it. The model does not give an error.
For example, there is no key `:payroll_tax_rate`. The correct key is
`:payroll_rate`. Check each key against this document. Then confirm that the
reform steady state is different from the baseline steady state.

### Demographics and life cycle

| Key | Default | Description |
|---|---|---|
| `:J_start` | 21 | The age at which a household enters the model. |
| `:J_retire` | 67 | The full retirement age. Labor supply is 0 from this age. Benefits start at this age. Change this key to model an increase in the retirement age. |
| `:J_max` | 100 | The maximum age. Survival is 0 after this age. |
| `:n_ages` | 80 | The number of age periods. This value must be equal to `J_max - J_start + 1`. The model does not calculate it. If you change `J_start` or `J_max`, change this key also. |
| `:n_years` | 35 | The number of years of highest indexed earnings in the average of AIME. This key is not the length of the projection. The length of the projection is the `n_years` keyword of `project_economy`. Decrease this key to make the averaging period shorter. Increase it, for example to 40, to make the period longer. |
| `:g_pop` | 0.000 | The growth of population in the cohort weights of the steady state. This key is not the same as the `g_pop` keyword of `project_economy`. |

Mortality is not a parameter. The model uses the SSA 2022 period life table.
The tables `SSA_Q_MALE` and `SSA_Q_FEMALE` are at the start of the source file.
Couples use the two tables for the change to a survivor state. Single households
use the average of the two tables. To change mortality, change these arrays.

### Preferences

| Key | Default | Description |
|---|---|---|
| `:beta` | 0.97 | The annual discount factor. |
| `:gamma` | 0.30 | The weight of leisure in the utility function. This key controls the elasticity of labor supply. |
| `:sigma` | 2.00 | The coefficient of relative risk aversion. It is also the inverse of the intertemporal elasticity of substitution. The projection uses it to calculate the rate of convergence. |

### Productivity and heterogeneity

| Key | Default | Description |
|---|---|---|
| `:z_perm_vals` | `[-0.45, 0.10, 0.85, 1.20]` | The permanent productivity types, in logs. |
| `:p_perm_vals` | `[0.40, 0.45, 0.12, 0.03]` | The population share of each type. Give values that sum to 1. The function `build_couple_types` divides the joint weights by their sum, so a different sum has no effect. The length must be equal to the length of `z_perm_vals`. |
| `:rho_assort` | 0.6 | The strength of assortative mating. A value of 0 gives random matching. A higher value gives more couples of the same type. The 4 types give the 16 couple types that the solver calculates in parallel. |
| `:rho_pers` | 0.990 | The AR(1) persistence of the productivity shock. |
| `:sigma_pers` | 0.007 | The variance of the innovation to the productivity shock. This value is a variance, not a standard deviation. The Tauchen method puts it on a grid of `n_z` points. |
| `:g_A_index` | 0.012 | The rate at which the model indexes past earnings for AIME. Current law indexes to the average wage at age 60. The code indexes to age 66 instead. Refer to the known problems in [EXTENDING.md](EXTENDING.md#known-problems). |
| `:aime_nwc_frac` | 0 | The fraction by which the model increases raw earnings before AIME. Use it to represent coverage of compensation that is not wages. |

Productivity by age is not a parameter. The function `age_efficiency` contains a
quadratic function of age.

### Production and capital

| Key | Default | Description |
|---|---|---|
| `:alpha` | 0.36 | The share of capital in the Cobb-Douglas production function. |
| `:delta` | 0.06 | The rate of depreciation. |
| `:zeta_income_corp` | 0.55 | The share of capital income in the corporate sector. |
| `:zeta_income_pass` | 0.45 | The share of capital income in pass-through entities. |
| `:nu1` | 3.5 | The curvature of the cost function for leverage. |
| `:nu2` | `nothing` | The scale of the cost of leverage. Do not set this key. The function `calibrate_nu2` calculates it at each iteration. |
| `:leverage_ratio_target` | 0.32 | The target ratio of debt to assets. The calibration of `nu2` uses this target. |

### Household taxes

| Key | Default | Description |
|---|---|---|
| `:ord_brackets`, `:ord_rates` | 2019 single schedule | The brackets and the marginal rates for ordinary income. The model applies them to single households and to survivors. The brackets are in nominal dollars. |
| `:pref_brackets`, `:pref_rates` | | The schedule for preferential income, such as capital gains and dividends, for single households. |
| `:ord_brackets_mfj`, `:ord_rates_mfj` | | The ordinary income schedule for married couples who file jointly. |
| `:pref_brackets_mfj`, `:pref_rates_mfj` | | The preferential schedule for married couples who file jointly. |
| `:payroll_rate` | 0.106 | The combined OASI payroll tax rate. Change this key to model a change to the payroll tax rate. |
| `:payroll_cap` | 176100.0 | The taxable maximum, in dollars. Set it to `Inf` to model the removal of the taxable maximum. Do not use `:ss_cap` or `:ss_cap_frac` for this purpose. |
| `:theta_lab_ORD` | 0.95 | The fraction of labor income that is taxable as ordinary income. |
| `:theta_corp_ORD`, `:theta_corp_PREF` | 0.60, 0.40 | The split of corporate distributions between ordinary treatment and preferential treatment. |
| `:theta_pass_ORD` | 1.00 | The fraction of pass-through income that the model taxes as ordinary income. |
| `:theta_ss_ORD` | 0.85 | A flat taxable share of benefits. The model uses this key only if `:ss_tax_use_provisional` is `false`. |
| `:ss_tax_use_provisional` | `true` | If `true`, the model uses the provisional income formula for the tax on benefits, with tiers of 50% and 85%. If `false`, it uses the flat share. |
| `:ss_tax_thresh1`, `:ss_tax_thresh2` | 32000, 44000 | The provisional income thresholds for the two tiers, for married couples who file jointly. Current law does not index these thresholds, and the model holds them at these nominal values. The projection keyword `thresh_sens` controls the effect of inflation on the taxable share. |
| `:tau_con` | 0.0 | The consumption tax rate. |
| `:tau_lumpsum` | 0.0 | A lump sum tax for each household, in dollars. |

### Corporate taxes

| Key | Default | Description |
|---|---|---|
| `:tau_statutory_corp` | 0.21 | The statutory corporate tax rate. |
| `:phi_exp_corp` | 0.50 | The share of investment that the firm expenses. |
| `:phi_int_corp` | 1.00 | The share of interest that the firm deducts. The calibration of `nu2` uses this key. |
| `:zeta_ded_corp` | 0.02 | Other deductions, as a fraction of the tax base. |

### Social Security benefits

| Key | Default | Description |
|---|---|---|
| `:first_rr`, `:second_rr`, `:third_rr` | 0.90, 0.32, 0.15 | The replacement factors of the PIA formula for the three AIME brackets. To model price indexing of benefits, decrease all three factors. **Caution:** the `pia_factor_indexing` branch of the projection contains the values 0.90, 0.32, and 0.15 in the code. Do not use both mechanisms together. |
| `:ss_bend1`, `:ss_bend2` | 1226 × 12, 7391 × 12 | The bend points of the PIA formula, in annual dollars. |
| `:ss_cap_frac` | 2.796 | The cap on earnings for AIME, as a multiple of the average wage. This cap applies to creditable earnings. It is not the taxable maximum for FICA, which is `:payroll_cap`. |
| `:spousal_cap` | calculated | The maximum spousal supplement under `benefits_capped_spousal`. The intent is the PIA of a full time worker at the minimum wage. The code divides the annual earnings by 12, so the value is much lower than that PIA. Refer to the known problems in [EXTENDING.md](EXTENDING.md#known-problems). The solver recalculates it at each iteration. |
| `:min_benefit_floor_pct` | 1.25 | The minimum benefit floor, as a multiple of the federal poverty line. The function `benefits_price_index_floor` uses this key. |
| `:fpl_single_2025` | 15650.0 | The federal poverty line for one person, from the 2025 HHS guideline. |
| `:min_benefit_floor` | calculated | The floor in model units. The solver recalculates it at each iteration. |
| `:benefit_fn` | `benefits_current_law` | The benefit rule. `create_params` stores it here, but no function reads it. Give the rule to `solve_steady_state` instead, as `solve_steady_state(par; benefit_fn = ...)`. Refer to [EXTENDING.md](EXTENDING.md). |
| `:aime_modifier` | `nothing` | An optional function that changes the earnings history before the model calculates AIME. The model supplies `caregiver_modifier`. `create_params` stores it here, but no function reads it. Give it to `solve_steady_state` instead. |

The function `caregiver_modifier` reads five more keys: `:cg_start_age` (25),
`:cg_end_age` (38), `:cg_credit_level` (0.50 of the average wage),
`:cg_n_threshold` (0.10) and `:cg_max_years` (10). If hours of work are less
than `:cg_n_threshold`, the model counts the year as care of a family member.

### Government

| Key | Default | Description |
|---|---|---|
| `:r_G` | 0.03 | The real interest rate on government debt. |
| `:debt_to_gdp` | 0.78 | The ratio of debt to GDP in the steady state. It sets the supply of government debt in household portfolios. |
| `:G_residual_share` | 0.18 | Government spending other than Social Security, as a share of GDP. |

### Numerical settings and calibration

| Key | Default | Description |
|---|---|---|
| `:n_a` | 30 | The number of points on the asset grid. |
| `:n_z` | 5 | The number of points on the productivity grid. The time to solve increases with the square of this value for each couple type. |
| `:a_min`, `:a_max` | 0.0, 50.0 | The limits of the asset grid, in model units. Make sure that households at the top of the wealth distribution do not collect at `a_max`. |
| `:tol_equil` | 1e-4 | The tolerance for convergence of the ratio K/L. |
| `:max_iter_equil` | 15 | The maximum number of equilibrium iterations. Increase this value if the solver prints `Did not converge`. |
| `:mu_dollar` | 50000.0 | The number of dollars in one model unit. This value is only a first estimate. The solver recalculates it at each iteration. |
| `:avg_wage_mu` | 1.0 | Average earnings in model units. The AIME cap and the level of the credit for care of a family member use it. This value is only a first estimate. Do not set it. The solver replaces it at each iteration. |
| `:target_avg_earn` | 63000.0 | The target for average annual earnings, in dollars. This target sets `mu_dollar`. |

The solver writes eight keys into `par` while it runs. Do not set them
yourself. They are `:nu2`, `:mu_dollar`, `:avg_wage_mu`, `:spousal_cap`,
`:min_benefit_floor`, `:pi_approx`, `:couple_types`, and `:couple_weights`. The
last three are not present until you solve a steady state.

The solver damps the update of K/L. The new value is `0.7 × old + 0.3 × new`.
Therefore 15 iterations are not many. A reform that changes capital by a large
amount can need more iterations.

### Parameters that do nothing

The function `create_params()` sets these keys, but no function reads them. If
you change one, nothing happens. This list is here to prevent lost time.

`:benefit_fn` and `:aime_modifier` also belong in this list. They are the most
harmful two, because `create_params` accepts them as keyword arguments and a
caller expects them to work. Give both to `solve_steady_state` instead. The
other unread keys are:

`:eta`, `:h`, `:tau_top_pit`, `:phi_exp_pass`, `:phi_int_pass`,
`:zeta_other_pass`, `:zeta_taxbase_corp`, `:zeta_cred_corp`,
`:zeta_other_corp`, `:sigma_trans`, `:ss_thresh_indexing`, `:ss_bend1_frac`,
`:ss_bend2_frac`, `:ss_cap`,
`:closure_year`, `:r_K_world`, `:zeta_debt_foreign_takeup`,
`:zeta_capital_foreign_takeup`, `:tau_corp_FOR`, `:tau_pass_FOR`,
`:T_residual_share`, `:theta_lab_PT`, `:tol_vfi`, `:max_iter_vfi`.

Some of these keys are for functions that the model does not have. Examples are
an open economy and a value function iteration solver. The model solves the
household problem with the endogenous grid method, which does not iterate to a
tolerance. To connect one of these keys, you must change the model code.

---

## 2. Projection keywords

The two functions do not accept the same keywords. Each table below gives the
limits. If you give a keyword to a function that does not accept it, Julia stops
with an error for an unsupported keyword argument.

### The reform

| Keyword | Default | Description |
|---|---|---|
| `ss_reform` | `nothing` | The steady state of the reform. Give this keyword and `reform_year` together. The projection then moves from baseline behavior to reform behavior. For a baseline projection, give neither keyword. |
| `reform_year` | `nothing` | The calendar year in which the reform starts. |
| `label` | `"Current Law"` | The name of the run in the printed output. |
| `baseline_proj` | `nothing` | The result of a baseline projection. If you give this keyword, the function prints the score table. The table shows the change in the actuarial balance and the share of the deficit that the reform removes. Without this keyword, the function prints levels only. |

In `project_economy`, the function `compute_convergence_rate` controls the change
from the baseline aggregates to the reform aggregates. The rate is the stable
eigenvalue of the linear approximation to the Ramsey system. For the cached
baseline the half life is about 11 years. The function `static_score_economy` does not do this
calculation. In a static score, the policy changes at `reform_year`, but the
aggregates do not move.

### Horizon and macroeconomic paths

| Keyword | Default | Description |
|---|---|---|
| `n_years` | 75 | The length of the projection, in years. SSA uses 75 years. |
| `start_year` | 2025 | The first calendar year. |
| `g_A` | 0.012 | The annual growth of real productivity. |
| `g_pop` | 0.005 | The annual growth of population. |
| `inflation` | 0.025 | The annual rate of headline CPI. |
| `cola_inflation` | `nothing` | The measure of inflation for the COLA, if it is different from `inflation`. To score a reform that decreases the COLA, set this keyword below `inflation`. Chained CPI is an example. |
| `dep_path` | `nothing` | A vector of the ratio of retirees to workers, with one value per projection year. Give `dep_fit`. If the vector is shorter than `n_years`, the function repeats the last value. If it is longer, the function removes the extra values. The default builds a general logistic path that does not agree with the Trustees Report. To model a reform that changes the number of retirees, scale this path. The file `analysis/new_policies/raise_retirement_age.jl` gives an example. |
| `scenario_periods` | `nothing` | A vector of dictionaries that replace values for a range of years. An example is `[Dict(:year_start=>2030, :year_end=>2040, :g_A=>0.005, :inflation=>0.05)]`. The keys `:g_A`, `:inflation`, `:g_pop`, and `:cola_inflation` are available. The counterfactual scripts use this keyword. |

### Social Security accounting

| Keyword | Default | Description |
|---|---|---|
| `ss_cola` | `"wage"` | The growth of the first benefit of a new claimant. The value `"wage"` gives growth with productivity, which is current law. The value `"cpi"` gives no growth in real terms. A `Float64` gives a fixed real rate. This keyword is not the COLA on benefits that are already in payment. The keyword `cola_inflation` controls that COLA. |
| `cap_base_dollars` | 176100.0 | The taxable maximum in `start_year`. After that year it grows with average wages. |
| `trust_fund_init` | 2.3e12 | The balance of the trust fund in `start_year`. The analysis scripts give 2.8e12. Use the same value for the baseline and for the reform. |
| `trust_fund_rate` | 0.035 | The nominal interest rate on the balance of the trust fund. The function also uses it as the discount rate for the present value of the actuarial balance. The analysis scripts give 0.047. |
| `thresh_sens` | 0.3 | The sensitivity of the taxable share of benefits to inflation. The thresholds for the tax on benefits are fixed in nominal dollars, so the taxable share increases as the price level increases. |
| `pia_factor_indexing` | `false` | If `true`, the projection multiplies the PIA factors by the ratio of prices to wages after `reform_year`. This models price indexing of benefits inside the projection. The code contains the factors 0.90, 0.32, and 0.15. Therefore do not use this keyword together with a change to `:first_rr` and the other two factors. |
| `cola_cap` | `false` | For `project_economy` only. If `true`, the projection limits the size in dollars of the COLA increase each year. The limit is `chained_cpi × cola_cap_pct × FPL(t)`. This is the COLA cap of the Committee for a Responsible Federal Budget. The limit applies to the increase, not to the level of the benefit. |
| `chained_cpi` | 0.021 | For `project_economy` only. The rate of chained CPI for the COLA cap and for the path of the federal poverty line. |
| `fpl_single_2025` | 15650.0 | For `project_economy` only. The federal poverty line in `start_year`, for the COLA cap. |
| `cola_cap_pct` | 1.25 | For `project_economy` only. The multiple of the federal poverty line that sets the COLA cap. |
| `eliminates_cap` | `false` | For `static_score_economy` only. If `true`, the function sets the share of earnings above the taxable maximum to 0 after the reform. A static score has no change in behavior, so it cannot calculate the new taxable share. |

The function `static_score_economy` does not have the COLA cap. It does not
accept `cola_cap`, `chained_cpi`, `fpl_single_2025`, or `cola_cap_pct`.
Therefore you cannot make a static score of a COLA cap reform. Use
`project_economy` for that reform.

### Anchors for scale

| Keyword | Default | Description |
|---|---|---|
| `gdp_anchor` | `nothing` | The target GDP in `start_year`, in dollars. The scripts give `28e12`. This keyword sets the number of households, and therefore all levels in dollars. Without it, the function sets the scale from 180 million covered workers. |
| `wages_to_gdp` | 0.46 | Total covered payroll, as a share of GDP. |
| `bracket_indexing` | `"cpi"` | The function accepts this keyword but does not use it. |

---

## 3. Projection outputs

Both functions return a NamedTuple. Most series in nominal dollars are in
billions. There are three exceptions. The series `avg_ben_per_retiree_nom` is in
thousands. The series `avg_earnings_nominal` and `payroll_cap_nom` are in
dollars. The other series are years, ratios, index numbers, or percentages.

**Series.** Each series has `n_years` elements: `year`, `t`, `dep_ratio`,
`price_level`, `cum_productivity`, `cum_population`, `g_A_annual`,
`inflation_annual`, `cola_inflation_annual`, `GDP_real`, `GDP_nominal`,
`avg_earnings_nominal`, `payroll_cap_nom`, `total_pay_nom`, `pct_above_nom`,
`taxable_payroll_nom`, `fica_revenue_nom`, `ss_outlays_nom`,
`avg_ben_per_retiree_nom`, `ss_benefit_tax_oasi_nom`, `ss_benefit_tax_hi_nom`,
`ss_benefit_tax_total_nom`, `trust_fund_interest_nom`, `ss_cash_flow_nom`,
`ss_total_income_nom`, `ss_balance_nom`, `trust_fund_nom`, `fica_rate_eff`,
`ss_cost_rate`, `ss_cash_flow_pct`, `ss_balance_pct`, `total_tax_rev_nom`,
`govt_spending_nom`, `debt_nom`, `debt_to_gdp`.

The three series that end in `_pct` are percentages of taxable payroll. SSA uses
this denominator. The series `ss_cash_flow_nom` does not include interest on the
trust fund. The series `ss_balance_nom` includes it.

**Single values.**

| Field | Description |
|---|---|
| `actuarial_balance_pct` | The 75 year actuarial balance in present value, as a percentage of taxable payroll. This value is the primary result. |
| `annual_balance_75_pct` | The cash flow balance in the last year, as a percentage of taxable payroll. This value is the second result that SSA reports. |
| `deficit_eliminated` | For `project_economy` only. `static_score_economy` does not return this field. It is `true` only if `actuarial_balance_pct` is 0 or more. A reform can remove most of the deficit and this field can still be `false`. |
| `change_lab`, `change_ab75` | The change in the two primary results, against `baseline_proj`. Both are `nothing` if you do not give `baseline_proj`. |
| `score_lab`, `score_ab75` | The same two changes, as a percentage of the baseline deficit. These fields give the share of the deficit that the reform removes. Both are `nothing` if the baseline has no deficit. |
| `depletion_year` | The first year in which the trust fund is negative, or `nothing`. |
| `cashflow_deficit_year` | The first year in which cash flow is negative, or `nothing`. |
| `convergence_rate` | The eigenvalue of the transition. It is 0.0 in a baseline projection. |
| `label` | The value of the `label` keyword. |

## Steady state outputs

The function `solve_steady_state` returns a NamedTuple with these fields: `par`,
`grids`, `hh_list`, `x_list`, `agg`, `prices`, `K`, `L`, `Y`, `D`, `w`, `r_K`,
`KL`, `ps`, `dep_ratio`, `worker_mass`, `retiree_mass`, and `benefit_fn`.

The fields `hh_list` and `x_list` hold the policy functions and the distributions
for each couple type. Start an analysis of distribution here. The file
`analysis/new_policies/individual_provision_score.jl` gives an example. It takes
quantiles of AIME from `ss.hh_list`.

The function `extract_ss_ratios(ss)` returns the aggregates that a projection
uses. They include `Y`, `K`, `L`, `w`, `C`, `A`, `fica_rev`, `ss_ben`,
`ss_balance`, `dep_ratio`, `avg_earnings`, `taxable_payroll`, `pct_above_cap`,
`ss_ben_per_retiree`, `mu_dollar`, and related ratios for each person.

To confirm that a reform has the effect that you intend, compare
`extract_ss_ratios(ss)` with `extract_ss_ratios(ss_reform)`. This check is quick.
