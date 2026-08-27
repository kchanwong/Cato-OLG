# A template for a new analysis. Copy this file into analysis/new_policies/ or
# analysis/counterfactuals/. Then change the copy. Run it with this command:
#
#     julia --threads=auto analysis/new_policies/my_reform.jl
#
# Refer to docs/EXTENDING.md to find the level of your reform. A reform can be
# a projection keyword, a parameter change, a new benefit rule, or a change to
# the model code.

include(joinpath(@__DIR__, "..", "source", "setup.jl"))   # Add one ".." per subdirectory.
using DataFrames, XLSX

# The macroeconomic assumptions. The baseline and the reform must use the same
# assumptions. If they do not, the difference between them is not an effect of
# the policy.
const ASSUMPTIONS = (
    n_years         = 75,
    start_year      = 2025,
    g_A             = 0.0114,
    g_pop           = 0.005,
    inflation       = 0.024,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    trust_fund_rate = 0.047,
    gdp_anchor      = 28e12,
)

const REFORM_YEAR = 2034

# Step 1. Define the reform.
par_reform = create_params()
par_reform[:J_retire] = 70          # Put your reform here.

# Step 2. Solve the steady state of the reform. This takes about 1 minute.
ss_reform = solve_steady_state(par_reform)

# Confirm that the reform changed the steady state.
let b = extract_ss_ratios(ss), r = extract_ss_ratios(ss_reform)
    @printf("benefit/retiree: %.4f -> %.4f   taxable payroll: %.4f -> %.4f\n",
            b.ss_ben_per_retiree, r.ss_ben_per_retiree,
            b.taxable_payroll,    r.taxable_payroll)
end

# Step 3. Calculate the baseline, then the dynamic score, then the static score.
proj_base = project_economy(ss; ASSUMPTIONS..., label = "Current Law")

proj_reform = project_economy(ss;
    ss_reform     = ss_reform,
    reform_year   = REFORM_YEAR,
    baseline_proj = proj_base,          # This makes the score table print.
    label         = "My Reform",
    ASSUMPTIONS...)

proj_static = static_score_economy(ss;
    ss_reform     = ss_reform,
    reform_year   = REFORM_YEAR,
    baseline_proj = proj_base,
    label         = "My Reform (static)",
    ASSUMPTIONS...)

# Step 4. Write the results to a spreadsheet.
proj_to_df(p) = DataFrame(
    year                    = collect(p.year),
    taxable_payroll_nom     = p.taxable_payroll_nom,
    fica_revenue_nom        = p.fica_revenue_nom,
    ss_outlays_nom          = p.ss_outlays_nom,
    ss_cash_flow_nom        = p.ss_cash_flow_nom,
    ss_cost_rate            = p.ss_cost_rate,
    ss_cash_flow_pct        = p.ss_cash_flow_pct,
    ss_balance_pct          = p.ss_balance_pct,
    trust_fund_nom          = p.trust_fund_nom,
    avg_ben_per_retiree_nom = p.avg_ben_per_retiree_nom,
    GDP_nominal             = p.GDP_nominal,
)

outpath = joinpath(REPO, "output", "new_policies", "my_reform.xlsx")
XLSX.openxlsx(outpath, mode = "w") do xf
    for (name, p) in [("Baseline", proj_base),
                      ("Dynamic",  proj_reform),
                      ("Static",   proj_static)]
        sheet = XLSX.addsheet!(xf, name)
        df = proj_to_df(p)
        for (c, col) in enumerate(names(df)); sheet[1, c] = col; end
        for (r, row) in enumerate(eachrow(df)), (c, v) in enumerate(row)
            sheet[r + 1, c] = v
        end
    end
end
println("wrote $outpath")
