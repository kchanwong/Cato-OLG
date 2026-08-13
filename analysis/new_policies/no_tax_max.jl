# Import Function 
include("C:/Users/kchanwong/Documents/PWBM/julia_port/functions_pwbm_w_spouse.jl")
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using JLD2;
@load "C:/Users/kchanwong/Documents/PWBM/julia_port/cached_objects.jld2" ss dep_fit
# Increase Payroll Tax Rate #
par_reform = create_params()
par_reform[:payroll_tax_rate] = 100000000
ss_reform  = solve_steady_state(par_reform);
proj_baseline = project_economy(ss,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0114,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.3e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);
proj_no_tax_max = project_economy(ss,
    ss_reform        = ss_reform,
    reform_year      = 2034,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0114,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.3e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);
proj_static = static_score_economy(ss,
    ss_reform        = ss_reform,
    reform_year      = 2034,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0114,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.3e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);
comparison = compare_projections(ss, ss_reform, 2034);
function proj_to_df(p)
    DataFrame(
        year                    = collect(p.year),
        ss_cash_flow_nom        = p.ss_cash_flow_nom,
        taxable_payroll_nom     = p.taxable_payroll_nom,
        ss_outlays_nom          = p.ss_outlays_nom,
        fica_revenue_nom        = p.fica_revenue_nom,
        ss_cost_rate            = p.ss_cost_rate,
        ss_cash_flow_pct        = p.ss_cash_flow_pct,
        ss_balance_pct          = p.ss_balance_pct,
        trust_fund_nom          = p.trust_fund_nom,
        avg_ben_per_retiree_nom = p.avg_ben_per_retiree_nom,
        GDP_nominal             = p.GDP_nominal
    )
end;
compare_projections(proj_static, proj_no_tax_max, proj_baseline)
const outpath_rr = "C:/Users/kchanwong/Documents/PWBM/julia_port/ROMINA_IVANE_PWBM_PAPER/new_policies/revenue_raisers.xlsx"
XLSX.openxlsx(outpath_rr, mode="rw") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "proj_no_tax_max")
    df = proj_to_df(proj_no_tax_max)
    for (c, col) in enumerate(names(df))
        sheet[1, c] = col
    end
    for (r, row) in enumerate(eachrow(df))
        for (c, val) in enumerate(row)
            sheet[r + 1, c] = val
        end
    end
end
println("revenue_raisers.xlsx written → $outpath_rr")