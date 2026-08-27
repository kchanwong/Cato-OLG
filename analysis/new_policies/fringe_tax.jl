# Import Function 
include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using JLD2;
# Reform
par_reform = create_params()
par_reform[:aime_nwc_frac] = 0.2
par_reform[:payroll_rate]  *= (1 + 0.2)
ss_reform  = solve_steady_state(par_reform);
proj_raise_fringe = project_economy(ss,
    ss_reform        = ss_reform,
    reform_year      = 2029,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0114,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.8e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);

# --- Write to revenue_raisers.xlsx ---
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
end

const outpath_rr = joinpath(REPO, "output", "new_policies", "revenue_raisers.xlsx")
XLSX.openxlsx(outpath_rr, mode="rw") do xf
    sheet = XLSX.addsheet!(xf, "proj_raise_fringe")
    df = proj_to_df(proj_raise_fringe)
    for (c, col) in enumerate(names(df))
        sheet[1, c] = col
    end
    for (r, row) in enumerate(eachrow(df))
        for (c, val) in enumerate(row)
            sheet[r + 1, c] = val
        end
    end
end
