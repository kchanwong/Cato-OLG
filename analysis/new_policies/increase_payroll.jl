# Import Function 
include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using JLD2;
# Increase Payroll Tax Rate #
par_reform = create_params()
par_reform[:payroll_rate] = 0.17
ss_reform  = solve_steady_state(par_reform);

proj_baseline = project_economy(ss,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0114,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.8e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12)
proj_raise_payroll = project_economy(ss,
    ss_reform        = ss_reform,
    reform_year      = 2026,
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
proj_raise_payroll.avg_earnings_nominal
proj_baseline.avg_earnings_nominal
ss_reform

print(proj_raise_payroll.ss_cash_flow_pct)

const outpath_rr = joinpath(REPO, "output", "new_policies", "revenue_raisers.xlsx")
XLSX.openxlsx(outpath_rr, mode="rw") do xf
    sheet = XLSX.addsheet!(xf, "proj_raise_payroll")
    df = proj_to_df(proj_raise_payroll)
    for (c, col) in enumerate(names(df))
        sheet[1, c] = col
    end
    for (r, row) in enumerate(eachrow(df))
        for (c, val) in enumerate(row)
            sheet[r + 1, c] = val
        end
    end
end
