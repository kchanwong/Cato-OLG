# Import Function 
include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using JLD2;
proj = project_economy(
    ss,
    n_years         = 75,
    start_year      = 2025,
    g_A             = 0.0114,
    g_pop           = 0.005,
    inflation       = 0.024,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    trust_fund_rate = 0.047,
    gdp_anchor      = 28e12
);
proj_ADD05 = project_economy(
    ss,
    n_years         = 75,
    start_year      = 2025,
    g_A             = 0.0114,
    g_pop           = 0.005,
    inflation       = 0.029,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    trust_fund_rate = 0.047,
    gdp_anchor      = 28e12
);
proj_ADD1 = project_economy(
    ss,
    n_years         = 75,
    start_year      = 2025,
    g_A             = 0.0114,
    g_pop           = 0.005,
    inflation       = 0.034,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    trust_fund_rate = 0.047,
    gdp_anchor      = 28e12);
plot(proj.year, 100 * proj.ss_cash_flow_nom./proj.taxable_payroll_nom, label = "Baseline", xlabel = "Year", 
ylabel = "% of Taxable Payroll", title = "Projected Outlays",
ylim = (-7, 0))

plot!(proj.year, 100 * proj_ADD05.ss_cash_flow_nom./proj_ADD05.taxable_payroll_nom, lwd = 3, label = "Add 0.5% Inflation")

XLSX.openxlsx(joinpath(REPO, "output", "projections_inflation.xlsx"), mode="w") do xf
    for (name, p) in [("Baseline", proj), ("Add0_5pct", proj_ADD05), ("Add1pct", proj_ADD1)]
        sheet = XLSX.addsheet!(xf, name)
        # Header row
        cols = [:year, :ss_cash_flow_nom, :taxable_payroll_nom, :ss_outlays_nom,
                :fica_revenue_nom, :ss_cost_rate, :ss_cash_flow_pct, :ss_balance_pct,
                :trust_fund_nom, :avg_ben_per_retiree_nom, :GDP_nominal]
        for (j, col) in enumerate(cols)
            sheet[1, j] = string(col)
        end
        # Data rows
        vals = [getfield(p, col) for col in cols]
        for i in 1:length(p.year)
            for (j, v) in enumerate(vals)
                sheet[i+1, j] = v[i]
            end
        end
    end
end