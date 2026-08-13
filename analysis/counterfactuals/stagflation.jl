# Import Function 
include("C:/Users/kchanwong/Documents/PWBM/julia_port/functions_pwbm_w_spouse.jl")
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using JLD2;
@load "C:/Users/kchanwong/Documents/PWBM/julia_port/cached_objects.jld2" ss dep_fit
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
proj_STAG1 = project_economy(
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
    gdp_anchor      = 28e12,
    scenario_periods =
    [Dict(:year_start => 2033, :year_end => 2042,
             :g_A => -0.01/2, :inflation => 0.024)]
);
# Historical stagflation episode (1973-1982) mapped to 2032-2041
# Only g_A (AWI % chg) and inflation (CPI-W % chg) are changed
# g_A = real wage growth = (1 + AWI%) / (1 + CPI%) - 1
# inflation = CPI-W % change (nominal)
stag_hist = [
    Dict(:year_start => 2032, :year_end => 2032, :g_A =>  0.0006, :inflation => 0.0620),  # 1973
    Dict(:year_start => 2033, :year_end => 2033, :g_A => -0.0457, :inflation => 0.1101),  # 1974
    Dict(:year_start => 2034, :year_end => 2034, :g_A => -0.0155, :inflation => 0.0916),  # 1975
    Dict(:year_start => 2035, :year_end => 2035, :g_A =>  0.0107, :inflation => 0.0577),  # 1976
    Dict(:year_start => 2036, :year_end => 2036, :g_A => -0.0046, :inflation => 0.0648),  # 1977
    Dict(:year_start => 2037, :year_end => 2037, :g_A =>  0.0038, :inflation => 0.0753),  # 1978
    Dict(:year_start => 2038, :year_end => 2038, :g_A => -0.0241, :inflation => 0.1144),  # 1979
    Dict(:year_start => 2039, :year_end => 2039, :g_A => -0.0394, :inflation => 0.1348),  # 1980
    Dict(:year_start => 2040, :year_end => 2040, :g_A => -0.0016, :inflation => 0.1025),  # 1981
    Dict(:year_start => 2041, :year_end => 2041, :g_A => -0.0046, :inflation => 0.0600),  # 1982
]
proj_STAG2 = project_economy(
    ss,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0114,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.8e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12,
    scenario_periods = stag_hist
);
plot(proj.year, 100 * proj.ss_cash_flow_nom./proj.taxable_payroll_nom, label = "Baseline", xlabel = "Year", 
ylabel = "% of Taxable Payroll", title = "Projected Outlays",
ylim = (-7, 0))
plot!(proj.year, 100 * proj_STAG1.ss_cash_flow_nom./proj_STAG1.taxable_payroll_nom, lwd = 3, label = "Add 0.5% Inflation");
plot!(proj.year, 100 * proj_STAG2.ss_cash_flow_nom./proj_STAG2.taxable_payroll_nom, lwd = 3, label = "Add 1% Inflation");
XLSX.openxlsx("C:/Users/kchanwong/Documents/PWBM/julia_port/projections_stagflation.xlsx", mode="w") do xf
for (name, p) in [("Baseline", proj), ("STAG1", proj_STAG1), ("STAG2", proj_STAG2)]
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