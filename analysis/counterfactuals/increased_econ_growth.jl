# Import Function 
include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using Optim;
# Save ss and dep_fit to cache for use by all other scripts

proj = project_economy(
    ss,
    n_years         = 75,
    start_year      = 2025,
    g_A             = 0.05,
    g_pop           = 0.005,
    inflation       = 0.024,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    trust_fund_rate = 0.047,
    gdp_anchor      = 28e12
);
vals = [
    13.34,  # 2025
    13.55,  # 2026
    13.71,  # 2027
    13.93,  # 2028
    14.03,  # 2029
    14.13,  # 2030
    14.19,  # 2031
    14.30,  # 2032
    14.39,  # 2033
    14.45,  # 2034
    14.51,  # 2035
    14.51,  # 2036
    14.64,  # 2037
    14.64,  # 2038
    14.74,  # 2039
    14.77,  # 2040
    14.81,  # 2041
    14.81,  # 2042
    14.84,  # 2043
    14.86,  # 2044
    14.86,  # 2045
    14.86,  # 2046
    14.86,  # 2047
    14.86,  # 2048
    14.94,  # 2049
    14.98,  # 2050
    15.02,  # 2051
    15.06,  # 2052
    15.11,  # 2053
    15.15,  # 2054
    15.21,  # 2055
    15.27,  # 2056
    15.35,  # 2057
    15.44,  # 2058
    15.52,  # 2059
    15.59,  # 2060
    15.73,  # 2061
    15.88,  # 2062
    16.01,  # 2063
    16.08,  # 2064
    16.14,  # 2065
    16.21,  # 2066
    16.34,  # 2067
    16.34,  # 2068
    16.48,  # 2069
    16.48,  # 2070
    16.48,  # 2071
    16.57,  # 2072
    16.57,  # 2073
    16.68,  # 2074
    16.68,  # 2075
    16.74,  # 2076
    16.82,  # 2077
    16.85,  # 2078
    16.87,  # 2079
    16.87,  # 2080
    16.87,  # 2081
    16.83,  # 2082
    16.76,  # 2083
    16.76,  # 2084
    16.64,  # 2085
    16.57,  # 2086
    16.57,  # 2087
    16.57,  # 2088
    16.46,  # 2089
    16.38,  # 2090
    16.38,  # 2091
    16.38,  # 2092
    16.28,  # 2093
    16.28,  # 2094
    16.27,  # 2095
    16.27,  # 2096
    16.27,  # 2097
    16.27,  # 2098
    16.27,  # 2099
]



100 * proj.ss_benefit_tax_total_nom./proj.taxable_payroll_nom
;
proj_ADD0_5 = project_economy(
    ss,
    n_years         = 75,
    start_year      = 2025,
    g_A             = 0.0163,
    g_pop           = 0.005,
    inflation       = 0.024,
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
    g_A             = 0.0213,
    g_pop           = 0.005,
    inflation       = 0.024,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    trust_fund_rate = 0.047,
    gdp_anchor      = 28e12
);
### Export ###
XLSX.openxlsx(joinpath(REPO, "output", "projections_gdp_increase.xlsx"), mode="w") do xf
    for (name, p) in [("Baseline", proj), ("Add0_5pct", proj_ADD0_5), ("Add1pct", proj_ADD1)]
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