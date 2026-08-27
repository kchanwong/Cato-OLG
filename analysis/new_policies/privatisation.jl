# privatisation.jl
# Simulate partial privatisation: divert X% of FICA revenues to private accounts
# for the first `diversion_years` years. Outlays unchanged; trust fund re-simulated.

include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))
using DataFrames, XLSX, CSV, Printf, Statistics, JLD2

# ── baseline projection (current law, no reform) ─────────────────────────────
proj_base = project_economy(ss,
    n_years         = 75,
    start_year      = 2025,
    g_A             = 0.0113,
    g_pop           = 0.005,
    inflation       = 0.024,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    gdp_anchor      = 28e12,
    trust_fund_rate = 0.047,
    label           = "Current Law (privatisation baseline)");
plot(proj_base.taxable_payroll_nom)
# ── simulation function ───────────────────────────────────────────────────────
function simulate_privatisation(proj;
    diversion_frac::Float64  = 0.20,
    diversion_years::Int     = 30,
    trust_fund_init::Float64 = 2.76e12,
    trust_fund_rate::Float64 = 0.047,
    label::String            = "Privatisation $(round(Int, diversion_frac*100))%"
)
    n = length(proj.year)

    fica_revenue_nom        = zeros(n)
    diverted_nom            = zeros(n)
    trust_fund_interest_nom = zeros(n)
    ss_cash_flow_nom        = zeros(n)
    ss_total_income_nom     = zeros(n)
    ss_balance_nom          = zeros(n)
    trust_fund_nom          = zeros(n)

    trust_fund = trust_fund_init / 1e9   # work in $B

    for t in 1:n
        divert     = t <= diversion_years ? diversion_frac : 0.0
        fica_t     = proj.fica_revenue_nom[t] * (1.0 - divert)   # $B after diversion
        diverted_t = proj.fica_revenue_nom[t] * divert            # $B redirected
        btax_t     = proj.ss_benefit_tax_oasi_nom[t]              # $B (unchanged)
        outlays_t  = proj.ss_outlays_nom[t]                       # $B (unchanged)

        tf_interest_t  = trust_fund > 0 ? trust_fund * trust_fund_rate : 0.0
        ss_cash_flow_t = fica_t + btax_t - outlays_t
        ss_total_inc_t = fica_t + btax_t + tf_interest_t
        ss_balance_t   = ss_total_inc_t - outlays_t
        trust_fund    += ss_balance_t

        fica_revenue_nom[t]        = fica_t
        diverted_nom[t]            = diverted_t
        trust_fund_interest_nom[t] = tf_interest_t
        ss_cash_flow_nom[t]        = ss_cash_flow_t
        ss_total_income_nom[t]     = ss_total_inc_t
        ss_balance_nom[t]          = ss_balance_t
        trust_fund_nom[t]          = trust_fund
    end

    taxable_pay      = max.(proj.taxable_payroll_nom, 1e-6)
    ss_cash_flow_pct = ss_cash_flow_nom ./ taxable_pay .* 100
    ss_balance_pct   = ss_balance_nom   ./ taxable_pay .* 100
    ss_cost_rate     = proj.ss_outlays_nom ./ taxable_pay .* 100

    dep_yr = findall(v -> v < 0, trust_fund_nom)
    depl   = isempty(dep_yr) ? nothing : proj.year[dep_yr[1]]
    cf_yr  = findall(v -> v < 0, ss_cash_flow_nom)
    cf_def = isempty(cf_yr)  ? nothing : proj.year[cf_yr[1]]

    @printf("\n=== %s ===\n", label)
    @printf("  Diversion: %d%% of FICA for years %d-%d\n",
            round(Int, diversion_frac*100), proj.year[1], proj.year[min(diversion_years,n)])
    @printf("  TF depletion:      %s\n",
            isnothing(depl)   ? "Solvent through $(proj.year[end])" : string(depl))
    @printf("  Cash-flow deficit: %s\n",
            isnothing(cf_def) ? "None" : string(cf_def))
    @printf("  %-6s %10s %10s %10s %10s\n", "Window","FICA(\$B)","Outlays(\$B)","Balance(\$B)","TF(\$T)")
    for w in filter(v -> v <= n, [10, 30, 75])
        @printf("  %-6d %10.1f %10.1f %10.1f %10.2f\n",
                w,
                sum(fica_revenue_nom[1:w]),
                sum(proj.ss_outlays_nom[1:w]),
                sum(ss_balance_nom[1:w]),
                trust_fund_nom[w] / 1000)
    end

    (year                    = proj.year,
     fica_revenue_nom        = fica_revenue_nom,
     diverted_nom            = diverted_nom,
     ss_outlays_nom          = proj.ss_outlays_nom,
     taxable_payroll_nom     = proj.taxable_payroll_nom,
     ss_benefit_tax_oasi_nom = proj.ss_benefit_tax_oasi_nom,
     trust_fund_interest_nom = trust_fund_interest_nom,
     ss_cash_flow_nom        = ss_cash_flow_nom,
     ss_total_income_nom     = ss_total_income_nom,
     ss_balance_nom          = ss_balance_nom,
     trust_fund_nom          = trust_fund_nom,
     avg_ben_per_retiree_nom = proj.avg_ben_per_retiree_nom,
     GDP_nominal             = proj.GDP_nominal,
     ss_cost_rate            = ss_cost_rate,
     ss_cash_flow_pct        = ss_cash_flow_pct,
     ss_balance_pct          = ss_balance_pct,
     depletion_year          = depl,
     cashflow_deficit_year   = cf_def,
     label                   = label)
end

# ── run three scenarios ───────────────────────────────────────────────────────
proj_priv_20 = simulate_privatisation(proj_base,
    diversion_frac=0.20, diversion_years=30,
    trust_fund_init=2.8e12, trust_fund_rate=0.047, label="Privatisation 20%");

proj_priv_30 = simulate_privatisation(proj_base,
    diversion_frac=0.30, diversion_years=30,
    trust_fund_init=2.8e12, trust_fund_rate=0.047, label="Privatisation 30%");

proj_priv_50 = simulate_privatisation(proj_base,
    diversion_frac=0.50, diversion_years=30,
    trust_fund_init=2.8e12, trust_fund_rate=0.047, label="Privatisation 50%");

# ── export to XLSX ────────────────────────────────────────────────────────────
function priv_to_df(p; n=30)
    DataFrame(
        year                    = collect(p.year[1:n]),
        ss_cash_flow_nom        = p.ss_cash_flow_nom[1:n],
        taxable_payroll_nom     = p.taxable_payroll_nom[1:n],
        ss_outlays_nom          = p.ss_outlays_nom[1:n],
        fica_revenue_nom        = p.fica_revenue_nom[1:n],
        ss_cost_rate            = p.ss_cost_rate[1:n],
        ss_cash_flow_pct        = p.ss_cash_flow_pct[1:n],
        ss_balance_pct          = p.ss_balance_pct[1:n],
        trust_fund_nom          = p.trust_fund_nom[1:n],
        avg_ben_per_retiree_nom = p.avg_ben_per_retiree_nom[1:n],
        GDP_nominal             = p.GDP_nominal[1:n]
    )
end

scenarios = [
    ("privatisation_20pct", proj_priv_20),
    ("privatisation_30pct", proj_priv_30),
    ("privatisation_50pct", proj_priv_50),
]

outpath = joinpath(
    joinpath(REPO, "output", "new_policies"),
    "privatization_scenarios.xlsx")

XLSX.openxlsx(outpath, mode="w") do xf
    for (i, (name, proj)) in enumerate(scenarios)
        df    = priv_to_df(proj)
        sheet = i == 1 ? xf[1] : XLSX.addsheet!(xf, name)
        i == 1 && XLSX.rename!(sheet, name)
        for (c, col) in enumerate(names(df))
            sheet[1, c] = col
        end
        for (r, row) in enumerate(eachrow(df))
            for (c, val) in enumerate(row)
                sheet[r + 1, c] = val
            end
        end
    end
end

println("\nprivatization_scenarios.xlsx written → $outpath")
println("Sheets: ", join(first.(scenarios), ", "));
