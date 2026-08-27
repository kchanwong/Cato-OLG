# Import Function 
include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using JLD2;

function benefits_no_credit(AIME1::Float64, AIME2::Float64,
                             par::Dict, status::Symbol)
    # AIME is in model units (normalized by mu_dollar)
    # original cap in model units = ss_cap_dollars / mu_dollar
    # ss_cap_dollars ≈ 176100 / year = 14675 / month
    # in model units = (176100/12) / mu_dollar... 
    # but actually ss_bend2 and ss_cap are already stored in dollars
    # so we need: cap_in_model_units = par[:original_cap_dollars] / par[:mu_dollar]
    cap_mu = get(par, :original_cap_dollars, 176100.0) / par[:mu_dollar]
    AIME1_capped = min(AIME1, cap_mu)
    AIME2_capped = min(AIME2, cap_mu)
    benefits_current_law(AIME1_capped, AIME2_capped, par, status)
end
# Solve baseline first — mu_dollar calibrates freely
par_nc = create_params(benefit_fn = benefits_no_credit)
par_nc[:payroll_cap]         = Inf
ss_nc = solve_steady_state(par_nc, benefit_fn = benefits_no_credit);
par_wc = create_params()
par_wc[:payroll_cap]          = Inf
par_wc[:ss_aime_cap_dollars]  = Inf
par_wc[:ss_cap_frac] = Inf
ss_wc = solve_steady_state(par_wc);
# No Credits #
comparison_wc = compare_projections(ss, ss_wc, 2034);
comparison_nc = compare_projections(ss, ss_nc, 2025);
function proj_to_df(p)
    DataFrame(
        year                    = collect(p.year),
        ss_cash_flow_nom        = p.ss_cash_flow_nom,
        taxable_payroll_nom     = p.taxable_payroll_nom,
        ss_outlays_nom          = p.ss_outlays_nom,
        revenue_nom             = p.fica_revenue_nom + p.ss_benefit_tax_total_nom,
        ss_cost_rate            = p.ss_cost_rate,
        ss_cash_flow_pct        = p.ss_cash_flow_pct,
        ss_balance_pct          = p.ss_balance_pct,
        trust_fund_nom          = p.trust_fund_nom,
        avg_ben_per_retiree_nom = p.avg_ben_per_retiree_nom,
        GDP_nominal             = p.GDP_nominal
    )
end;
proj_static_wc = static_score_economy(ss,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0113,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.3e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);
proj_dynamic_nc = project_economy(ss,
    ss_reform        = ss_nc,
    reform_year      = 2025,
    n_years          = 81,
    start_year       = 2033,
    g_A              = 0.0113,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.3e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);
proj_static_nc = static_score_economy(ss,
    ss_reform        = ss_nc,
    reform_year      = 2034,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0113,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.3e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);
proj_dynamic_wc = project_economy(ss,
    ss_reform        = ss_wc,
    reform_year      = 2026,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0113,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.3e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);
proj_static_wc = static_score_economy(ss,
    ss_reform        = ss_wc,
    reform_year      = 2034,
    n_years          = 75,
    start_year       = 2025,
    g_A              = 0.0113,
    g_pop            = 0.005,
    inflation        = 0.024,
    dep_path         = dep_fit,
    ss_cola          = "wage",
    trust_fund_init  = 2.3e12,
    trust_fund_rate  = 0.047,
    gdp_anchor       = 28e12);
const outpath_sd = joinpath(REPO, "output", "new_policies", "static_dynamic_projections.xlsx")

# collect named DataFrames to write
sheets = Dict(
    "proj_dynamic_nc" => proj_to_df(proj_dynamic_nc),
    "proj_static_nc"  => proj_to_df(proj_static_nc),
    "proj_dynamic_wc" => proj_to_df(proj_dynamic_wc),
    "proj_static_wc"  => proj_to_df(proj_static_wc)
)

XLSX.openxlsx(outpath_sd, mode="w") do xf
    for (name, df) in sheets
        XLSX.addsheet!(xf, name)
        sheet = xf[name]
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

println("static_dynamic_projections.xlsx written → $outpath_sd")
