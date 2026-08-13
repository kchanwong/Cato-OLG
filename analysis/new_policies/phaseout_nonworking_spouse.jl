# Import Functions
include("C:/Users/kchanwong/Documents/PWBM/julia_port/functions_pwbm_w_spouse.jl")
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using Statistics
using JLD2;
@load "C:/Users/kchanwong/Documents/PWBM/julia_port/cached_objects.jld2" ss dep_fit

# ============================================================
# NEW BENEFIT FUNCTION: Phase out spousal benefit for top income quintile
# Households whose combined AIME >= 80th percentile get no spousal top-up.
# Survivor benefits are unchanged.
# ============================================================

function benefits_phaseout_top_quintile(AIME1::Float64, AIME2::Float64,
                                         par::Dict, status::Symbol)
    PIA1   = compute_PIA(AIME1, par)
    PIA2   = compute_PIA(AIME2, par)
    thresh = get(par, :top_quintile_aime_threshold, Inf)

    if status == :couple
        if (AIME1 + AIME2) >= thresh
            # Top quintile: no spousal top-up, each gets only own PIA
            return (PIA1, PIA2)
        else
            b1 = PIA1 + max(0.0, 0.5*PIA2 - PIA1)
            b2 = PIA2 + max(0.0, 0.5*PIA1 - PIA2)
            return (b1, b2)
        end
    elseif status == :survivor1   # spouse 1 alive, spouse 2 dead
        return (max(PIA1, PIA2), 0.0)
    elseif status == :survivor2   # spouse 2 alive, spouse 1 dead
        return (0.0, max(PIA2, PIA1))
    end
    (PIA1, PIA2)
end

# ============================================================
# SOLVE BASELINE
# ============================================================
# ============================================================
# COMPUTE 80th PERCENTILE OF JOINT AIME FROM BASELINE
# ============================================================
par = create_params()
nz = par[:n_z]
all_joint_aimes = [
    ss.hh_list[ict].AIME1[iz1, iz2] + ss.hh_list[ict].AIME2[iz1, iz2]
    for ict in 1:length(ss.hh_list)
    for iz1 in 1:nz, iz2 in 1:nz
]
quintile_threshold = quantile(all_joint_aimes, 0.80);
# ============================================================
# SOLVE REFORM STEADY STATE
# ============================================================
par_reform = create_params(benefit_fn = benefits_phaseout_top_quintile)
par_reform[:n_years] = 40
par_reform[:J_retire] = 70
par_reform[:max_iter_equil] = 300
par_reform[:top_quintile_aime_threshold] = quintile_threshold
ss_reform = solve_steady_state(par_reform; benefit_fn = benefits_phaseout_top_quintile);
# ============================================================
# PROJECT ECONOMY (75-year window)
# ============================================================
proj_base = project_economy(
    ss,
    n_years         = 75,
    start_year      = 2025,
    g_A             = 0.0113,
    g_pop           = 0.005,
    inflation       = 0.024,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    trust_fund_rate = 0.047,
    gdp_anchor      = 28e12
);
year_cal_full = collect(2025:2099)
fra_path = [yr < 2026 ? 67.0 :
            yr < 2037 ? 67.0 + 3.0 * (yr - 2026) / 8 :
            70.0 + (yr - 2037) / 24
            for yr in year_cal_full]
dep_path_le = [dep_fit[t] *
               ((100 - fra_path[t]) * (67 - 21)) /
               ((100 - 67)          * (fra_path[t] - 21))
               for t in 1:75];
function project_economy_w_sl(
    ss;
    n_years::Int          = 75,
    start_year::Int       = 2025,
    g_A::Float64          = 0.012,
    g_pop::Float64        = 0.005,
    inflation::Float64    = 0.025,
    cola_inflation::Union{Float64,Nothing} = nothing,  # chained CPI etc.; defaults to `inflation`
    dep_path              = nothing,
    ss_cola               = "wage",
    cap_base_dollars::Float64 = 176100.0,
    bracket_indexing::String  = "cpi",
    trust_fund_init::Float64  = 2.3e12,
    trust_fund_rate::Float64  = 0.035,
    ss_reform             = nothing,
    reform_year           = nothing,
    label::String         = "Current Law",
    scenario_periods      = nothing,
    thresh_sens::Float64  = 0.3,
    pia_factor_indexing::Bool = false,
    gdp_anchor            = nothing,
    wages_to_gdp::Float64 = 0.46,
    # --- Mandatory coverage of non-covered state/local workers ---
    mandate_coverage_year     = nothing,
    n_uncovered_base::Float64 = 6.5e6,
    new_covered_earn_ratio::Float64 = 1.0,
    pia_concavity::Float64    = 0.6,
    baseline_proj         = nothing
)
    @printf("=== Projecting %d years: %d-%d [%s] ===\n",
            n_years, start_year, start_year+n_years-1, label)

    base = extract_ss_ratios(ss)
    mu   = base.mu_dollar; par = ss.par

    ssa_covered_workers = 180e6
    N_hh = ssa_covered_workers / base.worker_mass
    if !isnothing(gdp_anchor)
        N_hh = Float64(gdp_anchor) / (base.Y * mu)
        @printf("  N_hh: %.1fM (anchored to GDP anchor \$%.2fT)\n",
                N_hh/1e6, gdp_anchor/1e12)
    else
        @printf("  N_hh: %.1fM (anchored to %sM covered workers)\n",
                N_hh/1e6, format_comma(round(ssa_covered_workers/1e6)))
    end
    @printf("  Implied GDP: \$%.1fT\n",
            base.Y * mu * N_hh / 1e12)

    has_reform = !isnothing(ss_reform) && !isnothing(reform_year)
    ref = nothing
    conv_rate = 0.0
    alpha = par[:alpha]

    K_base = ss.K;  L_base = ss.L
    K_ref  = 0.0;   L_ref  = 0.0

    if has_reform
        ref = extract_ss_ratios(ss_reform)
        K_ref = ss_reform.K;  L_ref = ss_reform.L
        @printf("  Reform at %d\n", reform_year)
        conv_rate = compute_convergence_rate(ss, g_A=g_A)
    end

    has_mandate = !isnothing(mandate_coverage_year)
    if has_mandate
        @printf("  Mandatory coverage at %d: %.1fM workers, earn_ratio=%.2f\n",
                mandate_coverage_year, n_uncovered_base/1e6, new_covered_earn_ratio)
    end

    if isnothing(dep_path)
        dep_path = build_dependency_path(n_years, dep_start=base.dep_ratio)
    end
    length(dep_path) < n_years &&
        (dep_path = vcat(dep_path, fill(dep_path[end], n_years - length(dep_path))))

    year_cal  = collect(start_year:(start_year+n_years-1))
    g_A_vec   = fill(g_A,       n_years)
    inf_vec   = fill(inflation,  n_years)
    g_pop_vec = fill(g_pop,     n_years)

    if !isnothing(scenario_periods)
        for ep in scenario_periods
            idx = findall(y -> y >= ep[:year_start] && y <= ep[:year_end], year_cal)
            isempty(idx) && continue
            haskey(ep,:g_A)       && (g_A_vec[idx]   .= ep[:g_A])
            haskey(ep,:inflation) && (inf_vec[idx]   .= ep[:inflation])
            haskey(ep,:g_pop)     && (g_pop_vec[idx] .= ep[:g_pop])
            @printf("  Scenario episode %d-%d: g_A=%.3f inf=%.3f\n",
                    ep[:year_start], ep[:year_end],
                    get(ep,:g_A,g_A), get(ep,:inflation,inflation))
        end
    end

    # COLA inflation measure (e.g., C-CPI-U ≈ CPI-U − 0.3pp)
    # When cola_inflation is nothing, COLA uses same measure as headline inflation.
    # When set (e.g. 0.021 for chained CPI vs 0.024 CPI-U), existing retirees'
    # real benefits erode by the wedge each year.
    cola_inf_vec = fill(isnothing(cola_inflation) ? inflation : Float64(cola_inflation), n_years)
    if !isnothing(scenario_periods)
        for ep in scenario_periods
            idx = findall(y -> y >= ep[:year_start] && y <= ep[:year_end], year_cal)
            haskey(ep, :cola_inflation) && (cola_inf_vec[idx] .= ep[:cola_inflation])
        end
    end

    cum_A     = cumprod([1.0; (1.0 .+ g_A_vec[2:end])])
    cum_pop   = cumprod([1.0; (1.0 .+ g_pop_vec[2:end])])
    cum_price = cumprod([1.0; (1.0 .+ inf_vec[2:end])])
    cum_wage  = cumprod([1.0; (1.0 .+ g_A_vec[2:end] .+ inf_vec[2:end])])
    cum_price_lag = [1.0; cum_price[1:end-1]]

    cap_nom_vec = awi_cap_path(cap_base_dollars, cum_wage)
    @printf("  Cap: \$%s -> \$%s (yr 10) -> \$%s (yr %d)\n",
            format_comma(cap_base_dollars),
            format_comma(cap_nom_vec[min(10,n_years)]),
            format_comma(cap_nom_vec[n_years]), n_years)

    cum_new_ben = ones(n_years)
    if ss_cola == "wage"
        cum_new_ben = copy(cum_A)
    elseif ss_cola == "cpi"
        cum_new_ben = ones(n_years)
    else
        rate = Float64(ss_cola)
        for t in 2:n_years; cum_new_ben[t] = cum_new_ben[t-1]*(1+rate); end
    end

    working_years = par[:J_retire] - par[:J_start]

    btax_rate_oasi_base = base.ss_benefit_tax_oasi_per_retiree /
                          max(base.ss_ben_per_retiree, 1e-12)
    btax_rate_hi_base   = base.ss_benefit_tax_hi_per_retiree /
                          max(base.ss_ben_per_retiree, 1e-12)
    btax_rate_oasi_ref = 0.0; btax_rate_hi_ref = 0.0
    if has_reform
        btax_rate_oasi_ref = ref.ss_benefit_tax_oasi_per_retiree /
                             max(ref.ss_ben_per_retiree, 1e-12)
        btax_rate_hi_ref   = ref.ss_benefit_tax_hi_per_retiree /
                             max(ref.ss_ben_per_retiree, 1e-12)
    end

    GDP_real                 = zeros(n_years)
    GDP_nominal              = zeros(n_years)
    avg_earnings_nominal     = zeros(n_years)
    payroll_cap_nom          = zeros(n_years)
    taxable_payroll_nom      = zeros(n_years)
    fica_revenue_nom         = zeros(n_years)
    ss_outlays_nom           = zeros(n_years)
    avg_ben_per_retiree_nom  = zeros(n_years)
    ss_benefit_tax_oasi_nom  = zeros(n_years)
    ss_benefit_tax_hi_nom    = zeros(n_years)
    ss_benefit_tax_total_nom = zeros(n_years)
    trust_fund_interest_nom  = zeros(n_years)
    ss_cash_flow_nom         = zeros(n_years)
    ss_total_income_nom      = zeros(n_years)
    ss_balance_nom           = zeros(n_years)
    trust_fund_nom_v         = zeros(n_years)
    fica_rate_eff            = zeros(n_years)
    ss_cost_rate             = zeros(n_years)
    ss_cash_flow_pct         = zeros(n_years)
    ss_balance_pct           = zeros(n_years)
    total_tax_rev_nom        = zeros(n_years)
    govt_spending_nom        = zeros(n_years)
    debt_nom                 = zeros(n_years)
    debt_to_gdp_v            = zeros(n_years)
    total_pay_nom            = zeros(n_years)
    pct_above_nom            = zeros(n_years)
    mandate_fica_nom         = zeros(n_years)
    mandate_outlays_nom      = zeros(n_years)

    # --- Precompute AIME and bend points for PIA factor indexing ---
    b1_mu = par[:ss_bend1] / mu
    b2_mu = par[:ss_bend2] / mu
    avg_aime_mu = b2_mu
    pia_base_level = base.ss_ben_per_retiree
    for _ in 1:20
        pia_try = 0.90*min(avg_aime_mu, b1_mu) +
                  0.32*max(0.0, min(avg_aime_mu, b2_mu) - b1_mu) +
                  0.15*max(0.0, avg_aime_mu - b2_mu)
        err = pia_try - pia_base_level
        mr  = avg_aime_mu <= b1_mu ? 0.90 : avg_aime_mu <= b2_mu ? 0.32 : 0.15
        avg_aime_mu -= err / max(mr, 0.01)
        avg_aime_mu  = max(avg_aime_mu, 0.01)
    end
    pia_current_law = 0.90*min(avg_aime_mu, b1_mu) +
                      0.32*max(0.0, min(avg_aime_mu, b2_mu) - b1_mu) +
                      0.15*max(0.0, avg_aime_mu - b2_mu)

    trust_fund   = trust_fund_init
    avg_ben_real = base.ss_ben_per_retiree

    for t in 1:n_years
        yr  = year_cal[t]

        # === Transition path (Mankiw-Weinzierl eqs 34-35) ===
        omega = 0.0
        if has_reform && yr >= reform_year
            tau   = Float64(yr - reform_year)
            decay = exp(conv_rate * tau)
            omega = 1.0 - decay

            K_mu = exp(log(K_ref) + (log(K_base) - log(K_ref)) * decay)
            L_mu = exp(log(L_ref) + (log(L_base) - log(L_ref)) * decay)
        else
            K_mu = K_base
            L_mu = L_base
        end

        Y_mu = K_mu^alpha * L_mu^(1-alpha)
        w_mu = (1-alpha) * (K_mu / L_mu)^alpha

        interp(bv, rv) = bv + omega * (rv - bv)
        fica_rate_t      = has_reform ? interp(base.payroll_rate,  ref.payroll_rate)  : base.payroll_rate
        pct_above_t      = has_reform ? interp(base.pct_above_cap, ref.pct_above_cap) : base.pct_above_cap
        inc_tax_mu       = has_reform ? interp(base.income_tax,    ref.income_tax)    : base.income_tax
        btax_oasi_rate_t = has_reform ? interp(btax_rate_oasi_base, btax_rate_oasi_ref) : btax_rate_oasi_base
        btax_hi_rate_t   = has_reform ? interp(btax_rate_hi_base,   btax_rate_hi_ref)   : btax_rate_hi_base
        new_ben_scale_t  = has_reform ?
            interp(base.ss_ben_per_retiree, ref.ss_ben_per_retiree) / base.ss_ben_per_retiree : 1.0

        zeta_corp = par[:zeta_income_corp]
        corp_profit_mu = zeta_corp * (Y_mu - w_mu * L_mu - par[:delta] * K_mu)
        corp_tax_mu = par[:tau_statutory_corp] *
                      max(0.0, corp_profit_mu) *
                      (1 - par[:phi_exp_corp]*0.3 - par[:zeta_ded_corp])

        dep_t     = dep_path[t]

        # --- PIA factor indexing ---
        if pia_factor_indexing && !isnothing(reform_year) && yr >= reform_year
            t_ref = findfirst(==(reform_year), year_cal)
            price_wage_ratio = isnothing(t_ref) ? 1.0 : cum_A[t_ref] / cum_A[t]
            f1 = 0.90 * price_wage_ratio
            f2 = 0.32 * price_wage_ratio
            f3 = 0.15 * price_wage_ratio
            pia_yr = f1*min(avg_aime_mu, b1_mu) +
                     f2*max(0.0, min(avg_aime_mu, b2_mu) - b1_mu) +
                     f3*max(0.0, avg_aime_mu - b2_mu)
            new_ben_scale_t = pia_yr / max(pia_current_law, 1e-10)
        end

        # === Scale model units → nominal dollars ===
        N_hh_t     = N_hh * cum_pop[t]
        GDP_real_t = Y_mu * cum_A[t] * mu * N_hh_t
        GDP_nom_t  = GDP_real_t * cum_price[t]
        avg_earn_t = (w_mu * L_mu / max(base.worker_mass, 1e-10)) *
                     cum_A[t] * mu * cum_price[t]

        # Payroll
        total_pay_nom_t   = wages_to_gdp * GDP_nom_t
        wage_cap_ratio    = cum_wage[t] / (cap_nom_vec[t] / cap_base_dollars)
        pct_above_adj     = max(0, pct_above_t * wage_cap_ratio)
        taxable_pay_nom_t = total_pay_nom_t * (1.0 - pct_above_adj)
        fica_rev_nom_t    = fica_rate_t * taxable_pay_nom_t

        # SS outlays
        phi_t          = max(0, 1.0 / (working_years * dep_t))
        new_ben_real_t = base.ss_ben_per_retiree * cum_new_ben[t] * new_ben_scale_t

        # Existing retirees' real benefits erode when COLA measure < headline CPI.
        # cola_erosion_t = 1.0 when cola_inflation == inflation (no wedge).
        # When cola_inflation < inflation (e.g. C-CPI-U), the ratio < 1 and
        # the stock of existing benefits loses real purchasing power each period.
        cola_erosion_t = (1.0 + cola_inf_vec[t]) / (1.0 + inf_vec[t])
        avg_ben_real   = (1.0 - phi_t) * avg_ben_real * cola_erosion_t + phi_t * new_ben_real_t

        avg_ben_nom    = avg_ben_real * mu * cum_price_lag[t]
        n_ret_t        = dep_t * ssa_covered_workers * cum_pop[t]
        ss_outlays_nom_t = avg_ben_nom * n_ret_t

        # === Mandatory coverage of non-covered state/local workers ===
        add_fica_t    = 0.0
        add_outlays_t = 0.0
        if has_mandate && yr >= mandate_coverage_year
            tau_m = yr - mandate_coverage_year

            # All newly covered workers pay FICA immediately
            n_new         = n_uncovered_base * cum_pop[t]
            earn_new      = avg_earn_t * new_covered_earn_ratio
            add_taxable   = n_new * earn_new * (1.0 - pct_above_adj)
            add_fica_t    = fica_rate_t * add_taxable

            # Benefits: AIME diluted by zero-earning pre-mandate years
            # career_frac = fraction of full career with covered earnings
            # PIA is concave in AIME → PIA(frac·AIME)/PIA(AIME) ≈ frac^pia_concavity
            years_covered = min(tau_m, working_years)
            career_frac   = years_covered / working_years
            pia_frac      = career_frac ^ pia_concavity

            n_new_ret     = n_new * dep_t
            add_outlays_t = avg_ben_nom * pia_frac * new_covered_earn_ratio * n_new_ret

            # Inject into flows
            fica_rev_nom_t    += add_fica_t
            taxable_pay_nom_t += add_taxable
            ss_outlays_nom_t  += add_outlays_t
        end

        # Benefit taxation
        base_frac_oasi   = max(btax_oasi_rate_t / 0.85, 1e-6)
        creep_oasi       = 1.0 / (1 + (1/base_frac_oasi - 1) *
                                  exp(-thresh_sens * log(cum_price[t])))
        btax_oasi_adj    = (creep_oasi / base_frac_oasi) * btax_oasi_rate_t

        base_frac_hi     = max(btax_hi_rate_t / 0.85, 1e-6)
        creep_hi         = 1.0 / (1 + (1/base_frac_hi - 1) *
                                  exp(-0.6 * thresh_sens * log(cum_price[t])))
        btax_hi_adj      = (creep_hi / base_frac_hi) * btax_hi_rate_t

        ss_btax_oasi_nom_t  = btax_oasi_adj * avg_ben_nom * n_ret_t
        ss_btax_hi_nom_t    = btax_hi_adj   * avg_ben_nom * n_ret_t
        ss_btax_total_nom_t = ss_btax_oasi_nom_t + ss_btax_hi_nom_t

        # Trust fund dynamics
        r_nom_t        = (1.0 + par[:r_G]) * (1.0 + inf_vec[t]) - 1.0
        tf_interest_t  = trust_fund > 0 ? trust_fund * trust_fund_rate : 0.0
        ss_cash_flow_t = fica_rev_nom_t + ss_btax_oasi_nom_t - ss_outlays_nom_t
        ss_total_inc_t = fica_rev_nom_t + ss_btax_oasi_nom_t + tf_interest_t
        ss_balance_t   = ss_total_inc_t - ss_outlays_nom_t
        trust_fund    += ss_balance_t

        # Government budget
        other_tax_nom_t = (inc_tax_mu + corp_tax_mu) * cum_A[t] * mu * cum_price[t] * N_hh_t
        total_tax_nom_t = fica_rev_nom_t + other_tax_nom_t
        G_nom_t         = par[:G_residual_share] * GDP_nom_t

        debt_t = t == 1 ?
            base.D * mu * N_hh * cum_price[t] :
            debt_nom[t-1] * 1e9 * (1 + r_nom_t) +
            (G_nom_t + ss_outlays_nom_t - total_tax_nom_t)

        GDP_real[t]                = GDP_real_t / 1e9
        GDP_nominal[t]             = GDP_nom_t  / 1e9
        avg_earnings_nominal[t]    = avg_earn_t
        payroll_cap_nom[t]         = cap_nom_vec[t]
        taxable_payroll_nom[t]     = taxable_pay_nom_t / 1e9
        fica_revenue_nom[t]        = fica_rev_nom_t    / 1e9
        ss_outlays_nom[t]          = ss_outlays_nom_t  / 1e9
        avg_ben_per_retiree_nom[t] = avg_ben_nom / 1e3
        ss_benefit_tax_oasi_nom[t] = ss_btax_oasi_nom_t  / 1e9
        ss_benefit_tax_hi_nom[t]   = ss_btax_hi_nom_t    / 1e9
        ss_benefit_tax_total_nom[t]= ss_btax_total_nom_t / 1e9
        trust_fund_interest_nom[t] = tf_interest_t / 1e9
        total_pay_nom[t]           = total_pay_nom_t / 1e9
        pct_above_nom[t]           = pct_above_t
        ss_cash_flow_nom[t]        = ss_cash_flow_t / 1e9
        ss_total_income_nom[t]     = ss_total_inc_t / 1e9
        ss_balance_nom[t]          = ss_balance_t   / 1e9
        trust_fund_nom_v[t]        = trust_fund     / 1e9
        fica_rate_eff[t]           = fica_rate_t * 100
        ss_cost_rate[t]            = ss_outlays_nom_t / max(taxable_pay_nom_t,1) * 100
        ss_cash_flow_pct[t]        = ss_cash_flow_t  / max(taxable_pay_nom_t,1) * 100
        ss_balance_pct[t]          = ss_balance_t    / max(taxable_pay_nom_t,1) * 100
        total_tax_rev_nom[t]       = total_tax_nom_t / 1e9
        govt_spending_nom[t]       = G_nom_t / 1e9
        debt_nom[t]                = debt_t  / 1e9
        debt_to_gdp_v[t]           = debt_t  / GDP_nom_t
        mandate_fica_nom[t]        = add_fica_t / 1e9
        mandate_outlays_nom[t]     = add_outlays_t / 1e9
    end

    dep_yr = findall(v -> v < 0, trust_fund_nom_v)
    depl   = isempty(dep_yr) ? nothing : year_cal[dep_yr[1]]
    cf_yr  = findall(v -> v < 0, ss_cash_flow_nom)
    cf_def = isempty(cf_yr)  ? nothing : year_cal[cf_yr[1]]

    # --- Print actuarial summary ---
    windows = filter(w -> w <= n_years, [10, 30, 75])
    println("\n--- Actuarial Summary [$label] ---")
    if !isnothing(cola_inflation)
        @printf("  COLA inflation measure: %.1f%% (headline CPI: %.1f%%)\n",
                cola_inflation*100, inflation*100)
    end
    @printf("  %-8s %10s %10s %10s %12s %10s\n",
            "Window","FICA(\$B)","Outlays(\$B)","Balance(\$B)","TF End(\$B)","CostRate")
    println(repeat("-", 68))
    for w in windows
        @printf("  %-8s %10.1f %10.1f %10.1f %12.1f %9.2f%%\n",
                "$(w)-yr",
                sum(fica_revenue_nom[1:w]),
                sum(ss_outlays_nom[1:w]),
                sum(ss_balance_nom[1:w]),
                trust_fund_nom_v[w],
                mean(ss_cost_rate[1:w]))
    end
    if !isnothing(depl)
        @printf("\n  TF depletion:        %d\n", depl)
    else
        @printf("\n  TF solvent through:  %d\n", start_year+n_years-1)
    end
    !isnothing(cf_def) && @printf("  Cash-flow deficit:   %d\n", cf_def)

    # --- Actuarial score (SSA-style) ---
    disc = [1.0 / (1.0 + trust_fund_rate)^(t-1) for t in 1:n_years]
    pv_cf      = sum(ss_cash_flow_nom   .* disc)
    pv_tp      = sum(taxable_payroll_nom .* disc)
    actuarial_balance_pct = pv_cf / max(pv_tp, 1e-10) * 100.0
    annual_balance_75_pct = ss_cash_flow_pct[n_years]

    change_lab  = nothing; change_ab75  = nothing
    score_lab   = nothing; score_ab75   = nothing
    if !isnothing(baseline_proj)
        change_lab  = actuarial_balance_pct - baseline_proj.actuarial_balance_pct
        change_ab75 = annual_balance_75_pct - baseline_proj.annual_balance_75_pct
        base_lab    = baseline_proj.actuarial_balance_pct
        base_ab75   = baseline_proj.annual_balance_75_pct
        score_lab   = base_lab  < 0.0 ? change_lab  / abs(base_lab)  * 100.0 : nothing
        score_ab75  = base_ab75 < 0.0 ? change_ab75 / abs(base_ab75) * 100.0 : nothing

        println("\n--- SSA Score [$label] ---")
        @printf("  %-35s %9s   %9s\n", "", "Chg(%pay)", "Elim(%)")
        println("  " * repeat("-", 57))
        @printf("  %-35s %+9.2f   %s\n", "Long-range actuarial balance:", change_lab,
                isnothing(score_lab)  ? "  n/a" : @sprintf("%8.1f%%", score_lab))
        @printf("  %-35s %+9.2f   %s\n", "Annual balance in 75th year:", change_ab75,
                isnothing(score_ab75) ? "  n/a" : @sprintf("%8.1f%%", score_ab75))
    else
        println("\n--- SSA Score [$label] ---")
        @printf("  Long-range actuarial balance: %+.2f%% of taxable payroll\n",
                actuarial_balance_pct)
        @printf("  Annual balance in 75th year:  %+.2f%% of taxable payroll\n",
                annual_balance_75_pct)
    end

    # --- Mandatory coverage summary ---
    if has_mandate
        mand_start_idx = findfirst(yr -> yr >= mandate_coverage_year, year_cal)
        if !isnothing(mand_start_idx)
            println("\n--- Mandatory Coverage Impact [$label] ---")
            show_m = sort(unique(filter(t -> t <= n_years && t >= mand_start_idx,
                [mand_start_idx, mand_start_idx+4, mand_start_idx+9,
                 mand_start_idx+19, mand_start_idx+29, mand_start_idx+44, n_years])))
            @printf("  %-6s %10s %10s %10s %8s\n",
                    "Year", "AddFICA\$B", "AddOut\$B", "Net\$B", "PIA_frac")
            println(repeat("-", 50))
            for t in show_m
                tau_m = year_cal[t] - mandate_coverage_year
                cf = min(tau_m, working_years) / working_years
                pf = cf ^ pia_concavity
                @printf("  %-6d %10.1f %10.1f %10.1f %7.1f%%\n",
                        year_cal[t],
                        mandate_fica_nom[t],
                        mandate_outlays_nom[t],
                        mandate_fica_nom[t] - mandate_outlays_nom[t],
                        pf * 100)
            end
            total_net = sum(mandate_fica_nom) - sum(mandate_outlays_nom)
            @printf("  Cumulative net: \$%.1fB over %d years\n",
                    total_net, n_years - mand_start_idx + 1)
        end
    end

    # --- Year-by-year table ---
    show_t = sort(unique(filter(t -> t <= n_years, [1,2,3,5,10,15,20,25,30,40,50,75])))
    println("\n--- Year-by-Year Projection [$label] (selected years) ---")
    @printf("  %-6s %-4s %5s %5s %9s %8s %8s %8s %8s %8s %8s\n",
            "Year","Dep","gA%","CPI%","Cap(\$K)","CostRate",
            "CshFlow","Balance","TF(\$T)","AvgBen\$K","GDP\$T")
    println(repeat("-", 96))
    for t in show_t
        @printf("  %-6d %4.3f %4.1f%% %4.1f%% %9.1f %7.2f%% %7.2f%% %7.2f%% %8.2f %8.1f %8.1f\n",
                year_cal[t], dep_path[t],
                g_A_vec[t]*100, inf_vec[t]*100,
                payroll_cap_nom[t]/1e3,
                ss_cost_rate[t], ss_cash_flow_pct[t], ss_balance_pct[t],
                trust_fund_nom_v[t]/1000,
                avg_ben_per_retiree_nom[t],
                GDP_nominal[t]/1000)
    end

    (year=year_cal, t=1:n_years,
     dep_ratio=dep_path, price_level=cum_price,
     cum_productivity=cum_A, cum_population=cum_pop,
     g_A_annual=g_A_vec, inflation_annual=inf_vec,
     cola_inflation_annual=cola_inf_vec,
     GDP_real=GDP_real, GDP_nominal=GDP_nominal,
     avg_earnings_nominal=avg_earnings_nominal,
     payroll_cap_nom=payroll_cap_nom,
     total_pay_nom=total_pay_nom,
     pct_above_nom=pct_above_nom,
     taxable_payroll_nom=taxable_payroll_nom,
     fica_revenue_nom=fica_revenue_nom,
     ss_outlays_nom=ss_outlays_nom,
     avg_ben_per_retiree_nom=avg_ben_per_retiree_nom,
     ss_benefit_tax_oasi_nom=ss_benefit_tax_oasi_nom,
     ss_benefit_tax_hi_nom=ss_benefit_tax_hi_nom,
     ss_benefit_tax_total_nom=ss_benefit_tax_total_nom,
     trust_fund_interest_nom=trust_fund_interest_nom,
     ss_cash_flow_nom=ss_cash_flow_nom,
     ss_total_income_nom=ss_total_income_nom,
     ss_balance_nom=ss_balance_nom,
     trust_fund_nom=trust_fund_nom_v,
     fica_rate_eff=fica_rate_eff,
     ss_cost_rate=ss_cost_rate,
     ss_cash_flow_pct=ss_cash_flow_pct,
     ss_balance_pct=ss_balance_pct,
     total_tax_rev_nom=total_tax_rev_nom,
     govt_spending_nom=govt_spending_nom,
     debt_nom=debt_nom,
     debt_to_gdp=debt_to_gdp_v,
     mandate_fica_nom=mandate_fica_nom,
     mandate_outlays_nom=mandate_outlays_nom,
     depletion_year=depl,
     cashflow_deficit_year=cf_def,
     actuarial_balance_pct=actuarial_balance_pct,
     annual_balance_75_pct=annual_balance_75_pct,
     change_lab=change_lab,
     change_ab75=change_ab75,
     score_lab=score_lab,
     score_ab75=score_ab75,
     label=label,
     convergence_rate=conv_rate)
end;
proj_reform = project_economy_w_sl(
    ss,
    ss_reform = ss_reform,
    mandate_coverage_year = 2026,
    reform_year = 2026,
    n_years          = 75,
    start_year       = 2025, 
    g_A              = 0.0113,
    g_pop            = 0.005,
    cola_inflation  = 0.021,       # C-CPI-U ≈ 0.3pp lower
    label           = "C-CPI-U COLA",
    inflation        = 0.024,
    dep_path         = dep_path_le,
    ss_cola          = "wage",
    trust_fund_init  = 2.8e12,
    trust_fund_rate  = 0.047,
    baseline_proj    = proj_base,
    gdp_anchor       = 28e12
);
# ============================================================
# EXPORT TO XLSX
# ============================================================

function proj_to_df(p)
    DataFrame(
        year                  = p.year,
        ss_cash_flow_B        = p.ss_cash_flow_nom,
        taxable_payroll_nom_B = p.taxable_payroll_nom,
        ss_outlays_nom_B      = p.ss_outlays_nom,
        fica_revenue_nom_B    = p.fica_revenue_nom,
        ss_cost_rate_pct      = p.ss_cost_rate,
        ss_cash_flow_pct      = p.ss_cash_flow_pct,
        ss_balance_pct        = p.ss_balance_pct,
        trust_fund_nom_B      = p.trust_fund_nom,
        avg_ben_per_retiree_K = p.avg_ben_per_retiree_nom,
        GDP_nominal_B         = p.GDP_nominal,
    )
end

XLSX.openxlsx("reforms_3_5_2026.xlsx", mode="w") do xf
    # Baseline sheet
    sheet = xf[1]
    XLSX.rename!(sheet, "baseline")
    df_base = proj_to_df(proj_base)
    XLSX.writetable!(sheet, df_base)

    # Reform sheet
    XLSX.addsheet!(xf, "reform")
    df_reform = proj_to_df(proj_reform)
    XLSX.writetable!(xf["reform"], df_reform)
end

println("reforms_3_5_2026.xlsx")