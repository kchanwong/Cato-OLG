# Import Function 
include(joinpath(@__DIR__, "..", "..", "source", "setup.jl"))
# Packages Needed
using DataFrames
using Plots
using XLSX
using CSV
using JLD2;
# Solve steady state for the reform scenario
par_reform = create_params()
par_reform[:first_rr] = 1.02
par_reform[:second_rr] = 0
par_reform[:third_rr] = 0
ss_reform  = solve_steady_state(par_reform);
# Custom Project Economy
function pia_factors(yr)
    yr < 2029 && return (0.90, 0.32, 0.15)
    yr >= 2036 && return (1.02, 0.00, 0.00)
    t = yr - 2029
    f1 = min(1.02, 0.90 + 0.02*t)
    f2 = max(0.00, 0.32 - 0.04*t)
    f3 = max(0.00, 0.15 - 0.03*t)
    (f1, f2, f3)
end

function project_economy_flat_benefits(
    ss;
    n_years::Int          = 75,
    start_year::Int       = 2025,
    g_A::Float64          = 0.0114,
    g_pop::Float64        = 0.005,
    inflation::Float64    = 0.024,
    dep_path              = nothing,
    ss_cola               = "wage",
    cap_base_dollars::Float64 = 176100.0,
    bracket_indexing::String  = "cpi",
    trust_fund_init::Float64  = 2.8e12,
    trust_fund_rate::Float64  = 0.047,
    ss_reform             = nothing,
    reform_year           = nothing,
    label::String         = "Current Law",
    scenario_periods      = nothing,
    thresh_sens::Float64  = 0.3,
    ccpiu_wedge::Float64  = 0.003,
    ccpiu_start_year::Int = 2036,
    pia_factor_indexing::Bool = false,
    gdp_anchor            = nothing,  # if set, rescale N_hh so year-1 GDP = gdp_anchor
    wages_to_gdp::Float64 = 0.46      # wages & salaries as share of GDP (empirical calibration)
)
    @printf("=== Projecting %d years: %d-%d [%s] ===\n",
            n_years, start_year, start_year+n_years-1, label)

    base = extract_ss_ratios(ss)
    mu   = base.mu_dollar; par = ss.par
    alpha = par[:alpha]

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

    # --- Steady-state aggregates for transition path ---
    K_base = ss.K;  L_base = ss.L
    K_ref  = 0.0;   L_ref  = 0.0

    has_reform = !isnothing(ss_reform) && !isnothing(reform_year)
    ref = nothing
    conv_rate = 0.0
    if has_reform
        ref = extract_ss_ratios(ss_reform)
        K_ref = ss_reform.K;  L_ref = ss_reform.L
        @printf("  Reform at %d\n", reform_year)
        conv_rate = compute_convergence_rate(ss, g_A=g_A)
    end
    if ccpiu_wedge > 0
        @printf("  C-CPI-U indexing from %d (wedge=%.1f bp)\n",
                ccpiu_start_year, ccpiu_wedge*10000)
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

    # Apply C-CPI-U wedge: slower new-beneficiary indexing after ccpiu_start_year
    if ccpiu_wedge > 0
        cum_new_ben_adj = copy(cum_new_ben)
        for t in 2:n_years
            if year_cal[t] >= ccpiu_start_year
                base_growth = cum_new_ben[t] / cum_new_ben[t-1]
                adj_growth  = base_growth / (1 + ccpiu_wedge)
                cum_new_ben_adj[t] = cum_new_ben_adj[t-1] * adj_growth
            else
                cum_new_ben_adj[t] = cum_new_ben[t]
            end
        end
        cum_new_ben = cum_new_ben_adj
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

    # --- Precompute bend points and average AIME ---
    b1_mu = par[:ss_bend1] / mu
    b2_mu = par[:ss_bend2] / mu

    avg_aime_mu = b2_mu
    pia_base_level = base.ss_ben_per_retiree
    for _ in 1:20
        pia_try = 0.90*min(avg_aime_mu, b1_mu) +
                  0.32*max(0.0, min(avg_aime_mu, b2_mu) - b1_mu) +
                  0.15*max(0.0, avg_aime_mu - b2_mu)
        err = pia_try - pia_base_level
        mr = avg_aime_mu <= b1_mu ? 0.90 :
             avg_aime_mu <= b2_mu ? 0.32 : 0.15
        avg_aime_mu -= err / max(mr, 0.01)
        avg_aime_mu = max(avg_aime_mu, 0.01)
    end

    pia_current_law = 0.90*min(avg_aime_mu, b1_mu) +
                      0.32*max(0.0, min(avg_aime_mu, b2_mu) - b1_mu) +
                      0.15*max(0.0, avg_aime_mu - b2_mu)

    @printf("  Avg AIME (model units): %.4f  PIA: %.4f\n", avg_aime_mu, pia_current_law)

    # Output storage
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

    trust_fund   = trust_fund_init
    avg_ben_real = base.ss_ben_per_retiree

    for t in 1:n_years
        yr  = year_cal[t]

        # === Mankiw-Weinzierl transition path ===
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

        # Derive macro aggregates from production function
        Y_mu = K_mu^alpha * L_mu^(1-alpha)
        w_mu = (1-alpha) * (K_mu / L_mu)^alpha

        # Interpolate distributional variables
        interp(bv, rv) = bv + omega * (rv - bv)
        fica_rate_t      = has_reform ? interp(base.payroll_rate,   ref.payroll_rate)   : base.payroll_rate
        pct_above_t      = has_reform ? interp(base.pct_above_cap,  ref.pct_above_cap)  : base.pct_above_cap
        inc_tax_mu       = has_reform ? interp(base.income_tax,     ref.income_tax)     : base.income_tax
        btax_oasi_rate_t = has_reform ? interp(btax_rate_oasi_base, btax_rate_oasi_ref) : btax_rate_oasi_base
        btax_hi_rate_t   = has_reform ? interp(btax_rate_hi_base,   btax_rate_hi_ref)   : btax_rate_hi_base

        # Corporate tax from transition path
        zeta_corp = par[:zeta_income_corp]
        corp_profit_mu = zeta_corp * (Y_mu - w_mu * L_mu - par[:delta] * K_mu)
        corp_tax_mu = par[:tau_statutory_corp] *
                      max(0.0, corp_profit_mu) *
                      (1 - par[:phi_exp_corp]*0.3 - par[:zeta_ded_corp])

        dep_t     = dep_path[t]

        # --- Year-specific PIA factors ---
        if pia_factor_indexing && !isnothing(reform_year) && yr >= reform_year
            t_ref = findfirst(==(reform_year), year_cal)
            price_wage_ratio = isnothing(t_ref) ? 1.0 : cum_A[t_ref] / cum_A[t]
            f1 = 0.90 * price_wage_ratio
            f2 = 0.32 * price_wage_ratio
            f3 = 0.15 * price_wage_ratio
        else
            f1, f2, f3 = pia_factors(yr)
        end

        # Bend points: AWI-indexed before C-CPI-U switch, CPI-indexed after
        if ccpiu_wedge > 0 && yr >= ccpiu_start_year
            b1_t = b1_mu / cum_A[t]
            b2_t = b2_mu / cum_A[t]
        else
            b1_t = b1_mu
            b2_t = b2_mu
        end

        pia_yr = f1*min(avg_aime_mu, b1_t) +
                 f2*max(0.0, min(avg_aime_mu, b2_t) - b1_t) +
                 f3*max(0.0, avg_aime_mu - b2_t)
        new_ben_scale_t = pia_yr / max(pia_current_law, 1e-10)

        # === Scale to nominal dollars ===
        N_hh_t     = N_hh * cum_pop[t]
        GDP_real_t = Y_mu * cum_A[t] * mu * N_hh_t
        GDP_nom_t  = GDP_real_t * cum_price[t]
        avg_earn_t = (w_mu * L_mu / max(base.worker_mass, 1e-10)) *
                     cum_A[t] * mu * cum_price[t]

        # Payroll from transition path
        total_pay_nom_t   = wages_to_gdp * GDP_nom_t
        wage_cap_ratio    = cum_wage[t] / (cap_nom_vec[t] / cap_base_dollars)
        pct_above_adj     = max(0, pct_above_t * wage_cap_ratio)
        taxable_pay_nom_t = total_pay_nom_t * (1.0 - pct_above_adj)
        fica_rev_nom_t    = fica_rate_t * taxable_pay_nom_t

        # SS outlays
        phi_t          = max(0, 1.0 / (working_years * dep_t))
        new_ben_real_t = base.ss_ben_per_retiree * cum_new_ben[t] * new_ben_scale_t
        avg_ben_real   = (1.0-phi_t)*avg_ben_real + phi_t*new_ben_real_t
        avg_ben_nom    = avg_ben_real * mu * cum_price_lag[t]
        n_ret_t        = dep_t * ssa_covered_workers * cum_pop[t]
        ss_outlays_nom_t = avg_ben_nom * n_ret_t

        # Benefit taxation with bracket creep
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

        # Trust fund
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

        # Store
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
    end

    dep_yr = findall(v -> v < 0, trust_fund_nom_v)
    depl   = isempty(dep_yr) ? nothing : year_cal[dep_yr[1]]
    cf_yr  = findall(v -> v < 0, ss_cash_flow_nom)
    cf_def = isempty(cf_yr)  ? nothing : year_cal[cf_yr[1]]

    # --- Print actuarial summary ---
    windows = filter(w -> w <= n_years, [10, 30, 75])
    println("\n--- Actuarial Summary [$label] ---")
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
     depletion_year=depl,
     cashflow_deficit_year=cf_def,
     convergence_rate=conv_rate,
     label=label)
end

###
proj = project_economy_flat_benefits(ss,
    ss_reform       = ss_reform,
    reform_year     = 2029,
    n_years         = 75,
    start_year      = 2025,
    gdp_anchor      = 28e12,
    g_A             = 0.0114,
    g_pop           = 0.005,
    inflation       = 0.024,
    dep_path        = dep_fit,
    ss_cola         = "wage",
    trust_fund_init = 2.8e12,
    trust_fund_rate = 0.047,
    ccpiu_wedge     = 0.003,
    ccpiu_start_year = 2036,
    label="Flat Benefit (102/0/0)");