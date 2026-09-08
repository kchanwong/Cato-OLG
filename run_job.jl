# AWS Batch entrypoint. One job = one score.
#
#   julia --threads=auto run_job.jl '<json>'
#   julia --threads=auto run_job.jl /path/to/job.json
#   julia --threads=auto run_job.jl s3://bucket/key.json
#
# The JSON object holds:
#   label           The name of the run.               (default "Reform")
#   reform_year     The year in which the reform starts. Needed for a reform.
#   params          Keys of create_params(). Refer to docs/PARAMETERS.md #1.
#   benefit_fn      A benefit rule, by name.           (default current law)
#   aime_modifier   "caregiver_modifier", or null.
#   projection      Keywords of project_economy and static_score_economy.
#                   Refer to docs/PARAMETERS.md #2. dep_path defaults to the
#                   cached dep_fit.
#   static          true to add a static score.        (default true)
#
# With no params, no benefit_fn and no aime_modifier, the job runs the baseline
# only. Results go to $OUTPUT_S3_URI, or to $OUTPUT_DIR (default /data).

include(joinpath(@__DIR__, "source", "setup.jl"))
using JSON, DataFrames, XLSX

const BENEFIT_FNS = Dict(string(f) => f for f in (
    benefits_current_law, benefits_capped_spousal, benefits_earnings_sharing,
    benefits_caregiver_credits, benefits_price_index_floor))

const AIME_MODIFIERS = Dict("caregiver_modifier" => caregiver_modifier)

# Keywords that only one of the two projection functions accepts.
const PROJECT_ONLY = (:cola_cap, :chained_cpi, :fpl_single_2025, :cola_cap_pct)
const STATIC_ONLY  = (:eliminates_cap,)
const INT_KEYWORDS = (:n_years, :start_year, :reform_year, :year_start, :year_end)

lookup(d, name, what) = get(d, name) do
    error("Unknown $what \"$name\". Use one of: $(join(sort(collect(keys(d))), ", "))")
end

# JSON has no Inf and no distinction between 70 and 70.0. Give each value the
# type that the model expects.
coerce(::Bool, v) = Bool(v)
coerce(::Integer, v) = isinteger(v) ? Int(v) : Float64(v)   # :aime_nwc_frac is an Int 0.
coerce(::AbstractFloat, v) = v isa AbstractString ? parse(Float64, v) : Float64(v)
coerce(::AbstractVector, v) = Float64[Float64(x) for x in v]
coerce(_, v) = v                       # Strings, nothing, functions.

coerce_keyword(k, v::Bool) = v
coerce_keyword(k, v::Real) = k in INT_KEYWORDS ? Int(v) : Float64(v)
coerce_keyword(k, v::AbstractVector) = k == :scenario_periods ?
    [Dict(Symbol(kk) => coerce_keyword(Symbol(kk), vv) for (kk, vv) in p) for p in v] :
    Float64[Float64(x) for x in v]
coerce_keyword(k, v::AbstractString) = v == "Inf" ? Inf : v
coerce_keyword(k, v) = v

function build_params(spec::AbstractDict)
    par = create_params()
    for (k, v) in spec
        key = Symbol(k)
        haskey(par, key) || error("Unknown parameter :$key. Refer to docs/PARAMETERS.md.")
        par[key] = coerce(par[key], v)
    end
    par
end

function projection_kwargs(spec::AbstractDict)
    kw = Dict{Symbol,Any}(:dep_path => dep_fit)
    for (k, v) in spec
        key = Symbol(k)
        kw[key] = coerce_keyword(key, v)
    end
    kw
end

drop(kw, keys) = Dict(k => v for (k, v) in kw if !(k in keys))

function read_job(arg::AbstractString)
    text = if startswith(arg, "s3://")
        read(`aws s3 cp $arg -`, String)
    elseif startswith(strip(arg), "{")
        arg
    else
        read(arg, String)
    end
    JSON.parse(text)
end

proj_to_df(p) = DataFrame(
    year                    = collect(p.year),
    taxable_payroll_nom     = p.taxable_payroll_nom,
    fica_revenue_nom        = p.fica_revenue_nom,
    ss_outlays_nom          = p.ss_outlays_nom,
    ss_cash_flow_nom        = p.ss_cash_flow_nom,
    ss_cost_rate            = p.ss_cost_rate,
    ss_cash_flow_pct        = p.ss_cash_flow_pct,
    ss_balance_pct          = p.ss_balance_pct,
    trust_fund_nom          = p.trust_fund_nom,
    avg_ben_per_retiree_nom = p.avg_ben_per_retiree_nom,
    GDP_nominal             = p.GDP_nominal,
)

summarize(p) = Dict(
    "label"                 => p.label,
    "actuarial_balance_pct" => p.actuarial_balance_pct,
    "annual_balance_75_pct" => p.annual_balance_75_pct,
    "change_lab"            => p.change_lab,
    "change_ab75"           => p.change_ab75,
    "score_lab"             => p.score_lab,
    "score_ab75"            => p.score_ab75,
    "depletion_year"        => p.depletion_year,
    "cashflow_deficit_year" => p.cashflow_deficit_year,
)

function write_results(runs, outdir, name)
    mkpath(outdir)
    xlsx = joinpath(outdir, "$name.xlsx")
    XLSX.openxlsx(xlsx, mode = "w") do xf
        for (sheet_name, p) in runs
            sheet = XLSX.addsheet!(xf, sheet_name)
            df = proj_to_df(p)
            for (c, col) in enumerate(names(df)); sheet[1, c] = col; end
            for (r, row) in enumerate(eachrow(df)), (c, v) in enumerate(row)
                sheet[r + 1, c] = v
            end
        end
    end
    summary = joinpath(outdir, "$name.json")
    open(summary, "w") do io
        JSON.print(io, Dict(k => summarize(p) for (k, p) in runs), 2)
    end
    println("wrote $xlsx and $summary")

    uri = get(ENV, "OUTPUT_S3_URI", "")
    if !isempty(uri)
        dest = endswith(uri, "/") ? uri : uri * "/"
        for f in (xlsx, summary)
            run(`aws s3 cp $f $(dest * basename(f))`)
        end
    end
end

function main(arg::AbstractString)
    job = read_job(arg)
    label       = get(job, "label", "Reform")
    reform_year = get(job, "reform_year", nothing)
    par_spec    = get(job, "params", Dict())
    benefit_fn  = lookup(BENEFIT_FNS, get(job, "benefit_fn", "benefits_current_law"),
                         "benefit rule")
    modifier    = let m = get(job, "aime_modifier", nothing)
        isnothing(m) ? nothing : lookup(AIME_MODIFIERS, m, "AIME modifier")
    end
    kw = projection_kwargs(get(job, "projection", Dict()))

    proj_base = project_economy(ss; drop(kw, STATIC_ONLY)..., label = "Current Law")
    runs = Pair{String,Any}["Baseline" => proj_base]

    is_reform = !isempty(par_spec) || benefit_fn !== benefits_current_law ||
                !isnothing(modifier)
    if is_reform
        isnothing(reform_year) && error("A reform needs \"reform_year\".")
        ss_reform = solve_steady_state(build_params(par_spec);
                                       benefit_fn = benefit_fn, aime_modifier = modifier)

        b, r = extract_ss_ratios(ss), extract_ss_ratios(ss_reform; benefit_fn = benefit_fn)
        @printf("benefit/retiree: %.4f -> %.4f   taxable payroll: %.4f -> %.4f\n",
                b.ss_ben_per_retiree, r.ss_ben_per_retiree,
                b.taxable_payroll,    r.taxable_payroll)

        reform_kw = (ss_reform = ss_reform, reform_year = Int(reform_year),
                     baseline_proj = proj_base)
        push!(runs, "Dynamic" => project_economy(ss; drop(kw, STATIC_ONLY)...,
                                                 reform_kw..., label = label))
        if get(job, "static", true)
            push!(runs, "Static" => static_score_economy(ss; drop(kw, PROJECT_ONLY)...,
                                                         reform_kw...,
                                                         label = "$label (static)"))
        end
    end

    write_results(runs, get(ENV, "OUTPUT_DIR", "/data"),
                  get(job, "output_prefix", "score"))
end

function selftest()
    par = build_params(Dict("J_retire" => 70, "payroll_rate" => 0.124,
                            "payroll_cap" => "Inf", "p_perm_vals" => [1, 0, 0, 0]))
    @assert par[:J_retire] === 70
    @assert par[:payroll_rate] === 0.124
    @assert par[:payroll_cap] === Inf
    @assert par[:p_perm_vals] == [1.0, 0.0, 0.0, 0.0]
    @assert par[:J_start] === create_params()[:J_start]     # Untouched keys keep defaults.

    kw = projection_kwargs(Dict("n_years" => 75, "g_A" => 1e-2, "cola_cap" => true,
                                "ss_cola" => "wage", "trust_fund_init" => 28_000,
                                "scenario_periods" => [Dict("year_start" => 2030,
                                                            "year_end" => 2040,
                                                            "g_A" => 0.005)]))
    @assert kw[:n_years] === 75 && kw[:g_A] === 0.01
    @assert kw[:cola_cap] === true && kw[:ss_cola] == "wage"
    @assert kw[:trust_fund_init] === 28000.0
    @assert kw[:scenario_periods][1][:year_start] === 2030
    @assert kw[:scenario_periods][1][:g_A] === 0.005
    @assert kw[:dep_path] === dep_fit                       # Cached path by default.
    @assert !haskey(drop(kw, PROJECT_ONLY), :cola_cap)

    @assert read_job("{\"label\": \"x\"}")["label"] == "x"
    try; build_params(Dict("payroll_tax_rate" => 0.1)); @assert false
    catch e; @assert occursin("Unknown parameter", sprint(showerror, e)); end
    try; lookup(BENEFIT_FNS, "nope", "benefit rule"); @assert false
    catch e; @assert occursin("Unknown benefit rule", sprint(showerror, e)); end
    println("selftest ok")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Give one argument: JSON, a file path, or an s3:// URI.")
    ARGS[1] == "--selftest" ? selftest() : main(ARGS[1])
end
