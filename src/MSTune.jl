
include("helper_functions.jl")

##################################### import packages #####################################
using CSV
using DataFrames
using LinearAlgebra
using Statistics
using Random
using Base.Threads
using BlackBoxOptim


##################################### parse command-line arguments #####################################
args = parse_args()
query_data = args["--query_data"]
reference_data = args["--reference_data"]
output = args["--output"]
platform = args["--platform"]
optimization_method = args["--optimization_method"]
metric = args["--metric"]
params_to_optimize = args["--params_to_optimize"]
spectrum_preprocessing_order = args["--spectrum_preprocessing_order"]
crossvalidation = parse(Bool, args["--crossvalidation"])
n_folds = parse(Int, args["--n_folds"])
filter_reference_candidates = parse(Bool, args["--filter_reference_candidates"])
bootstrap_query = parse(Bool, args["--bootstrap_query"])
random_seed = parse(Int, args["--random_seed"])
num_threads = parse(Int, args["--threads"])
max_steps = parse(Int, args["--max_steps"])
pop_size = parse(Int, args["--pop_size"])
LB_ws_matching = parse(Float64, args["--LB_ws_matching"])
UB_ws_matching = parse(Float64, args["--UB_ws_matching"])
LB_ws_centroiding = parse(Float64, args["--LB_ws_centroiding"])
UB_ws_centroiding = parse(Float64, args["--UB_ws_centroiding"])
LB_noise_thresh = parse(Float64, args["--LB_noise_thresh"])
UB_noise_thresh = parse(Float64, args["--UB_noise_thresh"])
n_grid_points = parse(Int, args["--n_grid_points"])
LB_LET_thresh = parse(Float64, args["--LB_LET_thresh"])
UB_LET_thresh = parse(Float64, args["--UB_LET_thresh"])
LB_noise_thresh = parse(Float64, args["--LB_noise_thresh"])
UB_noise_thresh = parse(Float64, args["--UB_noise_thresh"])
LB_wf_mz = parse(Float64, args["--LB_wf_mz"])
UB_wf_mz = parse(Float64, args["--UB_wf_mz"])
LB_wf_int = parse(Float64, args["--LB_wf_int"])
UB_wf_int = parse(Float64, args["--UB_wf_int"])
ws_matching = parse(Float64, args["--ws_matching"])
ws_centroiding = parse(Float64, args["--ws_centroiding"])
LET_thresh = parse(Float64, args["--LET_thresh"])
noise_thresh = parse(Float64, args["--noise_thresh"])
wf_int = parse(Float64, args["--wf_int"])
wf_mz = parse(Float64, args["--wf_mz"])

Random.seed!(random_seed)
rng = MersenneTwister(random_seed)
const K = n_folds

if nthreads() < num_threads
    println("Warning: Julia was started with $(nthreads()) threads, cannot increase to $num_threads")
end


if platform == "NRMS"
    const PARAM_NAMES = ["LET_thresh", "noise_thresh", "wf_int", "wf_mz"]
else
    const PARAM_NAMES = ["LET_thresh", "noise_thresh", "wf_int", "wf_mz", "ws_matching", "ws_centroiding"]
end

function parse_params_to_optimize(s::AbstractString, platform)
    s = strip(s)
    if s == "all"
        if platform == "NRMS"
            return collect(1:4)
        else
            return collect(1:6)
        end
    end
    parts = split(s, ',')
    parts = strip.(parts)
    if platform == "NRMS"
        name_to_idx = Dict("LET_thresh" => 1, "noise_thresh" => 2, "wf_int" => 3, "wf_mz" => 4)
    else 
        name_to_idx = Dict("LET_thresh" => 1, "noise_thresh" => 2, "wf_int" => 3, "wf_mz" => 4, "ws_matching" => 5, "ws_centroiding" => 6)
    end
    idxs = Int[]
    for p in parts
        haskey(name_to_idx, p) || error("Unknown param in --params_to_optimize: '$p'. Allowed: all, LET_thresh, noise_thresh, wf_int, wf_mz, ws_matching, ws_centroiding (comma-separated).")
        push!(idxs, name_to_idx[p])
    end
    idxs = unique(idxs)
    sort!(idxs)
    return idxs
end

opt_idxs = parse_params_to_optimize(params_to_optimize, platform)

if platform == "NRMS"
    p_fixed_full = Float64[LET_thresh, noise_thresh, wf_int, wf_mz]
    LB_full = Float64[LB_LET_thresh, LB_noise_thresh, LB_wf_int, LB_wf_mz]
    UB_full = Float64[UB_LET_thresh, UB_noise_thresh, UB_wf_int, UB_wf_mz]
else
    p_fixed_full = Float64[LET_thresh, noise_thresh, wf_int, wf_mz, ws_matching, ws_centroiding]
    LB_full = Float64[LB_LET_thresh, LB_noise_thresh, LB_wf_int, LB_wf_mz, ws_matching, ws_centroiding]
    UB_full = Float64[UB_LET_thresh, UB_noise_thresh, UB_wf_int, UB_wf_mz, ws_matching, ws_centroiding]
end
bounds_free = [(LB_full[i], UB_full[i]) for i in opt_idxs]

function expand_params(p_free::AbstractVector{<:Real}, p_fixed_full::Vector{Float64}, opt_idxs::Vector{Int})
    p = copy(p_fixed_full)
    @assert length(p_free) == length(opt_idxs)
    for (j, idx) in enumerate(opt_idxs)
        p[idx] = Float64(p_free[j])
    end
    return p
end


##################################### import data #####################################
df_query_raw_tmp = CSV.read(query_data, DataFrame; delim='\t')

if platform == "NRMS"
    if bootstrap_query
        n_tmp = nrow(df_query_raw_tmp)
        idxs = sample(1:n_tmp, n_tmp; replace=true)
        df_query_raw = df_query_raw_tmp[idxs,:]
    else
        df_query_raw = df_query_raw_tmp
    end
elseif platform == "HRMS"
    if bootstrap_query
        spectrum_col = :id
        spec_ids = unique(df_query_raw_tmp[!, spectrum_col])
        n_spec = length(spec_ids)
        sampled_spec_ids = sample(spec_ids, n_spec; replace=true)
        boot_blocks = DataFrame[]
        for (b, sid) in enumerate(sampled_spec_ids)
            block = copy(df_query_raw_tmp[df_query_raw_tmp[!, spectrum_col] .== sid, :])
            block[!, :bootstrap_spectrum_id] .= b
            push!(boot_blocks, block)
        end
        df_query_raw = vcat(boot_blocks...)
    else
        df_query_raw = df_query_raw_tmp
    end
end

df_ref_raw = CSV.read(reference_data, DataFrame; delim='\t', pool=false, missingstring=["","NA","N/A"], stringtype=String)
q_ids_all = string.(df_query_raw[!, :id])
r_ids_all = string.(df_ref_raw[!, :id])
println("Note: There are $(length(setdiff(unique(q_ids_all),unique(r_ids_all)))) query spectra whose ground-truth ID is not contained in the reference library.")


##################################### run analysis #####################################
callback = oc -> begin
    push!(fitness_progress_history, (evals = BlackBoxOptim.num_func_evals(oc), fitness = BlackBoxOptim.best_fitness(oc), params = copy(BlackBoxOptim.best_candidate(oc))))
    false
end

if platform == "NRMS"
    Q0 = Matrix{Float64}(df_query_raw[:, 2:end])
    R0 = Matrix{Float64}(df_ref_raw[:, 2:end])
    m = size(Q0, 2)
    @assert m == size(R0, 2) "Query and reference must have same number of columns"
    mzs = collect(1.0:m)

    if crossvalidation
        folds = make_folds_NRMS(size(Q0,1), K; rng=rng)
        if metric == "accuracy"
            objective_full = objective_acc_cv_NRMS
        else
            objective_full = objective_MRR_cv_NRMS
        end
    else
        if metric == "accuracy"
            objective_full = objective_acc_no_cv_NRMS
        else
            objective_full = objective_MRR_no_cv_NRMS
        end
    end

    objective_free(p_free) = objective_full(expand_params(p_free, p_fixed_full, opt_idxs))

    if optimization_method == "DE"
        fitness_progress_history = Vector{NamedTuple{(:evals, :fitness, :params_full), Tuple{Int, Float64, Vector{Float64}}}}()
        if params_to_optimize == "all"
            initial_guess = [LET_thresh, noise_thresh, wf_int, wf_mz]
        elseif params_to_optimize == "LET_thresh"
            initial_guess = [LET_thresh]
        elseif params_to_optimize == "noise_thresh"
            initial_guess = [noise_thresh]
        elseif params_to_optimize == "wf_int"
            initial_guess = [wf_int]
        elseif params_to_optimize == "wf_mz"
            initial_guess = [wf_mz]
        end

        callback = oc -> begin
            p_free_best = BlackBoxOptim.best_candidate(oc)
            p_full_best = expand_params(p_free_best, p_fixed_full, opt_idxs)
            push!(fitness_progress_history, (evals=BlackBoxOptim.num_func_evals(oc), fitness=BlackBoxOptim.best_fitness(oc), params_full=p_full_best))
            false
        end
        res = bboptimize(objective_free, initial_guess; SearchRange=bounds_free, NumDimensions=length(bounds_free), Method=:de_rand_1_bin, PopulationSize=pop_size, MaxSteps=max_steps, CallbackFunction=callback, CallbackInterval=0.0, rng=rng)

    elseif optimization_method == "grid"
        grid_full = [
            collect(range(LB_full[1], UB_full[1]; length=n_grid_points)),
            collect(range(LB_full[2], UB_full[2]; length=n_grid_points)),
            collect(range(LB_full[3], UB_full[3]; length=n_grid_points)),
            collect(range(LB_full[4], UB_full[4]; length=n_grid_points))
        ]
        grids_free = [grid_full[i] for i in opt_idxs]
        combos = collect(Iterators.product(grids_free...))
        n = length(combos)
        LET_out = Vector{Float64}(undef, n)
        noise_out = Vector{Float64}(undef, n)
        wf_int_out = Vector{Float64}(undef, n)
        wf_mz_out = Vector{Float64}(undef, n)
        fit_out = Vector{Float64}(undef, n)

        Threads.@threads for i in 1:n
            p_free = collect(combos[i])
            p_full = expand_params(p_free, p_fixed_full, opt_idxs)
            fit = objective_full(p_full)
            LET_out[i] = p_full[1]
            noise_out[i] = p_full[2]
            wf_int_out[i] = p_full[3]
            wf_mz_out[i] = p_full[4]
            fit_out[i] = fit
        end
        grid_results_tmp = DataFrame(LET_thresh=LET_out, noise_thresh=noise_out, wf_int=wf_int_out, wf_mz=wf_mz_out, fitness=fit_out)
        grid_results = sort!(grid_results_tmp, :fitness, rev=false)
    elseif optimization_method == "none"
        S = get_scores_NRMS(Q0, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
        df_scores = DataFrame(S, r_ids_all)
        insertcols!(df_scores, 1, :query_id => q_ids_all)
    end
else
    if crossvalidation
        folds = make_folds_HRMS(unique(q_ids_all), K; rng=rng)
        if metric == "accuracy"
            if filter_reference_candidates
                objective_full = objective_cv_acc_filtering_HRMS
            else
                objective_full = objective_cv_acc_no_filtering_HRMS
            end
        elseif metric == "MRR"
            if filter_reference_candidates
                objective_full = objective_cv_MRR_filtering_HRMS
            else
                objective_full = objective_cv_MRR_no_filtering_HRMS
            end
        end
    else
        if metric == "accuracy"
            if filter_reference_candidates
                objective_full = objective_no_cv_acc_filtering_HRMS
            else
                objective_full = objective_no_cv_acc_no_filtering_HRMS
            end
        elseif metric == "MRR"
            if filter_reference_candidates
                objective_full = objective_no_cv_MRR_filtering_HRMS
            else
                objective_full = objective_no_cv_MRR_no_filtering_HRMS
            end
        end
    end

    if optimization_method == "DE"
        if params_to_optimize == "all"
            initial_guess = [LET_thresh, noise_thresh, wf_int, wf_mz, ws_matching, ws_centroiding]
        elseif params_to_optimize == "LET_thresh"
            initial_guess = [LET_thresh]
        elseif params_to_optimize == "noise_thresh"
            initial_guess = [noise_thresh]
        elseif params_to_optimize == "wf_int"
            initial_guess = [wf_int]
        elseif params_to_optimize == "wf_mz"
            initial_guess = [wf_mz]
        elseif params_to_optimize == "ws_matching"
            initial_guess = [ws_matching]
        elseif params_to_optimize == "ws_centroiding"
            initial_guess = [ws_centroiding]
        end

        fitness_progress_history = Vector{NamedTuple{(:evals, :fitness, :params_full), Tuple{Int, Float64, Vector{Float64}}}}()
        objective_free(p_free) = objective_full(expand_params(p_free, p_fixed_full, opt_idxs))
        callback = oc -> begin
            p_free_best = BlackBoxOptim.best_candidate(oc)
            p_full_best = expand_params(p_free_best, p_fixed_full, opt_idxs)
            push!(fitness_progress_history, (evals=BlackBoxOptim.num_func_evals(oc), fitness=BlackBoxOptim.best_fitness(oc), params_full=p_full_best))
            false
        end
        res = bboptimize(objective_free, initial_guess; SearchRange=bounds_free, NumDimensions=length(bounds_free), Method=:de_rand_1_bin, PopulationSize=pop_size, MaxSteps=max_steps, CallbackFunction=callback, CallbackInterval=0.0, rng=rng)

    end
    if optimization_method == "grid"
        grid_LET = collect(range(LB_LET_thresh, UB_LET_thresh; length=n_grid_points))
        grid_noise = collect(range(LB_noise_thresh, UB_noise_thresh; length=n_grid_points))
        grid_wf_int = collect(range(LB_wf_int, UB_wf_int; length=n_grid_points))
        grid_wf_mz = collect(range(LB_wf_mz, UB_wf_mz; length=n_grid_points))
        grid_ws_matching = collect(range(LB_ws_matching, UB_ws_matching; length=n_grid_points))
        grid_ws_centroiding = collect(range(LB_ws_centroiding, UB_ws_centroiding; length=n_grid_points))
        combos = collect(Iterators.product(grid_LET, grid_noise, grid_wf_int, grid_wf_mz, grid_ws_matching, grid_ws_centroiding))
        n = length(combos)
        LET_out = Vector{Float64}(undef, n)
        noise_out = Vector{Float64}(undef, n)
        wf_int_out = Vector{Float64}(undef, n)
        wf_mz_out = Vector{Float64}(undef, n)
        ws_matching_out = Vector{Float64}(undef, n)
        ws_centroiding_out = Vector{Float64}(undef, n)
        obj_out = Vector{Float64}(undef, n)
        Threads.@threads for i in 1:n
            LET, noise, wf_int, wf_mz, ws_matching, ws_centroiding = combos[i]
            x = (LET, noise, wf_int, wf_mz, ws_matching, ws_centroiding)
            acc = objective_full(collect(x))
            LET_out[i] = LET
            noise_out[i] = noise
            wf_int_out[i] = wf_int
            wf_mz_out[i] = wf_mz
            ws_matching_out[i] = ws_matching
            ws_centroiding_out[i] = ws_centroiding
            obj_out[i] = acc
        end
        grid_results_tmp = DataFrame(LET_thresh=LET_out, noise_thresh=noise_out, wf_int=wf_int_out, wf_mz=wf_mz_out, ws_matching=ws_matching_out, ws_centroiding=ws_centroiding_out, fitness=obj_out, )
        grid_results = sort!(grid_results_tmp, :fitness, rev=false)
    elseif optimization_method == "none"
        df_scores = get_scores_HRMS(df_query_raw, df_ref_raw; order=spectrum_preprocessing_order, LET_threshold=LET_thresh, noise_threshold=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, window_size_matching=ws_matching, window_size_centroiding=ws_centroiding)
    end
end


if optimization_method == "DE"
    open(output, "w") do io
        if platform == "NRMS"
            header = ["step", "evals", "fitness", "LET_thresh", "noise_thresh", "wf_int", "wf_mz"]
        else
            header = ["step", "evals", "fitness", "LET_thresh", "noise_thresh", "wf_int", "wf_mz", "ws_matching", "ws_centroiding"]
        end
        write(io, join(header, '\t') * "\n")
        for (step, rec) in enumerate(fitness_progress_history)
            if platform == "NRMS"
                row = String[string(step-1), string(rec.evals), string(rec.fitness), string(rec.params_full[1]), string(rec.params_full[2]), string(rec.params_full[3]), string(rec.params_full[4])]
            else
                row = String[string(step-1), string(rec.evals), string(rec.fitness), string(rec.params_full[1]), string(rec.params_full[2]), string(rec.params_full[3]), string(rec.params_full[4]), string(rec.params_full[5]), string(rec.params_full[6])]
            end
            write(io, join(row, '\t') * "\n")
        end
    end
elseif optimization_method == "grid"
    CSV.write(output, grid_results; delim='\t')
elseif optimization_method == "none"
    CSV.write(output, df_scores; delim='\t')
end


