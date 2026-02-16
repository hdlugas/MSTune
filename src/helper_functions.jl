
using Base.Threads
using DataFrames
using Random
using StatsBase

function print_help()
    println("""
Help Message for OptiMS.jl

Usage:
  julia OptiMS.jl \
    --query_data <string> \
    --reference_data <string> \
    --output <string> \
    --optimization_method <string> \
    --metric <string> \
    --crossvalidation <boolean> \
    --n_folds <int> \
    --bootstrap_query <boolean> \
    --random_seed <int> \
    --params_to_optimize <string> \
    --spectrum_preprocessing_order <string> \
    --LB_LET_thresh <float> \
    --UB_LET_thresh <float> \
    --LB_noise_thresh <float> \
    --UB_noise_thresh <float> \
    --LB_wf_int <float> \
    --UB_wf_int <float> \
    --LB_wf_mz <float> \
    --UB_wf_mz <float> \
    --LET_thresh <float> \
    --noise_thresh <float> \
    --wf_int <float> \
    --wf_mz <float> \
    --threads <int> \
    --max_steps <int>
    --pop_size <int>
    --n_grid_points <int> \


Arguments:
  --query_data                    Path to input TXT file of query dataset (required).
  --reference_data                Path to input TXT file of reference dataset (required).
  --output                        Path to output TXT file (required).
  --optimization_method           Optimization approach (optional, options=[DE,grid,none], default=DE).
  --metric                        Quantity to maximize in the objective function (optional, options=[accuracy,MRR], default=accuracy).
  --crossvalidation               Boolean indicating whether or not to perform 5-fold cross-validation inside the objective function (optional, options=[true,false], default=true).
  --n_folds                       Number folds to use for cross-validation. Only applicable for crossvalidation is \'true\' and optimization_method is not \'none\' (optional, default=5).
  --bootstrap_query               Boolean indicating whether or not to construct query dataset of same size from resampling query spectra with replacement; only useful for computing confidence intervals of parameter estimates (optional, options=[true,false], default=false).
  --random_seed                   Random seed to be used in all computations with a stochastic component (optional, default=1).
  --params_to_optimize            String denoting the parameters to optimize (optional, options='all','LET_thresh','noise_thresh','wf_int','wf_mz', default='all').
  --spectrum_preprocessing_order  String denoting the order of spectrum preprocessing transformations; format must be a string with 0-3 characters of either L (low-entropy trannsformation) and/or W (weight-factor-transformation) and/or N (noise removal) (optional, default: 'NWL').
  --threads                       Number of threads to use (optional, default=1).
  --max_steps                     Maximum number of iterations allowed in differential evolution optimization; only applicable for DE optimization_method (optional, default=5).
  --pop_size                      Population size in differential evolution optimization; only applicable for DE optimization_method (optional, default=50).
  --n_grid_points                 Number of grid points to use for each parameter; only applicable for grid-based optimization (optional, default=2).
  --LB_LET_thresh                 Float denoting the lower bound of the low-entropy threshold parameter (optional, default=0.0).
  --UB_LET_thresh                 Float denoting the upper bound of the low-entropy threshold parameter (optional, default=5.0).
  --LB_noise_thresh               Float denoting the lower bound of the noise removal threshold parameter (optional, default=0.0).
  --UB_noise_thresh               Float denoting the upper bound of the noise removal threshold parameter (optional, default=1.0).
  --LB_wf_int                     Float denoting the lower bound of the intensity weight factor parameter (optional, default=0.0).
  --UB_wf_int                     Float denoting the upper bound of the intensity weight factor parameter (optional, default=5.0).
  --LB_wf_mz                      Float denoting the lower bound of the mass/charge weight factor parameter (optional, default=0.0).
  --UB_wf_mz                      Float denoting the upper bound of the mass/charge weight factor parameter (optional, default=5.0).
  --LET_thresh                    Float denoting the low-entropy threshold parameter; not applicable for optimization_method = grid (optional, default=0.0).
  --noise_thresh                  Float denoting the noise removal threshold parameter; not applicable for optimization_method = grid (optional, default=0.1).
  --wf_int                        Float denoting the intensity weight factor parameter; not applicable for optimization_method = grid (optional, default=1.0).
  --wf_mz                         Float denoting the mass/charge weight factor parameter; not applicable for optimization_method = grid(optional, default=0.0).
  --help                          Show this help message.
""")
end


function parse_args()
    args = Dict{String, String}()
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--help" || ARGS[i] == "-h"
            print_help()
            exit(0)
        elseif startswith(ARGS[i], "--")
            if i == length(ARGS)
                error("Missing value for argument: ", ARGS[i])
            end
            args[ARGS[i]] = ARGS[i+1]
            i += 2
        else
            error("Unknown argument: ", ARGS[i])
        end
    end
    for req in ["--query_data", "--reference_data", "--output"]
        haskey(args, req) || error("Missing required argument: $req. Use --help for usage.")
    end
    get!(args, "--metric", "accuracy")
    get!(args, "--crossvalidation", "true")
    get!(args, "--n_folds", "5")
    get!(args, "--bootstrap_query", "false")
    get!(args, "--random_seed", "1")
    get!(args, "--params_to_optimize", "all")
    get!(args, "--spectrum_preprocessing_order", "NWL")
    get!(args, "--threads", "1")
    get!(args, "--n_grid_points", "2")
    get!(args, "--max_steps", "5")
    get!(args, "--pop_size", "50")
    get!(args, "--LB_LET_thresh", "0.0")
    get!(args, "--UB_LET_thresh", "5.0")
    get!(args, "--LB_noise_thresh", "0.0")
    get!(args, "--UB_noise_thresh", "5.0")
    get!(args, "--LB_wf_int", "0.0")
    get!(args, "--UB_wf_int", "5.0")
    get!(args, "--LB_wf_mz", "0.0")
    get!(args, "--UB_wf_mz", "1.3")
    get!(args, "--LET_thresh", "0.0")
    get!(args, "--noise_thresh", "0.0")
    get!(args, "--wf_int", "1.0")
    get!(args, "--wf_mz", "0.0")

    if !(args["--metric"] in ["accuracy","MRR"])
        println("Warning: metric must be either 'accuracy' or 'MRR'.")
    end

    if !(args["--spectrum_preprocessing_order"] in ["","L","N","W","LN","NL","LW","WL","NW","WN","LNW","LWN","NLW","NWL","LNW","LWN"])
        println("Warning: spectrum_preprocessing_order must be either '','L','N','W','LN','NL','LW','WL','NW','WN','LNW','LWN','NLW','NWL','LNW','LWN'.")
    end

    if !(args["--params_to_optimize"] in ["all","LET_thresh","noise_thresh","wf_int","wf_mz"])
        println("Error: invalid params_to_optimize parameter. Run <julia OptiMS.jl --help> for usage instructions")
    end

    if !(args["--optimization_method"] in ["DE","grid","none"])
        println("Warning: optimization_method must be either 'DE', 'grid', or 'none'.")
    end

    if !(args["--crossvalidation"] in ["true","false"])
        println("Warning: crossvalidation must be either 'true' or 'false'.")
    end

    if !(args["--bootstrap_query"] in ["true","false"])
        println("Warning: bootstrap_query must be either 'true' or 'false'.")
    end

    return args
end


function make_folds(n::Int, K::Int; rng::AbstractRNG)
    idx = collect(1:n)
    shuffle!(rng, idx)
    folds = [Int[] for _ in 1:K]
    for (t, i) in enumerate(idx)
        push!(folds[mod1(t, K)], i)
    end
    return folds
end



function wf_transformation(X::AbstractMatrix, wf_mz::Real, wf_int::Real; mzs::AbstractVector)
    wrow = reshape(mzs .^ wf_mz, 1, :)
    return (X .^ wf_int) .* wrow
end


function LE_transformation(X::AbstractMatrix, LET_thresh)
    T = promote_type(eltype(X), typeof(LET_thresh))
    Xf = Array{T}(X)
    rs = sum(Xf, dims=2)
    P  = similar(Xf); fill!(P, zero(T))
    nz = vec(rs) .> zero(T)
    P[nz, :] .= Xf[nz, :] ./ rs[nz, :]
    Tmat = zeros(T, size(P))
    idx = P .> zero(T)
    @inbounds Tmat[idx] .= P[idx] .* log.(P[idx])
    Sv = vec(-sum(Tmat, dims=2))
    lt_mask = (Sv .> zero(T)) .& (Sv .< LET_thresh)
    if any(lt_mask)
        w = (one(T) .+ Sv) ./ (one(T) .+ LET_thresh)
        out = copy(P)
        @threads for i in eachindex(Sv)
            if lt_mask[i]
                @inbounds out[i, :] .= P[i, :].^w[i]
            end
        end
        return out
    else
        return P
    end
end


function remove_noise(X::AbstractMatrix, noise_thresh::Real)
    Xf = Array{Float64}(X)
    rowmax = maximum(Xf, dims=2)
    cutoff = noise_thresh .* rowmax
    out = copy(Xf)
    @threads for i in 1:size(out, 1)
        @inbounds @views begin
            r = out[i, :]
            c = cutoff[i]
            for j in eachindex(r)
                if r[j] < c
                    r[j] = 0.0
                end
            end
        end
    end
    return out
end


function row_l2_normalize!(X::AbstractMatrix)
    norms = sqrt.(sum(abs2, X; dims=2))
    nz = vec(norms) .> 0.0
    X[nz, :] .= X[nz, :] ./ norms[nz, :]
    return X
end


function get_acc(Q_mat::AbstractMatrix, R_mat::AbstractMatrix, q_ids_all::AbstractVector{<:AbstractString}, r_ids_all::AbstractVector{<:AbstractString})
    Qn = copy(Q_mat); Rn = copy(R_mat)
    row_l2_normalize!(Qn); row_l2_normalize!(Rn)
    S = Qn * Rn'
    preds = Vector{String}(undef, size(S, 1))
    hits  = falses(size(S, 1))
    @threads for i in 1:size(S, 1)
        srow = @view S[i, :]
        _, j = findmax(srow)
        preds[i] = String(r_ids_all[j])
        hits[i] = (String(q_ids_all[i]) == preds[i])
    end
    acc = count(hits) / length(hits)
    return acc
end


function get_MRR(Q_mat::AbstractMatrix, R_mat::AbstractMatrix, q_ids_all::AbstractVector{<:AbstractString}, r_ids_all::AbstractVector{<:AbstractString})
    Qn = copy(Q_mat); Rn = copy(R_mat)
    row_l2_normalize!(Qn); row_l2_normalize!(Rn)
    S = Qn * Rn'
    ref_index = Dict{String, Int}(String(r_ids_all[j]) => j for j in eachindex(r_ids_all))
    rr = Vector{Float64}(undef, size(S, 1))
    @threads for i in 1:size(S, 1)
        srow = @view S[i, :]
        _, jmax = findmax(srow)
        qid = String(q_ids_all[i])
        jtrue = ref_index[qid]
        ord = sortperm(srow; rev=true)
        rnk = findfirst(==(jtrue), ord)
        rr[i] = 1.0 / rnk
    end
    MRR = sum(rr) / length(rr)
    return MRR
end


function apply_pipeline(Q::AbstractMatrix, R::AbstractMatrix; order::AbstractString, LET_thresh, noise_thresh, wf_int, wf_mz, mzs::AbstractVector)
    T = promote_type(eltype(Q), eltype(R), typeof(LET_thresh), typeof(noise_thresh), typeof(wf_int), typeof(wf_mz))
    Q_mat = Array{T}(Q)
    R_mat = Array{T}(R)
    for c in order
        if c == 'W'
            Q_mat = wf_transformation(Q_mat, wf_mz, wf_int; mzs=mzs)
            R_mat = wf_transformation(R_mat, wf_mz, wf_int; mzs=mzs)
        elseif c == 'L'
            Q_mat = LE_transformation(Q_mat, LET_thresh)
            R_mat = LE_transformation(R_mat, LET_thresh)
        elseif c == 'N'
            Q_mat = remove_noise(Q_mat, noise_thresh)
            R_mat = remove_noise(R_mat, noise_thresh)
        else
            error("Unknown transform code '$c'. Use only 'W','L','N'.")
        end
    end
    return Q_mat, R_mat
end


function get_scores(Q0, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
    Qn = copy(Q0); Rn = copy(R0)
    Qp, Rp = apply_pipeline(Qn, Rn; order=order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
    row_l2_normalize!(Qp); row_l2_normalize!(Rp)
    S = Qp * Rp'
    return S
end


function objective_acc_cv(x::Vector)
    LET_thresh, noise_thresh, wf_int, wf_mz = x
    min_acc = 99999
    for k in 1:K
        val_idx = folds[k]
        Q0_val = @view Q0[val_idx, :]
        q_ids_val = q_ids_all[val_idx]
        Qp_val, Rp = apply_pipeline(Q0_val, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
        acc_k = get_acc(Qp_val, Rp, q_ids_val, r_ids_all)
        min_acc = min(min_acc, acc_k)
        if min_acc == 0.0
            break
        end
    end
    return 1.0 - min_acc
end


function objective_MRR_cv(x)
    LET_thresh, noise_thresh, wf_int, wf_mz = x
    min_MRR = 99999
    for k in 1:K
        val_idx = folds[k]
        Q0_val = @view Q0[val_idx, :]
        q_ids_val = q_ids_all[val_idx]
        Qp_val, Rp = apply_pipeline(Q0_val, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
        MRR_k = get_MRR(Qp_val, Rp, q_ids_val, r_ids_all)
        min_MRR = min(min_MRR, MRR_k)
        if min_MRR == 0.0
            break
        end
    end
    return 1.0 - min_MRR
end


function objective_acc_no_cv(x)
    LET_thresh, noise_thresh, wf_int, wf_mz = x
    Qp, Rp = apply_pipeline(Q0, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs = mzs)
    acc = get_acc(Qp, Rp, q_ids_all, r_ids_all)
    return 1.0 - acc
end


function objective_MRR_no_cv(x)
    LET_thresh, noise_thresh, wf_int, wf_mz = x
    Qp, Rp = apply_pipeline(Q0, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs = mzs)
    MRR = get_MRR(Qp, Rp, q_ids_all, r_ids_all)
    return 1.0 - MRR
end


