
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
    --platform <string> \
    --optimization_method <string> \
    --metric <string> \
    --crossvalidation <boolean> \
    --filter_reference_candidates <boolean> \
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
    --max_steps <int> \
    --pop_size <int> \
    --n_grid_points <int>


Arguments:
  --query_data                    Path to input TXT file of query dataset (required).
  --reference_data                Path to input TXT file of reference dataset (required).
  --output                        Path to output TXT file (required).
  --platform                      Chromatography platform (required, options=[HRMS,NRMS]).
  --optimization_method           Optimization approach (optional, options=[DE,grid,none], default=DE).
  --metric                        Quantity to maximize in the objective function (optional, options=[accuracy,MRR], default=accuracy).
  --crossvalidation               Boolean indicating whether or not to perform 5-fold cross-validation inside the objective function (optional, options=[true,false], default=true).
  --filter_reference_candidates   Boolean indicating whether or not to filter on precursor ion mass/charge; only applicable to HRMS data (optional, options=[true,false], default=true).
  --n_folds                       Number folds to use for cross-validation. Only applicable for crossvalidation is \'true\' and optimization_method is not \'none\' (optional, default=5).
  --bootstrap_query               Boolean indicating whether or not to construct query dataset of same size from resampling query spectra with replacement; only useful for computing confidence intervals of parameter estimates (optional, options=[true,false], default=false).
  --random_seed                   Random seed to be used in all computations with a stochastic component (optional, default=1).
  --params_to_optimize            String denoting the parameters to optimize (optional, options='all','LET_thresh','noise_thresh','wf_int','wf_mz', default='all').
  --spectrum_preprocessing_order  String denoting the order of spectrum preprocessing transformations; format must be a string with 2-6 characters chosen from C, F, M, N, L, W representing centroiding, filtering based on mass/charge and intensity values, matching, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively; note that if an argument is passed for HRMS, then \'M\' must be contained in the argument since matching is a required preprocessing step in spectral library matching of HRMS data; furthermore, \'C\' must be performed before matching since centroiding can change the number of ion fragments in a given spectrum; note that C and M are not applicable to NRMS data (optional, default: NCMWL for HRMS and NWL for NRMS').
  --threads                       Number of threads to use (optional, default=1).
  --max_steps                     Maximum number of iterations allowed in differential evolution optimization; only applicable for DE optimization_method (optional, default=5).
  --pop_size                      Population size in differential evolution optimization; only applicable for DE optimization_method (optional, default=50).
  --n_grid_points                 Number of grid points to use for each parameter; only applicable for grid-based optimization (optional, default=2).
  --LB_ws_matching                Float denoting the lower bound of the matching window-size parameter (optional, default=0.01).
  --UB_ws_matching                Float denoting the upper bound of the matching window-size parameter (optional, default=0.9).
  --LB_ws_centroiding             Float denoting the lower bound of the centroiding window-size parameter (optional, default=0.01).
  --UB_ws_centroiding             Float denoting the upper bound of the centroiding window-size parameter (optional, default=0.9).
  --LB_LET_thresh                 Float denoting the lower bound of the low-entropy threshold parameter (optional, default=0.0).
  --UB_LET_thresh                 Float denoting the upper bound of the low-entropy threshold parameter (optional, default=5.0).
  --LB_noise_thresh               Float denoting the lower bound of the noise removal threshold parameter (optional, default=0.0).
  --UB_noise_thresh               Float denoting the upper bound of the noise removal threshold parameter (optional, default=1.0).
  --LB_wf_int                     Float denoting the lower bound of the intensity weight factor parameter (optional, default=0.0).
  --UB_wf_int                     Float denoting the upper bound of the intensity weight factor parameter (optional, default=5.0).
  --LB_wf_mz                      Float denoting the lower bound of the mass/charge weight factor parameter (optional, default=0.0).
  --UB_wf_mz                      Float denoting the upper bound of the mass/charge weight factor parameter (optional, default=5.0).
  --ws_matching                   Float denoting the matching window-size parameter; only applicable for optimization_method = none (optional, default=0.05).
  --ws_centroiding                Float denoting the centroiding window-size parameter; only applicable for optimization_method = none (optional, default=0.05).
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
    for req in ["--query_data", "--reference_data", "--output", "--platform"]
        haskey(args, req) || error("Missing required argument: $req. Use --help for usage.")
    end
    if args["--platform"] == "HRMS"
        get!(args, "--spectrum_preprocessing_order", "NCMWL")
    elseif args["--platform"] == "NRMS"
        get!(args, "--spectrum_preprocessing_order", "NWL")
    else
        println("Error: platform must be either 'HRMS' or 'NRMS'.")
        exit(1)
    end
    get!(args, "--metric", "accuracy")
    get!(args, "--crossvalidation", "true")
    get!(args, "--filter_reference_candidates", "true")
    get!(args, "--n_folds", "5")
    get!(args, "--bootstrap_query", "false")
    get!(args, "--random_seed", "1")
    get!(args, "--params_to_optimize", "all")
    get!(args, "--spectrum_preprocessing_order", "NWL")
    get!(args, "--threads", "1")
    get!(args, "--n_grid_points", "2")
    get!(args, "--max_steps", "5")
    get!(args, "--pop_size", "50")
    get!(args, "--LB_ws_matching", "0.01")
    get!(args, "--UB_ws_matching", "0.09")
    get!(args, "--LB_ws_centroiding", "0.01")
    get!(args, "--UB_ws_centroiding", "0.09")
    get!(args, "--LB_LET_thresh", "0.0")
    get!(args, "--UB_LET_thresh", "5.0")
    get!(args, "--LB_noise_thresh", "0.0")
    get!(args, "--UB_noise_thresh", "5.0")
    get!(args, "--LB_wf_int", "0.0")
    get!(args, "--UB_wf_int", "5.0")
    get!(args, "--LB_wf_mz", "0.0")
    get!(args, "--UB_wf_mz", "1.3")
    get!(args, "--ws_matching", "0.05")
    get!(args, "--ws_centroiding", "0.05")
    get!(args, "--LET_thresh", "0.0")
    get!(args, "--noise_thresh", "0.0")
    get!(args, "--wf_int", "1.0")
    get!(args, "--wf_mz", "0.0")

    if !(args["--metric"] in ["accuracy","MRR"])
        println("Warning: metric must be either 'accuracy' or 'MRR'.")
    end

    #if !(args["--spectrum_preprocessing_order"] in ["","L","N","W","LN","NL","LW","WL","NW","WN","LNW","LWN","NLW","NWL","LNW","LWN"])
    #    println("Warning: spectrum_preprocessing_order must be either '','L','N','W','LN','NL','LW','WL','NW','WN','LNW','LWN','NLW','NWL','LNW','LWN'.")
    #end

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


function make_folds_NRMS(n::Int, K::Int; rng::AbstractRNG)
    idx = collect(1:n)
    shuffle!(rng, idx)
    folds = [Int[] for _ in 1:K]
    for (t, i) in enumerate(idx)
        push!(folds[mod1(t, K)], i)
    end
    return folds
end


function wf_transformation_NRMS(X::AbstractMatrix, wf_mz::Real, wf_int::Real; mzs::AbstractVector)
    wrow = reshape(mzs .^ wf_mz, 1, :)
    return (X .^ wf_int) .* wrow
end


function LE_transformation_NRMS(X::AbstractMatrix, LET_thresh)
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


function remove_noise_NRMS(X::AbstractMatrix, noise_thresh::Real)
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


function row_l2_normalize_NRMS!(X::AbstractMatrix)
    norms = sqrt.(sum(abs2, X; dims=2))
    nz = vec(norms) .> 0.0
    X[nz, :] .= X[nz, :] ./ norms[nz, :]
    return X
end


function get_acc_NRMS(Q_mat::AbstractMatrix, R_mat::AbstractMatrix, q_ids_all::AbstractVector{<:AbstractString}, r_ids_all::AbstractVector{<:AbstractString})
    Qn = copy(Q_mat); Rn = copy(R_mat)
    row_l2_normalize_NRMS!(Qn); row_l2_normalize_NRMS!(Rn)
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


function get_MRR_NRMS(Q_mat::AbstractMatrix, R_mat::AbstractMatrix, q_ids_all::AbstractVector{<:AbstractString}, r_ids_all::AbstractVector{<:AbstractString})
    Qn = copy(Q_mat); Rn = copy(R_mat)
    row_l2_normalize_NRMS!(Qn); row_l2_normalize_NRMS!(Rn)
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


function apply_pipeline_NRMS(Q::AbstractMatrix, R::AbstractMatrix; order::AbstractString, LET_thresh, noise_thresh, wf_int, wf_mz, mzs::AbstractVector)
    T = promote_type(eltype(Q), eltype(R), typeof(LET_thresh), typeof(noise_thresh), typeof(wf_int), typeof(wf_mz))
    Q_mat = Array{T}(Q)
    R_mat = Array{T}(R)
    for c in order
        if c == 'W'
            Q_mat = wf_transformation_NRMS(Q_mat, wf_mz, wf_int; mzs=mzs)
            R_mat = wf_transformation_NRMS(R_mat, wf_mz, wf_int; mzs=mzs)
        elseif c == 'L'
            Q_mat = LE_transformation_NRMS(Q_mat, LET_thresh)
            R_mat = LE_transformation_NRMS(R_mat, LET_thresh)
        elseif c == 'N'
            Q_mat = remove_noise_NRMS(Q_mat, noise_thresh)
            R_mat = remove_noise_NRMS(R_mat, noise_thresh)
        else
            error("Unknown transform code '$c'. Use only 'W','L','N'.")
        end
    end
    return Q_mat, R_mat
end


function get_scores_NRMS(Q0, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
    Qn = copy(Q0); Rn = copy(R0)
    Qp, Rp = apply_pipeline_NRMS(Qn, Rn; order=order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
    row_l2_normalize_NRMS!(Qp); row_l2_normalize_NRMS!(Rp)
    S = Qp * Rp'
    return S
end


function objective_acc_cv_NRMS(x::Vector)
    LET_thresh, noise_thresh, wf_int, wf_mz = x
    min_acc = 99999
    for k in 1:K
        val_idx = folds[k]
        Q0_val = @view Q0[val_idx, :]
        q_ids_val = q_ids_all[val_idx]
        Qp_val, Rp = apply_pipeline_NRMS(Q0_val, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
        acc_k = get_acc_NRMS(Qp_val, Rp, q_ids_val, r_ids_all)
        min_acc = min(min_acc, acc_k)
        if min_acc == 0.0
            break
        end
    end
    return 1.0 - min_acc
end


function objective_MRR_cv_NRMS(x)
    LET_thresh, noise_thresh, wf_int, wf_mz = x
    min_MRR = 99999
    for k in 1:K
        val_idx = folds[k]
        Q0_val = @view Q0[val_idx, :]
        q_ids_val = q_ids_all[val_idx]
        Qp_val, Rp = apply_pipeline_NRMS(Q0_val, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs=mzs)
        MRR_k = get_MRR_NRMS(Qp_val, Rp, q_ids_val, r_ids_all)
        min_MRR = min(min_MRR, MRR_k)
        if min_MRR == 0.0
            break
        end
    end
    return 1.0 - min_MRR
end


function objective_acc_no_cv_NRMS(x)
    LET_thresh, noise_thresh, wf_int, wf_mz = x
    Qp, Rp = apply_pipeline_NRMS(Q0, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs = mzs)
    acc = get_acc_NRMS(Qp, Rp, q_ids_all, r_ids_all)
    return 1.0 - acc
end


function objective_MRR_no_cv_NRMS(x)
    LET_thresh, noise_thresh, wf_int, wf_mz = x
    Qp, Rp = apply_pipeline_NRMS(Q0, R0; order=spectrum_preprocessing_order, LET_thresh=LET_thresh, noise_thresh=noise_thresh, wf_int=wf_int, wf_mz=wf_mz, mzs = mzs)
    MRR = get_MRR_NRMS(Qp, Rp, q_ids_all, r_ids_all)
    return 1.0 - MRR
end







function make_folds_HRMS(ids::Vector{String}, K::Int; rng::AbstractRNG)
    ids_shuf = copy(ids)
    shuffle!(rng, ids_shuf)
    folds = [String[] for _ in 1:K]
    for (t, id) in enumerate(ids_shuf)
        push!(folds[mod1(t, K)], id)
    end
    return folds
end

function normalize_intensity_HRMS(ints::AbstractVector{<:Real})
    s = sum(ints)
    s > 0 ? Float64.(ints) ./ s : zeros(Float64, length(ints))
end


function wf_transform_HRMS!(spec::AbstractMatrix{<:Real}, wf_mz::Real, wf_int::Real)
    @views spec[:, 2] .= (Float64.(spec[:, 1]).^float(wf_mz)) .* (Float64.(spec[:, 2]).^float(wf_int))
    return spec
end


function LE_transform_HRMS!(spec::AbstractMatrix{<:Real}, thresh::Real)
    n = size(spec,1)
    n == 0 && return spec
    @views p = copy(spec[:, 2])
    s = sum(p)
    if s <= 0
        @views spec[:, 2] .= 0.0
        return spec
    end
    p ./= s
    S = 0.0
    @inbounds for v in p
        if v > 0
            S -= v * log(v)
        end
    end
    if 0.0 < S < float(thresh)
        w = (1 + S) / (1 + float(thresh))
        @views spec[:, 2] .= p .^ w
    end
    return spec
end


function remove_noise_HRMS!(spec::AbstractMatrix{<:Real}, nr::Union{Nothing,Real})
    (nr === nothing || size(spec,1) <= 1) && return spec
    @views m = maximum(spec[:, 2])
    cutoff = m * float(nr)
    @inbounds @views for i in axes(spec,1)
        if spec[i, 2] < cutoff
            spec[i, 2] = 0.0
        end
    end
    return spec
end


@views function centroid_spectrum_HRMS(spec::Matrix{Float64}, window_size::Float64)
    n = size(spec, 1)
    n == 0 && return spec
    spec = copy(spec)
    spec = spec[sortperm(spec[:, 1]), :]   # sort by m/z
    need_centroid = false
    if n > 1
        mindelta = minimum(spec[2:end, 1] .- spec[1:end-1, 1])
        need_centroid = (mindelta <= window_size)
    end
    !need_centroid && return spec
    intensity_order = sortperm(spec[:, 2]; rev=true)
    spec_new = Vector{NTuple{2,Float64}}()
    @inbounds for i in intensity_order
        if spec[i, 2] > 0
            mz0 = spec[i, 1]
            i_left = i - 1
            while i_left >= 1 && (mz0 - spec[i_left, 1]) <= window_size
                i_left -= 1
            end
            i_left += 1
            i_right = i + 1
            while i_right <= n && (spec[i_right, 1] - mz0) <= window_size
                i_right += 1
            end
            i_right -= 1
            slice = i_left:i_right
            vI = spec[slice, 2]
            vM = spec[slice, 1]
            I_sum   = sum(vI)
            mzI_sum = sum(vM .* vI)
            push!(spec_new, (mzI_sum / I_sum, I_sum))
            vI .= 0.0
        end
    end
    if isempty(spec_new)
        return zeros(1, 2)
    end
    out = Matrix{Float64}(undef, length(spec_new), 2)
    @inbounds for k in 1:length(spec_new)
        out[k, 1] = spec_new[k][1]
        out[k, 2] = spec_new[k][2]
    end
    out = out[sortperm(out[:, 1]), :]
    return size(out, 1) > 1 ? out : zeros(1, 2)
end



function apply_transformations_HRMS(spec::Matrix{Float64}, order::AbstractString; window_size_centroiding::Float64, noise_threshold::Union{Nothing,Float64}, LET_threshold::Float64, wf_mz::Float64, wf_int::Float64)
    spec = Matrix(spec)
    if size(spec,1) > 1
        spec = spec[sortperm(spec[:, 1]), :]
    end
    for ch in uppercase(order)
        if ch == 'N'
            remove_noise_HRMS!(spec, noise_threshold)
        elseif ch == 'C'
            spec = centroid_spectrum_HRMS(spec, window_size_centroiding)
        elseif ch == 'W'
            wf_transform_HRMS!(spec, wf_mz, wf_int)
        elseif ch == 'L'
            LE_transform_HRMS!(spec, LET_threshold)
        elseif ch == 'M'
            nothing
        else
            nothing
        end
    end
    return spec
end


function match_peaks_in_spectra_HRMS(spec_a::AbstractMatrix{<:Real}, spec_b::AbstractMatrix{<:Real}, window_size::Real)
    a = 1; b = 1
    na = size(spec_a, 1); nb = size(spec_b, 1)
    merged = Vector{NTuple{3,Float64}}()
    peak_b_int = 0.0
    while a <= na && b <= nb
        mass_delta = float(spec_a[a, 1]) - float(spec_b[b, 1])
        if mass_delta < -float(window_size)
            push!(merged, (float(spec_a[a, 1]), float(spec_a[a, 2]), peak_b_int))
            peak_b_int = 0.0; a += 1
        elseif mass_delta >  float(window_size)
            push!(merged, (float(spec_b[b, 1]), 0.0,                 float(spec_b[b, 2])))
            b += 1
        else
            peak_b_int += float(spec_b[b, 2])
            b += 1
        end
    end
    if a <= na && peak_b_int > 0.0
        push!(merged, (float(spec_a[a, 1]), float(spec_a[a, 2]), peak_b_int))
        a += 1
    end
    if b <= nb
        for i in b:nb
            push!(merged, (float(spec_b[i, 1]), 0.0, float(spec_b[i, 2])))
        end
    end
    if a <= na
        for i in a:na
            push!(merged, (float(spec_a[i, 1]), float(spec_a[i, 2]), 0.0))
        end
    end
    if isempty(merged)
        return [0.0 0.0 0.0]
    else
        M = Matrix{Float64}(undef, length(merged), 3)
        @inbounds for i in 1:length(merged)
            t = merged[i]; M[i,1]=t[1]; M[i,2]=t[2]; M[i,3]=t[3]
        end
        return M
    end
end


@inline function cosine_similarity_from_merged_HRMS(merged::Matrix{Float64})
    ia = view(merged, :, 2); ib = view(merged, :, 3)
    na = norm(ia); nb = norm(ib)
    (na == 0.0 || nb == 0.0) && return 0.0
    return dot(ia, ib) / (na * nb)
end


#=
function build_spectra_by_id_prec_HRMS(df::DataFrame; idcol::Symbol=:id, prec::Symbol=:precursor_ion_mz, mzcol::Symbol=:mz_ratio, intcol::Symbol=:intensity, keep_ids::AbstractVector{<:AbstractString}=String[])
    if !isempty(keep_ids)
        df = semijoin(df, DataFrame(id=keep_ids), on=idcol)
    end
    needed = (idcol, prec, mzcol, intcol)
    df = df[reduce(.&, (.!ismissing.(df[!, c]) for c in needed)), :]
    keys  = NamedTuple{(:id,:prec),Tuple{String,Float64}}[]
    precs = Float64[]
    specs = Matrix{Float64}[]
    for sub in groupby(df, [idcol, prec])
        id   = String(sub[1, idcol])
        p_mz = Float64(sub[1, prec])
        sort!(sub, mzcol)
        m = nrow(sub)
        mat = Matrix{Float64}(undef, m, 2)
        mat[:,1] = Float64.(sub[!, mzcol])
        mat[:,2] = Float64.(sub[!, intcol])
        push!(keys, (id=id, prec=p_mz))
        push!(precs, p_mz)
        push!(specs, mat)
    end
    return keys, precs, specs
end
=#

function build_spectra_by_id_prec_HRMS(df::DataFrame; idcol::Symbol = :id, prec::Symbol = :precursor_ion_mz, mzcol::Symbol = :mz_ratio, intcol::Symbol = :intensity, keep_ids::AbstractVector{<:AbstractString} = String[], bootstrap_query::Bool = false, bootstrap_spectrum_id_col::Symbol = :bootstrap_spectrum_id)
    if !isempty(keep_ids)
        df = semijoin(df, DataFrame(id = keep_ids), on = idcol)
    end
    needed = bootstrap_query ?
        (idcol, bootstrap_spectrum_id_col, prec, mzcol, intcol) :
        (idcol, prec, mzcol, intcol)
    df = df[reduce(.&, (.!ismissing.(df[!, c]) for c in needed)), :]
    keys  = NamedTuple{(:id, :prec), Tuple{String, Float64}}[]
    precs = Float64[]
    specs = Matrix{Float64}[]
    group_cols = bootstrap_query ? [idcol, bootstrap_spectrum_id_col, prec] : [idcol, prec]
    for sub in groupby(df, group_cols)
        id   = String(sub[1, idcol])
        p_mz = Float64(sub[1, prec])
        sort!(sub, mzcol)
        m = nrow(sub)
        mat = Matrix{Float64}(undef, m, 2)
        mat[:, 1] = Float64.(sub[!, mzcol])
        mat[:, 2] = Float64.(sub[!, intcol])
        push!(keys, (id = id, prec = p_mz))
        push!(precs, p_mz)
        push!(specs, mat)
    end
    return keys, precs, specs
end


@inline function cosine_similarity_HRMS(q_post_t::AbstractMatrix{<:Real}, r_post_t::AbstractMatrix{<:Real})
    @assert size(q_post_t, 2) == 2 "q_post_t must be N×2, got $(size(q_post_t))"
    @assert size(r_post_t, 2) == 2 "r_post_t must be N×2, got $(size(r_post_t))"
    n = size(q_post_t, 1)
    @assert size(r_post_t, 1) == n "row mismatch: $(size(q_post_t)) vs $(size(r_post_t))"
    num = 0.0
    na2 = 0.0
    nb2 = 0.0
    @inbounds for i in 1:n
        a = float(q_post_t[i, 2])
        b = float(r_post_t[i, 2])
        num += a * b
        na2 += a * a
        nb2 += b * b
    end
    denom = sqrt(na2) * sqrt(nb2)
    return denom == 0.0 ? 0.0 : (num / denom)
end


function objective_no_cv_acc_filtering_HRMS(x)
    LET_threshold = x[1]
    noise_threshold = x[2]
    wf_int = x[3]
    wf_mz = x[4]
    window_size_matching = x[5]
    window_size_centroiding  = x[6]
    q_keys, q_precs, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :int, keep_ids = String[])
    r_keys, r_precs, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :int, keep_ids = String[])
    neval    = zeros(Int, Threads.maxthreadid())
    ncorrect = zeros(Int, Threads.maxthreadid())
    @threads for i in eachindex(q_keys)
        tid = threadid()
        qprec = q_precs[i]
        mask = @. abs(r_precs - qprec) <=  0.5
        if !any(mask)
            continue
        end
        neval[tid] += 1
        qspec_t = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
        best_cos = -Inf
        best_j = 0
        @inbounds for j in eachindex(r_specs)
            mask[j] || continue
            rspec_t = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
            merged = match_peaks_in_spectra_HRMS(qspec_t, rspec_t, window_size_matching)
            cosval = cosine_similarity_from_merged_HRMS(merged)
            if cosval > best_cos
                best_cos = cosval
                best_j = j
            end
        end
        if best_j != 0 && r_keys[best_j].id == q_keys[i].id
            ncorrect[tid] += 1
        end
    end
    n_eval = sum(neval)
    n_correct = sum(ncorrect)
    return 1.0 - (n_correct / n_eval)
end


function objective_no_cv_acc_no_filtering_HRMS(x)
    LET_threshold = x[1]
    noise_threshold = x[2]
    wf_int = x[3]
    wf_mz = x[4]
    window_size_matching = x[5]
    window_size_centroiding  = x[6]
    q_keys, q_precs, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    r_keys, r_precs, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    neval    = zeros(Int, Threads.maxthreadid())
    ncorrect = zeros(Int, Threads.maxthreadid())
    @threads for i in eachindex(q_keys)
        tid = threadid()
        isempty(r_specs) && continue
        neval[tid] += 1
        qspec_t = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
        best_cos = -Inf
        best_j   = 0
        @inbounds for j in eachindex(r_specs)
            rspec_t = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
            merged = match_peaks_in_spectra_HRMS(qspec_t, rspec_t, window_size_matching)
            cosval = cosine_similarity_from_merged_HRMS(merged)
            if cosval > best_cos
                best_cos = cosval
                best_j   = j
            end
        end
        if best_j != 0 && r_keys[best_j].id == q_keys[i].id
            ncorrect[tid] += 1
        end
    end
    n_eval = sum(neval)
    n_correct = sum(ncorrect)
    return 1.0 - (n_correct / n_eval)
end


function objective_no_cv_MRR_filtering_HRMS(x)
    LET_threshold = x[1]
    noise_threshold = x[2]
    wf_int = x[3]
    wf_mz = x[4]
    window_size_matching = x[5]
    window_size_centroiding  = x[6]
    q_keys, q_precs, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    r_keys, r_precs, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    neval    = zeros(Int, Threads.maxthreadid())
    rr_sum = zeros(Float64, Threads.maxthreadid())
    @threads for i in eachindex(q_keys)
        tid = threadid()
        qprec = q_precs[i]
        mask = @. abs(r_precs - qprec) <= 0.5
        if !any(mask)
            continue
        end
        neval[tid] += 1
        qspec_t = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
        scores = Float64[]
        ids    = String[]
        @inbounds for j in eachindex(r_specs)
            mask[j] || continue
            rspec_t = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
            merged = match_peaks_in_spectra_HRMS(qspec_t, rspec_t, window_size_matching)
            cosval = cosine_similarity_from_merged_HRMS(merged)
            push!(scores, cosval)
            push!(ids, r_keys[j].id)
        end
        perm = sortperm(scores, rev = true)
        qid = q_keys[i].id
        for (rank, k) in enumerate(perm)
            if ids[k] == qid
                rr_sum[tid] += 1.0 / rank
                break
            end
        end
    end
    n_eval = sum(neval)
    return 1.0 - (sum(rr_sum) / n_eval)
end


function objective_no_cv_MRR_no_filtering_HRMS(x)
    LET_threshold = x[1]
    noise_threshold = x[2]
    wf_int = x[3]
    wf_mz = x[4]
    window_size_matching = x[5]
    window_size_centroiding  = x[6]
    q_keys, q_precs, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    r_keys, r_precs, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    neval    = zeros(Int, Threads.maxthreadid())
    rr_sum = zeros(Float64, Threads.maxthreadid())
    @threads for i in eachindex(q_keys)
        tid = threadid()
        isempty(r_specs) && continue
        neval[tid] += 1
        qspec_t = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
        scores = Vector{Float64}(undef, length(r_specs))
        @inbounds for j in eachindex(r_specs)
            rspec_t = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
            merged = match_peaks_in_spectra_HRMS(qspec_t, rspec_t, window_size_matching)
            scores[j] = cosine_similarity_from_merged_HRMS(merged)
        end
        perm = sortperm(scores, rev = true)
        qid = q_keys[i].id
        for (rank, j) in enumerate(perm)
            if r_keys[j].id == qid
                rr_sum[tid] += 1.0 / rank
                break
            end
        end
    end
    n_eval = sum(neval)
    return 1.0 - (sum(rr_sum) / n_eval)
end


function objective_cv_acc_filtering_HRMS(x)
    LET_threshold = x[1]
    noise_threshold = x[2]
    wf_int = x[3]
    wf_mz = x[4]
    window_size_matching = x[5]
    window_size_centroiding  = x[6]
    q_keys, q_precs, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    r_keys, r_precs, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    neval = zeros(Int, Threads.maxthreadid())
    ncorrect = zeros(Int, Threads.maxthreadid())
    @threads for i in eachindex(q_keys)
        tid = threadid()
        qprec = q_precs[i]
        mask = @. abs(r_precs - qprec) <= 0.5
        if !any(mask)
            continue
        end
        neval[tid] += 1
        qspec_t = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
        best_cos = -Inf
        best_j = 0
        @inbounds for j in eachindex(r_specs)
            mask[j] || continue
            rspec_t = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
            merged = match_peaks_in_spectra_HRMS(qspec_t, rspec_t, window_size_matching)
            cosval = cosine_similarity_from_merged_HRMS(merged)
            if cosval > best_cos
                best_cos = cosval
                best_j = j
            end
        end
        if best_j != 0 && r_keys[best_j].id == q_keys[i].id
            ncorrect[tid] += 1
        end
    end
    n_eval = sum(neval)
    n_correct = sum(ncorrect)
    return 1.0 - (n_correct / n_eval)
end


function objective_cv_acc_no_filtering_HRMS(x)
    LET_threshold = x[1]
    noise_threshold = x[2]
    wf_int = x[3]
    wf_mz = x[4]
    window_size_matching = x[5]
    window_size_centroiding  = x[6]
    q_keys, q_precs, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    r_keys, r_precs, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    neval = zeros(Int, Threads.maxthreadid())
    ncorrect = zeros(Int, Threads.maxthreadid())
    @threads for i in eachindex(q_keys)
        tid = threadid()
        isempty(r_specs) && continue
        neval[tid] += 1
        qspec_t = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
        best_cos = -Inf
        best_j   = 0
        @inbounds for j in eachindex(r_specs)
            rspec_t = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
            merged = match_peaks_in_spectra_HRMS(qspec_t, rspec_t, window_size_matching)
            cosval = cosine_similarity_from_merged_HRMS(merged)
            if cosval > best_cos
                best_cos = cosval
                best_j   = j
            end
        end
        if best_j != 0 && r_keys[best_j].id == q_keys[i].id
            ncorrect[tid] += 1
        end
    end
    n_eval = sum(neval)
    n_correct = sum(ncorrect)
    return 1.0 - (n_correct / n_eval)
end


function objective_cv_MRR_filtering_HRMS(x)
    LET_threshold = x[1]
    noise_threshold = x[2]
    wf_int = x[3]
    wf_mz = x[4]
    window_size_matching = x[5]
    window_size_centroiding  = x[6]
    q_keys, q_precs, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])

    r_keys, r_precs, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    neval = zeros(Int, Threads.maxthreadid())
    rr_sum = zeros(Float64, Threads.maxthreadid())
    @threads for i in eachindex(q_keys)
        tid = threadid()
        qprec = q_precs[i]
        mask = @. abs(r_precs - qprec) <= 0.5
        if !any(mask)
            continue
        end
        neval[tid] += 1
        qspec_t = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
        scores = Float64[]
        ids    = String[]
        @inbounds for j in eachindex(r_specs)
            mask[j] || continue
            rspec_t = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
            merged = match_peaks_in_spectra_HRMS(qspec_t, rspec_t, window_size_matching)
            cosval = cosine_similarity_from_merged_HRMS(merged)
            push!(scores, cosval)
            push!(ids, r_keys[j].id)
        end
        perm = sortperm(scores, rev = true)
        qid = q_keys[i].id
        for (rank, k) in enumerate(perm)
            if ids[k] == qid
                rr_sum[tid] += 1.0 / rank
                break
            end
        end
    end
    n_eval = sum(neval)
    return 1.0 - (sum(rr_sum) / n_eval)
end


function objective_cv_MRR_no_filtering_HRMS(x)
    LET_threshold = x[1]
    noise_threshold = x[2]
    wf_int = x[3]
    wf_mz = x[4]
    window_size_matching = x[5]
    window_size_centroiding  = x[6]
    q_keys, q_precs, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    r_keys, r_precs, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol = :id, prec = :precursor_ion_mz, mzcol = :mz_ratio, intcol = :intensity, keep_ids = String[])
    neval = zeros(Int, Threads.maxthreadid())
    rr_sum = zeros(Float64, Threads.maxthreadid())
    @threads for i in eachindex(q_keys)
        tid = threadid()
        isempty(r_specs) && continue
        neval[tid] += 1
        qspec_t = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
        scores = Vector{Float64}(undef, length(r_specs))
        @inbounds for j in eachindex(r_specs)
            rspec_t = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
            merged = match_peaks_in_spectra_HRMS(qspec_t, rspec_t, window_size_matching)
            scores[j] = cosine_similarity_from_merged_HRMS(merged)
        end
        perm = sortperm(scores, rev = true)
        qid = q_keys[i].id
        for (rank, j) in enumerate(perm)
            if r_keys[j].id == qid
                rr_sum[tid] += 1.0 / rank
                break
            end
        end
    end
    n_eval = sum(neval)
    return 1.0 - (sum(rr_sum) / n_eval)
end


function get_scores_HRMS(df_query_raw::DataFrame, df_ref_raw::DataFrame; order, window_size_matching, window_size_centroiding, noise_threshold, LET_threshold, wf_mz, wf_int)
    q_keys, _, q_specs = build_spectra_by_id_prec_HRMS(df_query_raw; idcol=:id, prec=:precursor_ion_mz, mzcol=:mz_ratio, intcol=:intensity, keep_ids=String[])
    r_keys, _, r_specs = build_spectra_by_id_prec_HRMS(df_ref_raw; idcol=:id, prec=:precursor_ion_mz, mzcol=:mz_ratio, intcol=:intensity, keep_ids=String[])
    nq, nr = length(q_specs), length(r_specs)
    q_specs_T = Vector{typeof(q_specs[1])}(undef, nq)
    @threads for i in 1:nq
        q_specs_T[i] = apply_transformations_HRMS(copy(q_specs[i]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
    end
    r_specs_T = Vector{typeof(r_specs[1])}(undef, nr)
    @threads for j in 1:nr
        r_specs_T[j] = apply_transformations_HRMS(copy(r_specs[j]), spectrum_preprocessing_order; window_size_centroiding = window_size_centroiding, noise_threshold = noise_threshold, LET_threshold = LET_threshold, wf_mz = wf_mz, wf_int = wf_int)
    end
    S = Matrix{Float64}(undef, nq, nr)
    @threads for i in 1:nq
        qi = q_specs_T[i]
        @inbounds for j in 1:nr
            merged = match_peaks_in_spectra_HRMS(qi, r_specs_T[j], window_size_matching)
            S[i,j] = cosine_similarity_from_merged_HRMS(merged)
        end
    end
    q_ids = String[k.id for k in q_keys]
    r_ids = String[k.id for k in r_keys]
    df = DataFrame(S, Symbol.(r_ids); copycols=false, makeunique=true)
    insertcols!(df, 1, :query_id => q_ids)
    return df
end

