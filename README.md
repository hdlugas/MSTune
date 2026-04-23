# OptiMS
OptiMS is a Julia-based command-line tool for tuning parameters involved in preprocessing (i) nominal-resolution mass spectrometry (NRMS) data such as gas chromatography - mass spectrometry data and (ii) high-resolution mass spectrometry (HRMS) data such as liquid chromatography - tandem mass spectrometry data. The three main functionalities of OptiMS are (i) identifying optimal parameters via differential evolution, (ii) identifying optimal parameters via exhaustive grid-based search, and (iii) simply running compound identification and recording similarity scores.


## Table of Contents
- [1. Install dependencies](#install-dependencies)
- [2. Parameter descriptions](#param-descriptions)
- [3. Functionality](#functionality)
   - [3.1 Optimize parameters via differential evolution](#DE)
   - [3.2 Optimize parameters via exhausive grid search](#grid-search)
   - [3.3 Run compound identification](#run-compound-identification)


<a name="install-dependencies"></a>
## 1. Install dependencies
The Julia packages required to run OptiMS are BlackBoxOptim, CSV, DataFrames, LinearAlgebra, Random, Statistics, and StatsBase. These dependencies can be installed with the Julia command:

```
using Pkg; Pkg.add(["BlackBoxOptim", "CSV", "DataFrames", "LinearAlgebra","Random", "Statistics", "StatsBase"])
```

<a name="param-descriptions"></a>
# 2. Parameter descriptions
The following two spectrum preprocessing transformations are offered:

-   Weight Factor Transformation: Given a pair of user-defined weight
    factor parameters $(\text{a,b})$ and spectrum $I$ with m/z values
    $(m_{1},m_{2},...,m_{n})$ and intensities $(x_{1},x_{2},...,x_{n})$,
    the transformed spectrum $I^{\star}$ has the same m/z values as $I$
    and has intensities given by
    $I^{\star}:=(m_{1}^{\text{a}}\cdot x_{1}^{\text{b}},m_{2}^{\text{a}}\cdot x_{2}^{\text{b}},...,m_{n}^{\text{a}}\cdot x_{n}^{\text{b}})$.

-   Noise Removal: Given a user-defined noise removal parameter $r$ and
    a spectrum $I$ with intensities $(x_{1},x_{2},...,x_{n})$, noise
    removal removes peaks from $I$ with
    $x_{j}< r\cdot\text{max}(\{x_{1},x_{2},...,x_{n}\})$ for
    $j\in\{1,2,...,n\}$.

-   Low-Entropy Transformation: Given a user-defined low-entropy
    threshold parameter $T$ and spectrum $I$ with intensities
    $(x_{1},x_{2},...,x_{n})$, $\sum_{i=1}^nx_i = 1$, and Shannon
    entropy $H_{Shannon}(I)=-\sum_{i=1}^{n}x_{i}\cdot ln(x_{i})$, the
    transformed spectrum intensities
    $I^{\star}=(x_{1}^{\star},x_{2}^{\star},...,x_{n}^{\star})$ are such
    that, for all $i\in\{1,2,...,n\}$, $x_{i}^{\star}=x_{i}$ if
    $H_{Shannon}(I)\geq T$ and
    $x_{i}^{\star}=x_{i}^{\frac{1+H_{Shannon}(I)}{1+T}}$ if
    $H_{Shannon}(I)<T$.

-   Matching (only applicable to HRMS data): Given a user-defined
    window-size parameter $w_{matching}$ and two spectra $I$, $J$ with
    m/z ratios $(a_{1},a_{2},...,a_{n}), (b_{1},b_{2},...,b_{m})$ and
    intensities $(x_{1},x_{2},...,x_{n}), (y_{1},y_{2},...,y_{m})$,
    respectively, of which we would like to measure the similarity
    between the matching procedure outputs two spectra
    $I^{\star},J^{\star}$ containing the same number of peaks with
    $I^{\star}$ and $J^{\star}$ having intensities and
    identical m/z ratios. Specifically, for a given peak $(a_{i},x_{i})$
    of $I$, if there are no peaks $(b_{j},y_{j})$ in $J$ with
    $|a_{i}-b_{j}|< w_{matching}$, then the peak $(a_{i},x_{i})$
    remains in $I^{\star}$ and the peak $(a_{i},0)$ is included in
    $J^{\star}$. If there is at least one peak $(b_{j},y_{j})$ with
    $|a_{i}-b_{j}|< w_{matching}$, then the peak $(a_{i},x_{i})$
    remains in $I^{\star}$ and the peak
    $(a_{i},\sum_{j\text{ such that }|a_{i}-b_{j}|< w_{matching}}b_{j})$
    is included in $J^{\star}$. This procedure is applied when
    transposing the roles of $I$ and $J$ as well.

-   Centroiding (only applicable to HRMS data): Given a user-defined
    window-size parameter $w_{centroiding}$ and a spectrum $I$ with m/z
    values $(m_{1},m_{2},...,m_{n})$ and intensities
    $(x_{1},x_{2},...,x_{n})$, the transformed spectrum $I^{\star}$
    merges adjacent peaks $(m_{i},x_{i}),(m_{i+1},x_{i+1})$ into the
    peak
    $(\frac{m_{i}\cdot x_{i}+m_{i+1}\cdot x_{i+1}}{x_{i}+x_{i+1}},x_{i}+x_{i+1})$
    if $|m_{i}-m_{i+1}|< w_{centroiding}$ for
    $i\in\{1,2,...,n-1\}$. This centroiding procedure generalizes to
    more than two peaks whose m/z values are within a distance
    $w_{centroiding}$ of each other.



Thus, there are two weight factor parameters, one noise removal threshold parameter, one low-entropy threshold parameter, and in the case of HRMS, two window-size parameters one can tweak.

Given a pair of processed spectra intensities
$I=(a_{1},a_{2},...,a_{n}), J=(b_{1},b_{2},...,b_{n})\in\mathbb{R}^{n}$
with $0\leq a_{i},b_{i}\leq 1$ for all $i\in\{1,2,...,n\}$ and
$\sum_{i=1}^{n}a_{i}=\sum_{i=1}^{n}b_{i}=1$, OptiMS computes the Cosine similarity between the two spectra:
```math
S_{Cosine}(I,J)=\frac{I\circ J}{|I|_{2}\cdot |J|_{2}}
```
where multiplication in the numerator refers to the dot product $I\circ J=a_{1}b_{1}+a_{2}b_{2}+...+a_{n}b_{n}$ of $I$ and $J$ and multiplication in the denominator refers to multiplication of the $L^{2}$-norms of $I$ and $J$, $\vert I\vert_{2}=\sqrt{a_{1}^{2}+a_{2}^{2}+...+a_{n}^{2}}, \vert J\vert_{2}=\sqrt{b_{1}^{2}+b_{2}^{2}+...+b_{n}^{2}}$.


<a name="functionality"></a>
## 3. Functionality

OptiMS has three main capabilities:
1. Tune parameters to maximize either cross-validated accuracy or cross-validated mean reciprocal rank (MRR) via differential evolution optimization.
2. Compute the cross-validated accuracy or cross-validated MRR for every set of parameters in a user-defined grid of possible parameters.
3. Compute the similarity scores between every query spectrum and every reference spectrum for a single user-defined choice of parameters.

Scripts which run toy examples illustrating each of these three methods are provided. These toy examples can be run by navigating to the necessary directory and executing the scripts:
```
cd toy_examples
./test_DE_optimization.sh
./test_grid_optimization.sh
./test_run_compound_identification.sh
```

To view the OptiMS usage instruction, one can run the following from the command-line (once the necessary Julia dependencies are installed):
```
julia src/OptiMS.jl --help
```

The complete usage instructions for OptiMS are:
```
Usage:
  julia OptiMS.jl --query_data <string> --reference_data <string> --output <string> --platform <string> --optimization_method <string> --metric <string> --crossvalidation <boolean> --filter_reference_candidates <boolean> --n_folds <int> --bootstrap_query <boolean> --random_seed <int> --params_to_optimize <string> --spectrum_preprocessing_order <string> --LB_LET_thresh <float> --UB_LET_thresh <float> --LB_noise_thresh <float> --UB_noise_thresh <float> --LB_wf_int <float> --UB_wf_int <float> --LB_wf_mz <float> --UB_wf_mz <float> --LET_thresh <float> --noise_thresh <float> --wf_int <float> --wf_mz <float> --threads <int> --max_steps <int> --pop_size <int> --n_grid_points <int>


Arguments:
  --query_data                    Path to input TXT file of query dataset (required).
  --reference_data                Path to input TXT file of reference dataset (required).
  --output                        Path to output TXT file (required).
  --platform                      Chromatography platform (required, options=[HRMS,NRMS]).
  --optimization_method           Optimization approach (optional, options=[DE,grid,none], default=DE).
  --metric                        Quantity to maximize in the objective function (optional, options=[accuracy,MRR], default=accuracy).
  --crossvalidation               Boolean indicating whether or not to perform 5-fold cross-validation inside the objective function (optional, options=[true,false], default=true).
  --filter_reference_candidates   Boolean indicating whether or not to filter on precursor ion mass/charge; only applicable to HRMS data (optional, options=[true,false], default=true).
  --n_folds                       Number folds to use for cross-validation. Only applicable for crossvalidation is 'true' and optimization_method is not 'none' (optional, default=5).
  --bootstrap_query               Boolean indicating whether or not to construct query dataset of same size from resampling query spectra with replacement; only useful for computing confidence intervals of parameter estimates (optional, options=[true,false], default=false).
  --random_seed                   Random seed to be used in all computations with a stochastic component (optional, default=1).
  --params_to_optimize            String denoting the parameters to optimize (optional, options='all','LET_thresh','noise_thresh','wf_int','wf_mz', default='all').
  --spectrum_preprocessing_order  String denoting the order of spectrum preprocessing transformations; format must be a string with 2-6 characters chosen from C, F, M, N, L, W representing centroiding, filtering based on mass/charge and intensity values, matching, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively; note that if an argument is passed for HRMS, then 'M' must be contained in the argument since matching is a required preprocessing step in spectral library matching of HRMS data; furthermore, 'C' must be performed before matching since centroiding can change the number of ion fragments in a given spectrum; note that C and M are not applicable to NRMS data (optional, default: NCMWL for HRMS and NWL for NRMS').
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
```

<a name="DE"></a>
### 3.1 Optimize parameters via differential evolution
To identify optimal parameters to maximize the metric (e.g. accuracy in this case) using differential evolution optimization with user-specified parameter bounds and maximum number of steps, one can run:
```
julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_DE_tuning.txt \
  --platform NRMS \
  --optimization_method DE \
  --params_to_optimize all \
  --metric accuracy \
  --crossvalidation true \
  --max_steps 20 \
  --pop_size 50 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0
```

<a name="grid-search"></a>
### 3.2 Optimize parameters via exhaustive grid search
To record the metric (e.g. MRR in this case) for each combination of parameters in a user-specified grid of parameters with user-specified parameter bounds, one can run:
```
julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_grid_tuning.txt \
  --platform NRMS \
  --optimization_method grid \
  --params_to_optimize all \
  --metric MRR \
  --n_grid_points 2 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0
```

<a name="run-compound-identification"></a>
### 3.3 Run compound identification
To simply run compound identification and record all similarity scores with user-specified parameters, one can run:
```
julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_similarity_scores.txt \
  --platform NRMS \
  --optimization_method none \
  --LET_thresh 3.0 \
  --noise_thresh 0.1 \
  --wf_int 1.5 \
  --wf_mz 0.5
```

