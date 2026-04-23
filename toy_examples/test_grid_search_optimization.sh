#!/bin/bash

cd ${PWD}/..
echo -e "\nTuning parameters via exhaustive grid-search...\n"

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_grid_tuning_NRMS.txt \
  --platform NRMS \
  --optimization_method grid \
  --params_to_optimize all \
  --metric MRR \
  --n_grid_points 3 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_grid_tuning_HRMS.txt \
  --platform HRMS \
  --optimization_method grid \
  --params_to_optimize all \
  --metric accuracy \
  --crossvalidation false \
  --filter_reference_candidates false \
  --bootstrap_query false \
  --n_grid_points 2 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0 \
  --LB_ws_matching 0.0 \
  --UB_ws_matching 0.5 \
  --LB_ws_centroiding 0.0 \
  --UB_ws_centroiding 0.5

