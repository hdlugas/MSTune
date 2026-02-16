#!/bin/bash

cd ${PWD}/..
echo -e "\nTuning parameters via differential evolution...\n"

julia --threads auto src/OptiMS.jl \
  --query_data toy_example/data/query_data.txt \
  --reference_data toy_example/data/reference_data.txt \
  --output toy_example/output_DE_tuning_no_cv.txt \
  --optimization_method DE \
  --params_to_optimize all \
  --metric accuracy \
  --crossvalidation false \
  --bootstrap_query false \
  --max_steps 50 \
  --pop_size 50 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_example/data/query_data.txt \
  --reference_data toy_example/data/reference_data.txt \
  --output toy_example/output_DE_tuning_cv.txt \
  --optimization_method DE \
  --params_to_optimize all \
  --metric accuracy \
  --crossvalidation true \
  --n_folds 5 \
  --bootstrap_query false \
  --max_steps 50 \
  --pop_size 50 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_example/data/query_data.txt \
  --reference_data toy_example/data/reference_data.txt \
  --output toy_example/output_DE_tuning_no_cv_with_bootstrapped_query_all_params.txt \
  --optimization_method DE \
  --params_to_optimize all \
  --metric accuracy \
  --crossvalidation false \
  --bootstrap_query true \
  --max_steps 50 \
  --pop_size 50 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_example/data/query_data.txt \
  --reference_data toy_example/data/reference_data.txt \
  --output toy_example/output_DE_tuning_no_cv_with_bootstrapped_query_LET_thresh.txt \
  --optimization_method DE \
  --params_to_optimize LET_thresh \
  --metric accuracy \
  --crossvalidation false \
  --bootstrap_query true \
  --max_steps 50 \
  --pop_size 50 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0

