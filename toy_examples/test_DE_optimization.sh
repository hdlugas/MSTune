#!/bin/bash

cd ${PWD}/..
echo -e "\nTuning parameters via differential evolution...\n"


julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_DE_tuning_acc_no_cv_no_filtering_no_bootstrap_all_params_HRMS.txt \
  --platform HRMS \
  --optimization_method DE \
  --params_to_optimize all \
  --metric accuracy \
  --filter_reference_candidates false \
  --crossvalidation false \
  --bootstrap_query false \
  --max_steps 20 \
  --pop_size 30 \
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

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_DE_tuning_MRR_no_cv_no_filtering_no_bootstrap_LET_thresh_HRMS.txt \
  --platform HRMS \
  --optimization_method DE \
  --params_to_optimize LET_thresh \
  --metric MRR \
  --filter_reference_candidates false \
  --crossvalidation false \
  --bootstrap_query false \
  --max_steps 20 \
  --pop_size 30 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_DE_tuning_acc_cv_no_filtering_no_bootstrap_noise_thresh_HRMS.txt \
  --platform HRMS \
  --optimization_method DE \
  --params_to_optimize noise_thresh \
  --metric accuracy \
  --filter_reference_candidates false \
  --crossvalidation true \
  --bootstrap_query false \
  --max_steps 20 \
  --pop_size 30 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_DE_tuning_acc_cv_filtering_no_bootstrap_wf_int_HRMS.txt \
  --platform HRMS \
  --optimization_method DE \
  --params_to_optimize wf_int \
  --metric accuracy \
  --filter_reference_candidates true \
  --crossvalidation true \
  --bootstrap_query false \
  --max_steps 20 \
  --pop_size 30 \
  --LB_wf_int 0.0 \
  --UB_wf_int 1.0
 
julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_DE_tuning_acc_cv_filtering_bootstrap_wf_mz_HRMS.txt \
  --platform HRMS \
  --optimization_method DE \
  --params_to_optimize wf_mz \
  --metric accuracy \
  --filter_reference_candidates true \
  --crossvalidation true \
  --bootstrap_query true \
  --max_steps 20 \
  --pop_size 30 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 1.0
 
julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_DE_tuning_MRR_cv_filtering_no_bootstrap_ws_matching_HRMS.txt \
  --platform HRMS \
  --optimization_method DE \
  --params_to_optimize ws_matching \
  --metric MRR \
  --filter_reference_candidates true \
  --crossvalidation true \
  --bootstrap_query false \
  --max_steps 20 \
  --pop_size 30 \
  --LB_ws_matching 0.0 \
  --UB_ws_matching 0.5

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_DE_tuning_accuracy_no_cv_no_filtering_no_bootstrap_ws_centroiding_HRMS.txt \
  --platform HRMS \
  --optimization_method DE \
  --params_to_optimize ws_matching \
  --metric MRR \
  --filter_reference_candidates false \
  --crossvalidation false \
  --bootstrap_query false \
  --max_steps 20 \
  --pop_size 30 \
  --LB_ws_centroiding 0.0 \
  --UB_ws_centroiding 0.5



julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_DE_tuning_no_cv_all_params_NRMS.txt \
  --platform NRMS \
  --optimization_method DE \
  --params_to_optimize all \
  --metric accuracy \
  --crossvalidation false \
  --bootstrap_query false \
  --max_steps 20 \
  --pop_size 30 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_DE_tuning_cv_all_params_NRMS.txt \
  --platform NRMS \
  --optimization_method DE \
  --params_to_optimize all \
  --metric accuracy \
  --crossvalidation true \
  --n_folds 5 \
  --bootstrap_query false \
  --max_steps 20 \
  --pop_size 30 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_DE_tuning_no_cv_with_bootstrapped_query_all_params_NRMS.txt \
  --platform NRMS \
  --optimization_method DE \
  --params_to_optimize all \
  --metric accuracy \
  --crossvalidation false \
  --bootstrap_query true \
  --max_steps 20 \
  --pop_size 30 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0 \
  --LB_noise_thresh 0.0 \
  --UB_noise_thresh 1.0 \
  --LB_wf_int 0.0 \
  --UB_wf_int 5.0 \
  --LB_wf_mz 0.0 \
  --UB_wf_mz 5.0

julia --threads auto src/OptiMS.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_DE_tuning_no_cv_with_bootstrapped_query_LET_thresh_NRMS.txt \
  --platform NRMS \
  --optimization_method DE \
  --params_to_optimize LET_thresh \
  --metric accuracy \
  --crossvalidation false \
  --bootstrap_query true \
  --max_steps 20 \
  --pop_size 30 \
  --LB_LET_thresh 0.0 \
  --UB_LET_thresh 5.0

