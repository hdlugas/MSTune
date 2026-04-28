#!/bin/bash

cd ${PWD}/..
echo -e "\nRunning compound identification...\n"

julia --threads auto src/MSTune.jl \
  --query_data toy_examples/data/gcms_query.txt \
  --reference_data toy_examples/data/gcms_reference.txt \
  --output toy_examples/output_similarity_scores_NRMS.txt \
  --platform NRMS \
  --optimization_method none \
  --LET_thresh 3.0 \
  --noise_thresh 0.1 \
  --wf_intensity 1.5 \
  --wf_mz 0.5

julia --threads auto src/MSTune.jl \
  --query_data toy_examples/data/lcmsms_query.txt \
  --reference_data toy_examples/data/lcmsms_reference.txt \
  --output toy_examples/output_similarity_scores_HRMS.txt \
  --platform HRMS \
  --optimization_method none \
  --LET_thresh 3.0 \
  --noise_thresh 0.1 \
  --wf_intensity 1.5 \
  --wf_mz 0.5 \
  --ws_matching 0.25 \
  --ws_centroiding 0.48

