#!/bin/bash

cd ${PWD}/..
echo -e "\nRunning compound identification...\n"

julia --threads auto src/OptiMS.jl \
  --query_data toy_example/data/query_data.txt \
  --reference_data toy_example/data/reference_data.txt \
  --output toy_example/output_similarity_scores.txt \
  --optimization_method none \
  --LET_thresh 3.0 \
  --noise_thresh 0.1 \
  --wf_intensity 1.5 \
  --wf_mz 0.5

