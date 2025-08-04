#!/bin/bash

# Function to extract GPU ID with the most free memory
get_gpu_with_most_free_memory() {
    # Query NVIDIA SMI for GPU index and free memory, sorting by free memory in descending order
    local gpu_info=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits | sort -t, -k2 -nr)

    # Get the first line (GPU with the most free memory)
    local top_gpu_info=$(echo "$gpu_info" | head -n 1)

    # Extract the GPU index (first column before the comma)
    local gpu_id=$(echo $top_gpu_info | cut -d, -f1)

    echo $gpu_id
}

# Main script execution
# Get the GPU ID with the most free memory
gpu_id=$(get_gpu_with_most_free_memory)

# Set the CUDA_VISIBLE_DEVICES environment variable
export CUDA_VISIBLE_DEVICES=$gpu_id

# Optionally, print the GPU ID for verification
echo "CUDA_VISIBLE_DEVICES set to GPU ID: $CUDA_VISIBLE_DEVICES"

# Note: This export will affect only sub-processes started from this script.
# To make it available globally or in the current shell, you need to source this script or manually set the variable.

# Modify these according to your server environments
MAINPATH="$(cd "$(dirname "$0")" && pwd)/supporting/nufft/finufft/build"
export FINUFFT_LIB=$MAINPATH
#export FINUFFT_LIB=/home/leeh5/recon/multitaskingReconTool/supporting/nufft/finufft/build
export CUDA_LIB=/usr/local/cuda/lib64

export LD_LIBRARY_PATH=$FINUFFT_LIB:$CUDA_LIB:$LD_LIBRARY_PATH

matlab -desktop 