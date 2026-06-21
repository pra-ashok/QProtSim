#!/bin/bash
# Wrapper script to run protqc with proper CUDA environment

export CUDA_PATH=/home/linuxuser/apps/local/miniforge3/envs/protqc
export CUDA_HOME=/home/linuxuser/apps/local/miniforge3/envs/protqc
export LD_LIBRARY_PATH=/home/linuxuser/apps/local/miniforge3/envs/protqc/lib:/home/linuxuser/apps/local/miniforge3/envs/protqc/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
export CPATH=/home/linuxuser/apps/local/miniforge3/envs/protqc/targets/x86_64-linux/include:$CPATH
export PATH=/home/linuxuser/apps/local/miniforge3/envs/protqc/bin:$PATH

eval "$(conda shell.bash hook)"
conda activate protqc

# Run the command
exec protqc "$@"
