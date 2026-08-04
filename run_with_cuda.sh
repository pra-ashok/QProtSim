#!/bin/bash
# Wrapper script to run qprotsim with proper CUDA environment

export CUDA_PATH=/home/linuxuser/apps/local/miniforge3/envs/qprotsim
export CUDA_HOME=/home/linuxuser/apps/local/miniforge3/envs/qprotsim
export LD_LIBRARY_PATH=/home/linuxuser/apps/local/miniforge3/envs/qprotsim/lib:/home/linuxuser/apps/local/miniforge3/envs/qprotsim/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
export CPATH=/home/linuxuser/apps/local/miniforge3/envs/qprotsim/targets/x86_64-linux/include:$CPATH
export PATH=/home/linuxuser/apps/local/miniforge3/envs/qprotsim/bin:$PATH

eval "$(conda shell.bash hook)"
conda activate qprotsim

# Run the command
exec qprotsim "$@"
