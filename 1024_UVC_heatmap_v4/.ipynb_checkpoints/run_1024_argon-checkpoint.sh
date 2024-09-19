#!/bin/bash
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -l walltime=70:00:00
#PBS -N run_1024_argon
#PBS -j oe
#PBS -o /rds/general/user/sk1021/home/1024_UVC_heatmap_v4/run_1024_argon.log

# load the necessary modules
module load anaconda3/personal
module load cuda/11.4.2

# activate the conda environment
eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate simulation

# Set the environment variable to disable file validation
export PYDEVD_DISABLE_FILE_VALIDATION=1

# navigate to the directory containing the notebook
cd /rds/general/user/sk1021/home/1024_UVC_heatmap_v4

# run notebook
python -Xfrozen_modules=off -m jupyter nbconvert --to notebook --execute 1024_argon.ipynb --output 1024_argon_executed.ipynb

# deactivate the conda environment
conda deactivate
