#!/bin/bash
#SBATCH -A p30041
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --mem-per-cpu=4G
#SBATCH -t 04:00:00
#SBATCH --job-name="rvs-3hpa"
#SBATCH --output=rvs-3hpa.txt

ulimit -c 0
module load python/anaconda3.6
source activate px-env
    
python 01_pathway_discovery_V4.0.py 1




