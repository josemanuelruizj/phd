#! /bin/bash
#SBATCH --partition=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/users/genomics/anabel/logs/main-%j.out
#SBATCH --error=/users/genomics/anabel/logs/main-%j.e

module load R/4.1.2

Rscript /users/genomics/anabel/scripts/mainscript.R $1 $2 $3 $4

