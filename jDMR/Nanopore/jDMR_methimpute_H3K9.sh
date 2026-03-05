#!/bin/bash -l
#SBATCH --job-name=jDMR_H3K9_R
#SBATCH --account=project_2000350
#SBATCH --output=/scratch/project_2000350/genomics/nanopore/sbatch_out/output_%A_%a.txt
#SBATCH --error=/scratch/project_2000350/genomics/nanopore/sbatch_out/errors_%A_%a.txt
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marvilla@jyu.fi

# Load r-env-singularity
module load r-env-singularity

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2000350" >> ~/.Renviron

# Run the R script
srun singularity_wrapper exec Rscript --no-save jDMR_methimpute_H3K9.R 

