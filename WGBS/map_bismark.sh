#!/bin/bash -l
#SBATCH --job-name=bismark_mapping
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0
#SBATCH --mail-type=END
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit

#Path to different programs
bismark=/projappl/project_2000350/Genomics/bismark
reference=/scratch/project_2000350/genomics/methylation/mreference/
#neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_Lambdaphage.fasta #Only ference directory is required for bismark

#samples=(MANC1 MG5L6 MG5L7 MG5L10 MG5L14 MG5L15 MG5L17 MG5L19 MG5L20 MG20L6 MG20L7 MG20L10 MG20L14 MG20L15 MG20L17 MG20L19 MG20L20 ML6G40 ML7G40 ML10G40 ML14G40 ML15G40 ML17G40 ML19G40 ML20G40) #Mat A samples with new library prep, G5, G20, and G40, missing L1 and L13
#samples=(MANC2a MG5L21 MG5L29 MG5L30 MG5L33 MG5L34 MG5L35 MG5L39 MG5L40 MG20L21 MG20L29 MG20L30 MG20L33 MG20L34 MG20L35 MG20L39 MG20L40 ML21G40 ML29G40 ML30G40 ML33G40 ML34G40 ML35G40 ML39G40 ML40G40) #Mat a samples with new library prep, G5, G20, and G40, missing L28 and L36
#samples=(MANC1_R2 MANC1_R3 MANC2a_R2 MANC2a_R3 MG5L1 MG5L10_R2 MG5L13 MG5L28 MG5L36 MG20L1 MG20L10_R2 MG20L13 MG20L28 MG20L36 ML10G40_R2 ML13G40 ML28G40 ML36G40)
samples=(ML1G40)

s=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

read1=/scratch/project_2000350/genomics/methylation/raw/trim.${s}_1.fq.gz
read2=/scratch/project_2000350/genomics/methylation/raw/trim.${s}_2.fq.gz

#Use bismark to map bisulfite reads to the genome

$bismark/bismark -X 700 --dovetail --genome $reference -1 $read1 -2 $read2
