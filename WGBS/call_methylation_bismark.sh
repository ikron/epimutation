#!/bin/bash -l
#SBATCH --job-name=bismark_mapping
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=12G
#SBATCH --array=0-8
#SBATCH --mail-type=END
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

#--array=0-24

module load biokit

#Path to different programs
bismark=/projappl/project_2000350/Genomics/bismark
reference=/scratch/project_2000350/genomics/methylation/mreference/
#neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_Lambdaphage.fasta #Only ference directory is required for bismark
alignments=/scratch/project_2000350/genomics/methylation/alignments_bismark
outdir=/scratch/project_2000350/genomics/methylation/methylation/

#samples=(MANC1 MG5L6 MG5L7 MG5L10 MG5L14 MG5L15 MG5L17 MG5L19 MG5L20 MG20L6 MG20L7 MG20L10 MG20L14 MG20L15 MG20L17 MG20L19 MG20L20 ML6G40 ML7G40 ML10G40 ML14G40 ML15G40 ML17G40 ML19G40 ML20G40) #Mat A samples with new library prep, G5, G20, and G40, missing L1 and L13
#samples=(MANC2a MG5L21 MG20L21 ML21G40 MG5L29 MG20L29 ML29G40)
#samples=(MG5L30 MG20L30 ML30G40 MG5L33 MG20L33 ML33G40 MG5L34 MG20L34 ML34G40)
samples=(MG5L35 MG20L35 ML35G40 MG5L39 MG20L39 ML39G40 MG5L40 MG20L40 ML40G40)
#samples=(MANC2a MG5L21 MG5L29 MG5L30 MG5L33 MG5L34 MG5L35 MG5L39 MG5L40 MG20L21 MG20L29 MG20L30 MG20L33 MG20L34 MG20L35 MG20L39 MG20L40 ML21G40 ML29G40 ML30G40 ML33G40 ML34G40 ML35G40 ML39G40 ML40G40) #Mat a samples with new library prep, G5, G20, and G40, missing L28 and L36

s=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

#path to input file
input=$alignments/trim.${s}_1_bismark_bt2_pe.bam

$bismark/bismark_methylation_extractor -p --ignore 4 -o $outdir --gzip --CX --buffer_size 10G --comprehensive --cytosine_report --genome_folder $reference $input


