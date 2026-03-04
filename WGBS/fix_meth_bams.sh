#!/bin/bash -l
#SBATCH --job-name=fix_bams
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=END
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit

#Paths to different programs and folders
reference=/scratch/project_2000350/genomics/methylation/mreference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_Lambdaphage.fasta
BisSNP=/projappl/project_2000350/Genomics/BisSNP
samplepath=/scratch/project_2000350/genomics/methylation/alignments
path_to_interval_list=/scratch/project_2000350/genomics/methylation/interval_list

samples=(MANC1 MG5L6 MG5L7 MG5L10 MG5L14 MG5L15 MG5L17 MG5L19 MG5L20 MG20L6 MG20L7 MG20L10 MG20L14 MG20L15 MG20L17 MG20L19 MG20L20 ML6G40 ML7G40 ML10G40 ML14G40 ML15G40 ML17G40 ML19G40 ML20G40) #Mat A samples with new library prep, G5, G20, and G40, missing L1 and L13

#Sort all files and index them again
for s in ${samples[@]}
do

#Sort the bamfile
samtools sort -o $samplepath/${s}.sorted.bam $samplepath/${s}.bam

#Make an index
samtools index $samplepath/${s}.sorted.bam

### Add or replace read-groups ###
#Note that sample names is used here in output and RGID!
picard AddOrReplaceReadGroups I=$samplepath/${s}.sorted.bam O=$samplepath/${s}.RG.bam SORT_ORDER=coordinate RGID=${s} RGLB=my_library RGPL=illumina RGSM=${s} RGPU=who_cares CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

##Remove temp files
rm $samplepath/${s}.sorted.bam

done

exit 0
