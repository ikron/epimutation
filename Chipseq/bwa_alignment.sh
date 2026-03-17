#!/bin/bash -l
#SBATCH --job-name=chipseq_bwa
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0-25
#SBATCH --mail-type=END
#SBATCH --mail-user=marvilla@jyu.fi

module load biokit

#samples 
samples=(A1_rep1_H3K9me3_S1_L001_R1_001 A1_rep2_H3K9me3_S2_L001_R1_001 A2_rep1_H3K9me3_S3_L001_R1_001 A2_rep1_H3K9me3_S4_L001_R1_001 A3_rep1_H3K9me3_S5_L001_R1_001 A4_rep1_H3K9me3_S6_L001_R1_001 M1_H3K9me3_S7_L001_R1_001 M2_H3K9me3_S8_L001_R1_001 M3_H3K9me3_S9_L001_R1_001 M4_H3K9me3_S10_L001_R1_001 M5_H3K9me3_S11_L001_R1_001 M6_H3K9me3_S12_L001_R1_001 M7_H3K9me3_S13_L001_R1_001 M8_H3K9me3_S14_L001_R1_001 M9_H3K9me3_S15_L001_R1_001 M10_H3K9me3_S16_L001_R1_001 M11_H3K9me3_S17_L001_R1_001 M12_H3K9me3_S18_L001_R1_001 M13_H3K9me3_S19_L001_R1_001 M14_H3K9me3_S20_L001_R1_001 M15_H3K9me3_S21_L001_R1_001 M16_H3K9me3_S22_L001_R1_001 M17_H3K9me3_S23_L001_R1_001 M18_H3K9me3_S24_L001_R1_001 M19_H3K9me3_S25_L001_R1_001 M20_H3K9me3_S26_L001_R1_001)

sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

#Define inut file
input=/scratch/project_2000350/genomics/chipseq/raw
#Define reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
#Define the output sample
output=/scratch/project_2000350/genomics/chipseq/bwa-mem


bwa mem $reference $input/$sample.fastq.gz | samtools sort > $output/$sample.bam 



