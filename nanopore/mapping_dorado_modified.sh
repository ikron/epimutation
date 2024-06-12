#!/bin/bash -l
#Define input file
input=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/demultiplex_dorado
#Define reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_lambaDNA.fasta
#Define the output folder
output=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alignments_modified
#Index the reference genome
lra index -ONT /scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_lambaDNA.fasta
#Path to Dorado
path2dorado=/projappl/project_2000350/Genomics/dorado-0.4.1-linux-x64/bin

#all_samples 
samples=(2489_MatA 2489_Mata G1_L11 G1_L23 G1_L31 G1_L2 G1_L25 G1_L5 G5_L11 G5_L23 G5_L31 G5_L2 G5_L25 G5_L5 G7_L11 G7_L23 G7_L5 G7_L2 G7_L25 G7_L31 G10_L11 G10_L25 G10_L2 G10_L31 G10_L23 G10_L5 G15_L31 G15_L23 G15_L11 G15_L25 G15_L2 G15_L5 G8_L23 G8_L5 G8_L11 G8_L25 G8_L2 G8_L31)
sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

### using dorado
$path2dorado/dorado aligner $reference $input/$sample.bam -t 4 > $output/$sample.sam

