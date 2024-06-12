#!/bin/bash -l

path2mod=/projappl/project_2000350/Genomics/dist #Path to modkit
alignments=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alignments_modified #input files
output=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alignments_modified/bedMethyl #output
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_lambaDNA.fasta #reference

#List of sample names
samples=(2489_MatA 2489_Mata G1_L11 G1_L23 G1_L31 G1_L2 G1_L25 G1_L5 G5_L11 G5_L23 G5_L31 G5_L2 G5_L25 G5_L5 G7_L11 G7_L23 G7_L5 G7_L2 G7_L25 G7_L31 G10_L11 G10_L25 G10_L2 G10_L31 G10_L23 G10_L5 G15_L31 G15_L23 G15_L11 G15_L25 G15_L2 G15_L5 G8_L23 G8_L5 G8_L11 G8_L25 G8_L2 G8_L31)
sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

#Run modkit pileup
$path2mod/modkit pileup $alignments/$sample.bam $output/$sample.bed --ref $reference --only-tabs --threads 4 
