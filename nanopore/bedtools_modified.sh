#!/bin/bash -l

module load bedtools/2.30.0

alignments=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/alignments_modified/bedMethyl #input files
output=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/bedMethyl_input #output
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_lambaDNA.fasta #reference
chrm_list=/scratch/project_2000350/genomics/Neurospora_reference/chrm_list #chromosome list
int=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/bedMethyl_int

#List of sample names
samples=(2489_MatA 2489_Mata G1_L11 G1_L23 G1_L31 G1_L2 G1_L25 G1_L5 G5_L11 G5_L23 G5_L31 G5_L2 G5_L25 G5_L5 G7_L11 G7_L23 G7_L5 G7_L2 G7_L25 G7_L31 G10_L11 G10_L25 G10_L2 G10_L31 G10_L23 G10_L5 G15_L31 G15_L23 G15_L11 G15_L25 G15_L2 G15_L5 G8_L23 G8_L5 G8_L11 G8_L25 G8_L2 G8_L31)
sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based


#Separate + and - strands and filter for a coverage of mmin 5
awk 'BEGIN {FS="\t";OFS="\t"};{if ($6 == "+" && $10 > 5) print $0}' $alignments/$sample.bed >  $int/$sample-foward.bed
awk 'BEGIN {FS="\t";OFS="\t"};{if ($6 == "-" && $10 > 5) print $0}' $alignments/$sample.bed >  $int/$sample-reverse.bed
#print the columns of interest and the end position to extract from fasta
awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$2+3,$3,$5,$6,$12,$13}' $int/$sample-foward.bed > $int/$sample-foward1.bed
awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$3-3,$3,$3,$5,$6,$12,$13}' $int/$sample-reverse.bed > $int/$sample-reverse1.bed

#concatenate files
cat $int/$sample-foward1.bed $int/$sample-reverse1.bed > $int/$sample-1.bed
#filter out the starts that are smaller that 0 and the coverage above 5
awk 'BEGIN {FS="\t";OFS="\t"}; {if ($2 > 1 && $5 > 5) print $0}' $int/$sample-1.bed  > $int/$sample-2.bed
bedtools sort -faidx $chrm_list -i $int/$sample-2.bed > $int/$sample-3.bed

#get trinucleotide with getfasta
bedtools getfasta -fi $reference -bed $int/$sample-3.bed -s -bedOut > $int/$sample-4.bed
#print in the right order to be like bismark output
awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$4,$6,$7,$8,$9,$9}' $int/$sample-4.bed > $int/$sample-5.bed

awk 'BEGIN {FS="\t";OFS="\t"}; { gsub(/CAA/, "CHH", $6); gsub(/CAG/, "CHG", $6); 
                                 gsub(/CAT/, "CHH", $6); gsub(/CAC/, "CHH", $6);
                               gsub(/CGA/, "CG", $6); gsub(/CGG/, "CG", $6);
                               gsub(/CGT/, "CG", $6); gsub(/CGC/, "CG", $6);
                               gsub(/CTA/, "CHH", $6); gsub(/CTG/, "CHG", $6);
                               gsub(/CTT/, "CHH", $6); gsub(/CTC/, "CHH", $6);
                               gsub(/CCA/, "CHH", $6); gsub(/CCG/, "CHG", $6);
                               gsub(/CCT/, "CHH", $6); gsub(/CCC/, "CHH", $6);
                               print }' $int/$sample-5.bed > $output/$sample.bed


