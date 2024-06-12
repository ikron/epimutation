# epimutation
Investigating epimutations in Neurospora crassa

In the nanopore folder the scripts for basecalling, demultiplex, mapping and methylation calling:
Run  in the following order
1.-basecalling_modified_bases_dorado.sh : Script for nanopore modified-base basecalling using Dorado. We used the model dna_r10.4.1_e8.2_400bps_sup@v4.2.0 v2. 
2.-demultiplex_dorado.sh : Script for demutliplex in the basecalled file using dorado.
3.-mapping_dorado.sh : Script for mapping to reference genome using dorado. 
4.-modkit_modified.sh : Script for calling methylated bases using modkit
5.-bedtools_modified.sh : Script to modify modkit output to make it similar to Bismark.
6.-methimpute.R : R script to run methimpute
