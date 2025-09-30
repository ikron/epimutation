# epimutation
Investigating epimutations in _Neurospora crassa_

We have previously performed a mutation accumulation experiment using the filamentous fungus _Neurospora crassa_. In this study we investigated spontaneous methylation changes in the MA-lines. First we performed bisulphite sequencing to determine methylated cytosines in these lines. We used the program Methimpute to to determine methylation status of each cytosine. Then we determined methylated regions in the genome using jDMR. To estimate rates of spontaneous methylation changes we ran the program AlphaBeta for both single cytosines and differentially methylated regions (DMRs).

## Processing bisulphite sequencing data

### Running Methimpute

### Running jDMR


## Processing nanopore data

In the nanopore folder the scripts for basecalling, demultiplex, mapping and methylation calling:
Run  in the following order:

1.basecalling_modified_bases_dorado.sh : Script for nanopore modified-base basecalling using Dorado. We used the model dna_r10.4.1_e8.2_400bps_sup@v4.2.0 v2. 

2.demultiplex_dorado.sh : Script for demutliplex in the basecalled file using dorado.

3.mapping_dorado.sh : Script for mapping to reference genome using dorado. 

4.modkit_modified.sh : Script for calling methylated bases using modkit

5.bedtools_modified.sh : Script to modify modkit output to make it similar to Bismark.

6.methimpute.R : R script to run methimpute

## Processing ChIP-seq

1. BWA.sh : Script to align the ChiP-seq data
   
3. macs2_advanced_sample.sh : Script to run macs2 (advanced version) with the metrics specific to each sample.
   
5. chrhmm : Script to run run chrmHMM (https://ernstlab.github.io/ChromHMM)
   
7. ChromTime : modifed script of ChromTime in python2 (https://github.com/ernstlab/ChromTime)
