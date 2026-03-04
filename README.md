# epimutation
Investigating epimutations in _Neurospora crassa_

We have previously performed a mutation accumulation experiment using the filamentous fungus _Neurospora crassa_. In this study we investigated spontaneous methylation changes in the MA-lines. First we performed bisulphite sequencing to determine methylated cytosines in these lines. We used the program Methimpute to to determine methylation status of each cytosine. Then we determined methylated regions in the genome using jDMR. To estimate rates of spontaneous methylation changes we ran the program AlphaBeta for both single cytosines and differentially methylated regions (DMRs).

Sequencing data related to this project has been deposited to European nucleotide archive, accession number: [PRJEB108830](https://www.ebi.ac.uk/ena/browser/view/PRJEB108830) and ChIP-seq data has been deposited to NCBI Gene Expression Omnibus database, accession number: [GSE313506](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313506)

The processing of different types of sequencing data are described below.

Note that the original scripts were mostly run on the Finnish CSC computer cluster that uses the SLURM batch job system. As such to make the scripts work on your computing environment you may need to modify the control commands accordingly. As well as input and ouput folders, program directories etc. of course.

## Processing bisulphite sequencing data

Scripts for processing bisulfite-seq data are in the folder WGBS

## Processing Nanopore sequencing data

Scripts for processing Nanopore sequencing data are in the folder nanopore


## Processing ChIP-seq

1. BWA.sh : Script to align the ChiP-seq data
   
3. macs2_advanced_sample.sh : Script to run macs2 (advanced version) with the metrics specific to each sample.
   
5. chrhmm : Script to run run chrmHMM (https://ernstlab.github.io/ChromHMM)
   
7. ChromTime : modifed script of ChromTime in python2 (https://github.com/ernstlab/ChromTime)
