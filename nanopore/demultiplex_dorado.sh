#!/bin/bash -l

#Demultiplex dorado
path2dorado=/projappl/project_2000350/Genomics/dorado-0.4.1-linux-x64/bin
inputdir=/scratch/project_2000350/genomics/nanopore/240207_ncrassa/basecalled_modified
outputdir=/scratch/project_2000350/genomics/nanopore/240207_ncrassa/demultiplex_dorado

$path2dorado/dorado demux --kit-name SQK-NBD114-24 --output-dir $outputdir $inputdir/basecalled_5mC_240207.sam
