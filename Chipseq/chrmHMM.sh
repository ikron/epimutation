#!/bin/bash -l
#SBATCH --job-name=chrmHMM
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END
#SBATCH --mail-user=marvilla@jyu.fi


module load biojava/17

#Define input file
input=/scratch/project_2000350/genomics/chipseq/bwa-mem
input_bed=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/chrmHMM
#control
control=/scratch/project_2000350/genomics/chipseq/bwa-mem
#Define the output directory
output=/scratch/project_2000350/genomics/chipseq/chrmHMM_bin
output_model=/scratch/project_2000350/genomics/chipseq/chrmHMM

chromosomelengthfile=/scratch/project_2000350/genomics/Neurospora_reference/chrm_size
cellmarkfiletable=/scratch/project_2000350/genomics/chipseq/scripts/cellmarkfiletable.txt
program=/scratch/project_2000350/genomics/chipseq/programs/ChromHMM
assembly=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
modelfile=/scratch/project_2000350/genomics/chipseq/chrmHMM/emissions_2modelfile=/scratch/project_2000350/genomics/chipseq/chrmHMM/emissions_2.txt



#java -mx4000M -jar $program/ChromHMM.jar BinarizeBam -b 200 -c $control $chromosomelengthfile $input $cellmarkfiletable $output

java -mx4000M -jar $program/ChromHMM.jar BinarizeBed -b 200 -center  $chromosomelengthfile $input_bed $cellmarkfiletable $output

java -mx4000M -jar $program/ChromHMM.jar LearnModel -b 200  $output $output_model 3 $assembly

java -mx4000M -jar $program/ChromHMM.jar MakeSegmentation $modelfile $output $output_model

