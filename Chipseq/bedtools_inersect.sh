#!/bin/bash -l
#SBATCH --job-name=bedtools_intersect
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2G
#SBATCH --array=0-11
#SBATCH --mail-type=END
#SBATCH --mail-user=marvilla@jyu.fi

#Does DMRS intersect with H3k9?

samples=(G5_L11 G8_L11 G10_L11 G5_L5 G8_L5 G10_L5 G5_L25 G8_L25 G10_L25 G5_L31 G8_L31 G10_L31)


sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based


module load bedtools/2.30.0

dirmeth=/scratch/project_2000350/genomics/nanopore/MAlines_5mC/chrmHMM
dirchip=/scratch/project_2000350/genomics/chipseq/macs
dirout=/scratch/project_2000350/genomics/chipseq/intersect_res

#the sites that overlap 100%
#bedtools intersect -b $dirchip/${sample}_2broad.bed -a $dirmeth/$sample.bed -f 1  > $dirout/${sample}_100.bed

#the sites that overlap 80%
#bedtools intersect -a $dirmeth/$sample.bed -b $dirchip/${sample}_2broad.bed  -f 0.80  > $dirout/${sample}_80.bed
#bedtools intersect -b $dirout/${sample}_100.bed -a $dirout/${sample}_80.bed -v  > $dirout/${sample}_only80.bed

# sites that overlap 50%
#bedtools intersect -b $dirchip/${sample}_2broad.bed -a $dirmeth/$sample.bed -f 0.50  > $dirout/${sample}_50.bed
#bedtools intersect -b $dirout/${sample}_80.bed -a $dirout/${sample}_50.bed -v  > $dirout/${sample}_only50.bed


#bedtools intersect -b $dirout/${sample}_50.bed -a $dirmeth/$sample.bed -v  > $dirout/${sample}_less50.bed

#bedtools intersect -a $dirmeth/$sample.bed -b $dirchip/${sample}_2broad.bed -v  > $dirout/${sample}_nooverlap.bed

#bedtools intersect -b $dirout/${sample}_nooverlap.bed -a $dirout/${sample}_less50.bed -v  > $dirout/${sample}_onlyless50.bed

bedtools intersect -b $dirmeth/$sample.bed -a $dirchip/${sample}_2broad.bed -v  > $dirout/${sample}_H3K9_nooverlap.bed
