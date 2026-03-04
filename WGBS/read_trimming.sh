#!/bin/bash

#For the PBAT libraries, I need to do some trimming, because read2 had some problems

#samples=(MANC1 ML1G40 ML10G40 ML15G40) #Samples first sequenced (old library prep method)
#samples=(MANC1 MG5L6 MG5L7 MG5L10 MG5L14 MG5L15 MG5L17 MG5L19 MG5L20 MG20L6 MG20L7 MG20L10 MG20L14 MG20L15 MG20L17 MG20L19 MG20L20 ML6G40 ML7G40 ML10G40 ML14G40 ML15G40 ML17G40 ML19G40 ML20G40 MANC2a MG5L21 MG5L29 MG5L30 MG5L33 MG5L34 MG5L35 MG5L39 MG5L40 MG20L21 MG20L29 MG20L30 MG20L33 MG20L34 MG20L35 MG20L39 MG20L40 ML21G40 ML29G40 ML30G40 ML33G40 ML34G40 ML35G40 ML39G40 ML40G40) #Mat A samples with new library prep, G5, G20, and G40, missing L1 and L13, #Mat a samples with new library prep, G5, G20, and G40, missing L28 and L36
samples=(M1133 MANC1_R2 MANC1_R3 MANC2a_R2 MANC2a_R3 MG5L1 MG5L10_R2 MG5L13 MG5L28 MG5L36 MG20L1 MG20L10_R2 MG20L13 MG20L28 MG20L36 ML10G40_R2 ML13G40 ML28G40 ML36G40) 



#s=MG5L6

for s in ${samples[@]}
do
read_1=~/Genomics/Neurospora/methylation/raw/${s}/${s}_1.fq.gz
read_2=~/Genomics/Neurospora/methylation/raw/${s}/${s}_2.fq.gz
tread_1=~/Genomics/Neurospora/methylation/raw/${s}/trim.${s}_1.fq.gz
tread_2=~/Genomics/Neurospora/methylation/raw/${s}/trim.${s}_2.fq.gz
#sample_name=${s}

#Trim reads, and remove 10 bp from 5' end of read2
~/fastp/fastp -i $read_1 -I $read_2 -o $tread_1 -O $tread_2 --trim_front2=10

done

exit 0
