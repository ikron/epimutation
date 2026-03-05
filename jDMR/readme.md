# Inferring differentially methylated regions (DMRs)

Inferring DMRs were done using jDMR

The file jDMR_mod.R contains some modifications to the original jDMR package to make it work on the CSC computer cluster

The scripts for processing WGBS data are in folder ./jDMR/WGBS and contain files
<ul>
  <li> jDMR_all.R and its batch control file jDMR_job.sh for inferring DMRs across the genome </li>
  <li> jDMR_euchromatin.R and jDMR_euchromatin_job.sh for inferring DMRs in euchromatic regions </li>
  <li> jDMR_centromericH3K9.R and jDMRcentromericH3K9_job.sh for inferring DMRs in centromeric regions </li>
  <li> jDMR_H3K27_ex_H3K9.R and jDMR_H3K27_ex_H3K9_job.sh for inferring DMRs in regions marked by H3K27 but excluding overlapping H3K9 regions</li>
  <li> jDMR_H3K9_ex_cent2.R and jDMR_H3K9_ex_cent2_job.sh for inferring DMRs in regions marked by H3K9 excluding centromeric H3K9 (i.e. interspersed heterochromatic regions)</li>
</ul>

The scripts for processing Nanopore data are in folder ./jDMR/Nanopore and contain files
<ul>
  <li> jDMR_all_nanopore.R and jDMR_nanopore_job.sh for inferring DMRs across the genome </li>
  <li> jDMR_methimpute_euchromatin.R and jDMR_euchromatin.sh for inferring DMRs across different regions of the genome (note that this file includes multiple regions) </li>
  <li> jDMR_methimpute_H3K9.R and jDMR_methimpute_H3K9.sh for inferring DMRs across different regions of the genome (file includes multiple regions) </li>
</ul>
