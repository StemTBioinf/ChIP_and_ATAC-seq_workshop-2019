#!/bin/bash
#A wrapper that loops over the samples and calls batch_maketagdirectory.sh for each sample.

samples="10340_DN1_H3K4me2 10500_DN3_H3K4me2 11046_DN1_H3K4me2 11307_DN1_PU1 11353_DN3_H3K4me2 11735_DN2b_PU1 DN1_input ThyDN3_input ATAC-seq_DN3_T_cells_rep1 ATAC-seq_DN3_T_cells_rep2 ATAC-seq_ETP_T_cells_rep1 ATAC-seq_ETP_T_cells_rep2"

for i in ${samples};do

	echo "sbatch --export=sample=${i} batch_maketagdirectory.sh"
	sbatch --export=sample=${i} batch_maketagdirectory.sh

done
