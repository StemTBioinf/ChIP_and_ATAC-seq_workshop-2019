#!/bin/bash
#A wrapper that loops over the ATAC samples and calls preprocessing_ATAC.sh for each sample.

#samples="ATAC-seq_ETP_T_cells_rep1 ATAC-seq_ETP_T_cells_rep2 ATAC-seq_DN3_T_cells_rep1 ATAC-seq_DN3_T_cells_rep2" 
samples="ATAC-seq_DN3_T_cells_rep2"
for i in ${samples};do

	echo "sbatch --export=sample=${i} preprocessing_ATAC.sh"
	sbatch --export=sample=${i} preprocessing_ATAC.sh

done
