#!/bin/bash
#A wrapper that loops over the ChIP samples and calls preprocessing_ChIP.sh for each sample.

samples="10340_DN1_H3K4me2 10500_DN3_H3K4me2 11046_DN1_H3K4me2 11307_DN1_PU1 11353_DN3_H3K4me2 11735_DN2b_PU1 DN1_input ThyDN3_input"

for i in ${samples};do

	echo "sbatch --export=sample=${i} preprocessing_ChIP.sh"
	sbatch --export=sample=${i} preprocessing_ChIP.sh

done
