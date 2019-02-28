#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 1 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 01:00:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'ChIP_preprocess' # name of the job
#SBATCH -o 'ChIP_preprocess%j.out' # stdout log file
#SBATCH -e 'ChIP_preprocess%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.


#Go to sample-directory
cd ~/NAS/ChIP-ATAC-workshop-2019-course-material/${sample}

#Catenate the fastq-file
cat *.fastq.gz > ${sample}.fq.gz 

#Load the bowtie2 module:
module load GCC/7.3.0-2.30  OpenMPI/3.1.1 Bowtie2/2.3.4.2
bowtie2 -x /projects/fs1/common/genome/lunarc/indicies/bowtie2/mouse/mm10/mm10 -q ${sample}.fq.gz > ${sample}.sam
module purge

#Load the samtools module
module load GCC/6.4.0-2.28  OpenMPI/2.1.2 SAMtools/1.7

#Make a bam-file and filter for a MAPQ-score of 10:
samtools view -q 10 -bS ${sample}.sam > ${sample}.bam

#Sort the bam-file:
samtools sort -T $SNIC_TMP ${sample}.bam > ${sample}.s.bam
module purge

#Load the bedtools module.
module load GCC/5.4.0-2.26  OpenMPI/1.10.3 icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 ifort/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 BEDTools/2.26.0

#Remove blacklisted regions:
bedtools intersect -v -abam ${sample}.s.bam -b /home/jonun/Documents/bioinfo_resources/mm10.blacklist.bed > ${sample}.filt.bam
module purge

# Clean up
rm ${sample}.fq.gz ${sample}.sam ${sample}.bam ${sample}.s.bam
