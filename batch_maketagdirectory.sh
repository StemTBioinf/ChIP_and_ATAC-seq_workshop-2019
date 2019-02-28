#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 4 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 00:20:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'makeTagDirectories' # name of the job
#SBATCH -o 'makeTagDirectories%j.out' # stdout log file
#SBATCH -e 'makeTagDirectories%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.


#Go to sample-directory
cd ~/NAS2/ChIP-ATAC-workshop-2019-course-material/${sample}
#Load the samtools module
module load GCC/7.3.0-2.30 SAMtools/1.9
#First create a sam-file
samtools view -@ 4 -h -o ${sample}.sam ${sample}.filt.bam
module purge
#Load the Homer module:
module load GCC/7.3.0-2.30 homer/4.10
#Create a tagdirectory
makeTagDirectory ${sample}_tagdir ${sample}.sam
rm ${sample}.sam
