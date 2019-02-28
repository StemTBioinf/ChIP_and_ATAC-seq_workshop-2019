#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 2 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 03:00:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'ATAC_preprocess' # name of the job
#SBATCH -o 'ATAC_preprocess%j.out' # stdout log file
#SBATCH -e 'ATAC_preprocess%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ATAC samples.


#Go to sample-directory
cd ~/NAS2/ChIP-ATAC-workshop-2019-course-material/${sample}

#Catenate the fastq-file
cat *.R1.fastq.gz > ${sample}.R1.fq.gz 
cat *.R2.fastq.gz > ${sample}.R2.fq.gz

#Load the bowtie2 module:
module load GCC/7.3.0-2.30  OpenMPI/3.1.1 Bowtie2/2.3.4.2
bowtie2 -x /projects/fs1/common/genome/lunarc/indicies/bowtie2/mouse/mm10/mm10 -p 2 -1 ${sample}.R1.fq.gz -2 ${sample}.R2.fq.gz > ${sample}.sam
module purge

rm ${sample}.R1.fq.gz ${sample}.R2.fq.gz
#Load the samtools module
module load GCC/6.4.0-2.28  OpenMPI/2.1.2 SAMtools/1.7

#Make a bam-file and filter for a MAPQ-score of 10:
samtools view  -@ 2 -q 10 -bS ${sample}.sam > ${sample}.bam

#Sort the bam-file:
samtools sort -@ 2 -T $SNIC_TMP ${sample}.bam > ${sample}.s.bam

#Index the sorted file so we can filter out chromosome MT.
samtools index -@ 2 ${sample}.s.bam


rm ${sample}.sam ${sample}.bam
#Remove MT chromosomes:
chromosomes=$(echo $(for j in {1..19}; do echo "chr"$j; done | xargs) "chrX" "chrY")

samtools view -@ 2 -b ${sample}.s.bam ${chromosomes} > ${sample}.no_MT.bam
module purge

#rm ${sample}.s.bam
#Load the Picard module:
module load picard/2.15.0-Java-1.8.0_92
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${sample}.no_MT.bam O=${sample}.rmdup.bam M=${sample}.mrkdupmetric.txt REMOVE_DUPLICATES=true

#rm ${sample}.no_MT.bam
module purge

#Load the bedtools module.
module load GCC/5.4.0-2.26  OpenMPI/1.10.3 icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 ifort/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 BEDTools/2.26.0

#Remove blacklisted regions:
bedtools intersect -v -abam ${sample}.rmdup.bam -b /home/jonun/Documents/bioinfo_resources/mm10.blacklist.bed > ${sample}.filt.bam
module purge

#rm ${sample}.rmdup.bam

#Optional: The final file above is normally fine for peak analysis but in order to identify the exact transponson insertion sites we need to shift each position a little. This is because the Tn5 transosome binds as a 9-nt dimer. In order to correct for this all reads aligning to the + strand were offset by +4 bps, and all reads aligning to the - strand were offset -5 bps.
#module load GCC/5.4.0-2.26  OpenMPI/1.10.3 icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 ifort/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 BEDTools/2.26.0
#bedtools bamtobed -i ${sample}.filt.bam > ${sample}.filt.bed
#Use awk for the above described filtration:
#awk -F '\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print ${sample}.filt.bed} > ${sample}.tn.5_shifted.bed




