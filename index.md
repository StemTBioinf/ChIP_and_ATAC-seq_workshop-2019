# ChIP and ATAC-seq workshop 2019

### -by Jonas "your-friendly-neighbourhood-chromatin-companion" Ungerbäck
##########################################################################
### Workshop information:

In this workshop you will learn to work with ChIP and ATAC-sequencing data, mainly via the command line interface. You will learn to call peaks, looking at overlapping peaks between different datasets, annotate peaks, motif analysis with Homer and the MEME suite and creating bigwigs and heatmaps with deepTools.

--------------------------------------------------------------
##### These are the tools you primarily will use for the workshop:

HOMER: http://homer.ucsd.edu/homer/

BedTools: https://bedtools.readthedocs.io/en/latest/

MACS2: https://github.com/taoliu/MACS

MEME: http://meme-suite.org/

DeepTools: https://deeptools.readthedocs.io/en/develop/

--------------------------------------------------------------
--------------------------------------------------------------
### The datasets
You will be working with ChIP -and ATAC-datasets from T-cell development in mouse. If you are interested you can read the original papers here:

https://www.ncbi.nlm.nih.gov/pubmed/22500808 (ChIP-seq datasets)
https://www.ncbi.nlm.nih.gov/pubmed/29466756 (ATAC-seq datasets)

#### The datasets are as follows:
##### ChIP (single-end sequencing data)
The transcription factor PU.1 as well as the active promoter/enhancer H3K4me2 marker.

#### ATAC (paired-end sequencing data)
ProT cells pre (ETPs) and post (DN3) T-cell lineage commitment

You will recieve data in the bed-format but the command for the preprocessing  can be found in the following scripts:

`preprocessing_ChIP.sh` `ChIP_wrapper.sh` `preprocessing_ATAC.sh` `ATAC_wrapper.sh` `batch_maketagdirectory.sh` `makeTag_wrapper.sh`.

#
--------------------------------------------------------------

### First and foremost:

We will use aurora and not Lsens2 for this workshops so after login into Aurora (aurora.lunarc.lu.se), go to the folder you will be working in /projects/fs3/_youruser_.

Copy the following folder to you user:

```bash
cp -r /tmp/ChIP-ATAC-workshop-2019-course-material .
#Go to the folder.
cd ChIP-ATAC-workshop-2019-course-material
```
###############################################################
###############################################################
### ChIP-seq assignments
**1.** We will user Homer to call ChIP-peaks but before we do that we have to create a file-format Homer can work with. This is called a Tag Directory.

Start by loading the required modules.
```bash
module load GCC/7.3.0-2.30 homer/4.10
```
Run the following to check what options makeTagDirectory has
```bash
makeTagDirectory
```
We will now use Homer's makeTagDirectory to create tag directories from all the filtered bam-files (including ATAC-samples since this will be used later). Sadly this will not immediately work since Homer requires samtools to this. Batch scripting to the rescue.

```bash
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
```
Create a wrapper script that loop over the samples.
```bash
#!/bin/bash
#A wrapper that loops over the samples and calls batch_maketagdirectory.sh for each sample.

samples="10340_DN1_H3K4me2 10500_DN3_H3K4me2 11046_DN1_H3K4me2 11307_DN1_PU1 11353_DN3_H3K4me2 11735_DN2b_PU1 DN1_input ThyDN3_input ATAC-seq_DN3_T_cells_rep1 ATAC-seq_DN3_T_cells_rep2 ATAC-seq_ETP_T_cells_rep1 ATAC-seq_ETP_T_cells_rep2"

for i in ${samples};do

        echo "sbatch --export=sample=${i} batch_maketagdirectory.sh"
        sbatch --export=sample=${i} batch_maketagdirectory.sh

done
```
Call the script
```bash
bash makeTag_wrapper.sh
```



Since this takes a while I have prepared it so in each folder you will find a subfolder called a something-something tagdir :). Alternatively, to circumvent samtools at this step, we could have created a bed-file and then made the tagdir from that.
You can read more about tag directories here:
http://homer.ucsd.edu/homer/ngs/tagDir.html

**We are now ready to call ChIP-peaks.**

**2.** **Peak finding**
Run the following to check what options findPeaks has:
```bash
module load GCC/7.3.0-2.30 homer/4.10
findPeaks
```
As you can see findPeaks can be customized with a large variety of options but by setting `-style` to either `factor` or `histone`, presets several options that are usually fine.



```bash
findPeaks 10340_DN1_H3K4me2/10340_DN1_H3K4me2_tagdir/ -style histone -o auto -i DN1_input/DN1_input_tagdir/
```
This will create a file in the tagdirectory called `peaks.txt` or `regions.txt`.

**Task: Do this for all ChIP-samples (differentiate between histones and factors). I would recommend using a sbatch script.**

How many peaks/regions did we get for each sample? **Hint:** Check the peak/region files in each tagdir with `head peak_file`. A good praxis is also to rename the peak-file to something descriptive.

> 10340_DN1_H3K4me2_regions.txt: Total peaks = 45085
> 11046_DN1_H3K4me2_regions.txt: Total peaks = 42673
> 10500_DN3_H3K4me2_regions.txt: Total peaks = 28009
> 11353_DN3_H3K4me2_regions.txt: Total peaks = 32967
> 11307_DN1_PU1_peaks.txt:       Total peaks = 51993
> 11735_DN2b_PU1_peaks.txt:      Total peaks = 17777

This is a good start and but no matter how good a peak caller is, extra filtering is often required. If replicates exist, IDR (https://github.com/nboley/idris) is the state-of-the-art method for peak reproducibility. It is however relatively cumbersome to use and outside the scope of this one-day-course. **What to do instead?**

**Task:**
**If NO replicates:** Filter for a Homer Normalized peak-score (column #6). Empirically I have found that when TF-ChIP a score of ≥ 10-15 is fine. The lower the score the higher the likelihood that we are picking up noise.
**Problem:** The hashtags (#) on the first lines of the Homer peak file. Possible to remove in Excel but faster to use the combined power of `sed` and `awk`:

Example:  `sed '/^#/d' 11735_DN2b_PU1_peaks.txt|awk 'BEGIN {FS=OFS="\t"} ; {if ($6 >= 15) print }' - > 11735_DN2b_PU1_peaks_15.txt `
$6 is the column where Homer stores the Normalized peak score in a peak file.

> 11307_DN1_PU1_peaks_15.txt : 37476 peaks
> 11735_DN2b_PU1_peaks_15.txt: 6688 peaks

**If replicates:** Later we will use BedTools but for now we will use Homer's `mergePeaks` in combination with the `-prefix` option.

Type:
```bash
mergePeaks
```
Create a new folder to merge the replicates in from the H3K4me2 samples and copy the peak files here:
```bash
mergePeaks 10340_DN1_H3K4me2_regions.txt 11046_DN1_H3K4me2_regions.txt -prefix DN1
	Max distance to merge: direct overlap required (-d given)
	Merging peaks...
	Comparing 10340_DN1_H3K4me2_regions.txt (45085 total) and 10340_DN1_H3K4me2_regions.txt (45085 total)
	Comparing 10340_DN1_H3K4me2_regions.txt (45085 total) and 11046_DN1_H3K4me2_regions.txt (42673 total)
	Comparing 11046_DN1_H3K4me2_regions.txt (42673 total) and 10340_DN1_H3K4me2_regions.txt (45085 total)
	Comparing 11046_DN1_H3K4me2_regions.txt (42673 total) and 11046_DN1_H3K4me2_regions.txt (42673 total)

10340_DN1_H3K4me2_regions.txt	11046_DN1_H3K4me2_regions.txt	Total	Name
	X	4530	11046_DN1_H3K4me2_regions.txt
X		7874	10340_DN1_H3K4me2_regions.txt
X	X	36998	10340_DN1_H3K4me2_regions.txt|11046_DN1_H3K4me2_regions.txt
```
```bash
mergePeaks 10500_DN3_H3K4me2_regions.txt 11353_DN3_H3K4me2_regions.txt -prefix DN3
	Max distance to merge: direct overlap required (-d given)
	Merging peaks...
	Comparing 10500_DN3_H3K4me2_regions.txt (28009 total) and 10500_DN3_H3K4me2_regions.txt (28009 total)
	Comparing 10500_DN3_H3K4me2_regions.txt (28009 total) and 11353_DN3_H3K4me2_regions.txt (32967 total)
	Comparing 11353_DN3_H3K4me2_regions.txt (32967 total) and 10500_DN3_H3K4me2_regions.txt (28009 total)
	Comparing 11353_DN3_H3K4me2_regions.txt (32967 total) and 11353_DN3_H3K4me2_regions.txt (32967 total)

10500_DN3_H3K4me2_regions.txt	11353_DN3_H3K4me2_regions.txt	Total	Name
	X	8823	11353_DN3_H3K4me2_regions.txt
X		3223	10500_DN3_H3K4me2_regions.txt
X	X	24031	10500_DN3_H3K4me2_regions.txt|11353_DN3_H3K4me2_regions.txt
```
Rename the overlapping peak files to something more convenient.

Collect all peak-files in one folder.


**3.** **Annotation of peak files:**
```bash
annotatePeaks.pl
```
This is Homer's workhorse so take some time to understand the options:
http://homer.ucsd.edu/homer/ngs/annotation.html and
http://homer.ucsd.edu/homer/ngs/quantification.html
For now we will just use the basic proximity based annotation. This requires a built-in Homer database (needs to be separately installed). In this case we will use _mm10_.

```bash
#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 4 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 00:20:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'annotation_mm10' # name of the job
#SBATCH -o 'annotation_mm10%j.out' # stdout log file
#SBATCH -e 'annotation_mm10%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.

#Load Homer module
module load GCC/7.3.0-2.30 homer/4.10
#Basic peak-filea annotation.
annotatePeaks.pl 11307_DN1_PU1_peaks_15.txt mm10 > 11307_DN1_PU1_peaks_15.mm10_annotated.txt

annotatePeaks.pl 11735_DN2b_PU1_peaks_15.txt mm10 > 11735_DN2b_PU1_peaks_15.mm10_annotated.txt

annotatePeaks.pl DN1_H3K4me2_regions.overlapping.txt mm10 > DN1_H3K4me2_regions.overlapping.mm10_annotated.txt

annotatePeaks.pl DN3_H3K4me2_regions.overlapping.txt mm10 > DN3_H3K4me2_regions.overlapping.mm10_annotated.txt

module purge
```
Copy the files to your own computer, open and look at them. Be warned, if you open the file in Excel some gene names will change to dates.

**Questions to think about/homework:**
* What are the ratios of the genomic regions bound? Create a pie-chart.
* What's the average peak-distance to a transcription-start-site (TSS)?

**4.** **Motif analysis**

When we have the peak file, basic motif-analysis in Homer is easy. Motif analysis is the reason Homer was developed and it is therefore very good at it with manually curated background models. To get all the options for basic motif-analysis, type:
`findMotifsGenome.pl`

**Task:** Perform motif analysis on the PU.1 ChIP-files. **WARNING!** A motif analysis can take a considerable amount of time if you have many peaks. Allocate accordingly.
```bash
#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 4 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 00:30:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'findMotifsGenome_mm10' # name of the job
#SBATCH -o 'findMotifsGenome_mm10%j.out' # stdout log file
#SBATCH -e 'findMotifsGenome_mm10%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.

#Load Homer module
module load GCC/7.3.0-2.30 homer/4.10
#Basic peak-filea annotation.
findMotifsGenome.pl 11307_DN1_PU1_peaks_15.mm10_annotated.txt mm10 ./Motifs-11307_DN1_PU1_peaks_15 -size 200 -mask -p 4 -preparsedDir /home/jonun/NAS2/ChIP-ATAC-workshop-2019-course-material/homer_background_preparsed_dir

findMotifsGenome.pl 11735_DN2b_PU1_peaks_15.mm10_annotated.txt mm10 ./Motifs-11735_DN2b_PU1_peaks_15 -size 200 -mask -p 4 -preparsedDir /home/jonun/NAS2/ChIP-ATAC-workshop-2019-course-material/homer_background_preparsed_dir

module purge
```
The `-preparsedDir` option is needed since we do not have write permissions in the Homer genome dir we need to specify where Homer builds the background model.

**Open the motif-folders and compare knownResults.html and homerResults.html. What is the difference and why?**

**What is the most enriched motif?**

_Motif analysis will take awhile so let's move on while we wait._

Another very common tool-suite is MEME. This can be run online or via the command-line. There are many tools in the meme suite but we will use **meme-chip** which is designed for ChIP data with thousands of sequences.

The input for meme is fasta-formatted sequences. This can be created from the peak file. MEME also works best when the sequences have the same length. Homer's peak size-estimate can be found in the tagInfo.txt-file in the tag directory.

> Let's use 200 bp as the estimate.

Load the Homer module:
```bash
module load GCC/7.3.0-2.30 homer/4.10
```

Adjust the PU.1 peak files to 200 bp:
```bash
adjustPeakFile.pl 11307_DN1_PU1_peaks_15.txt -size 200 > 11307_DN1_PU1_peaks_15_200.txt

adjustPeakFile.pl 11735_DN2b_PU1_peaks_15.txt -size 200 > 11735_DN2b_PU1_peaks_15_200.txt
```
`homeTools extract` can convert a peak file to a fasta-file.

```bash
homerTools extract 11307_DN1_PU1_peaks_15_200.txt /projects/fs1/common/genome/lunarc/genomes/mouse/mm10/mm10.fa -fa > 11307_DN1_PU1_peaks_15_200.fa

homerTools extract 11735_DN2b_PU1_peaks_15_200.txt /projects/fs1/common/genome/lunarc/genomes/mouse/mm10/mm10.fa -fa > 11735_DN2b_PU1_peaks_15_200.fa
```
**NOTE:** This can also be done with bedtools getfasta after converting the peak file to a bed file but more on that later.

Create and run a batch-script that calls chip-meme.
```bash
#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 1 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 02:00:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'chip_meme' # name of the job
#SBATCH -o 'chip_meme_mm10%j.out' # stdout log file
#SBATCH -e 'chip_meme%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.

#Load meme module
module load GCC/7.3.0-2.30  OpenMPI/3.1.1 MEME/5.0.4
#Run meme-chip on the fasta files
meme-chip -db /projects/fs3/jonun/ChIP-ATAC-workshop-2019-course-material/meme_motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme -meme-nmotifs 10 -o meme_11735_DN2b_PU1_peaks_15_200  11735_DN2b_PU1_peaks_15_200.fa

meme-chip -db /projects/fs3/jonun/ChIP-ATAC-workshop-2019-course-material/meme_motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme -meme-nmotifs 10 -o meme_11307_DN1_PU1_peaks_15_200  11307_DN1_PU1_peaks_15_200.fa

module purge
```
**OR**
these .fa-files can be uploaded to the online tool:
http://meme-suite.org/tools/meme-chip

Open the result-folder and study the output.

**Homework:** There are many more tools in the MEME suite. Explore them and see what they can do for you!

http://meme-suite.org/doc/overview.html?man_type=web

_Motif analysis will take awhile so let's move on while we wait._

### ATAC-seq assignments
We have now looked at some basic ChIP-seq analysis and we will come back to that again later. One "data type" that is very similar to ChIP-seq is ATAC-seq. In fact, it is analyzed in such a similar manner that an absolute majority of the analyzes are perfomed with ChIP-seq tools. ATAC-peak calling is, however, complicated by a few facts:

* ATAC-seq peaks are much more varying in size than intra-sample ChIP-peaks.
* Lack of normalization/input controls.
* Tn5 transposome sequence bias.
* Partly the large abundance of mitochondrial reads but correct pre-processing deals with this qute well.

In theory most peak callers can be used to find ATAC-peaks though they may require several rounds of peak calling and merging of the outputs. The most widely used is MACS2: https://github.com/taoliu/MACS

We will now use MACS2 to call peaks from our filtered bam ATAC-output. There are many ways to tweak ATAC-peak calling with MACS2 but we will for now stick to ENCODE's guidelines:

https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit

**1. MACS2 ATAC-peak calling**

From the MACS2 manual:

_Here are some examples for combining --shift and --extsize:
Option 1: To find enriched cutting sites such as some DNAse-Seq datasets. In this case, all 5' ends of sequenced reads should be extended in both direction to smooth the pileup signals. If the wanted smoothing window is 200bps, then use '--nomodel --shift -75 --extsize 150'._

_Option 2: For certain nucleosome-seq data, we need to pileup the centers of nucleosomes using a half-nucleosome size for wavelet analysis (e.g. NPS algorithm). Since the DNA wrapped on nucleosome is about 147bps, this option can be used: '--nomodel --shift 37 --extsize 73'._

**Homework:** Re-run the command below with option two to see how it performs against option 1.

First create a batch-script that perform the operation:

```bash
#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora projec
#SBATCH -n 1 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 00:30:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'macs2_callpeaks' # name of the job
#SBATCH -o 'macs2_callpeaks%j.out' # stdout log file
#SBATCH -e 'macs2_callpeaks%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.

#Load the MACS2 module
module load GCC/4.9.3-2.25  OpenMPI/1.10.2 icc/2015.3.187-GNU-4.9.3-2.25  impi/5.0.3.048 ifort/2015.3.187-GNU-4.9.3-2.25  impi/5.0.3.048 MACS2/2.1.0.20150731-Python-2.7.11

sample_path="/home/jonun/NAS2/ChIP-ATAC-workshop-2019-course-material"

#Call peaks
macs2 callpeak -t ${sample_path}/${sample}/${sample}.filt.bam -g mm -f AUTO -n MACS2_${sample} --nomodel --shift -75 --extsize 150 --keep-dup all --outdir ${sample_path}/${sample}/MACS2_${sample} -B --SPMR

module purge
```

Then create script that calls the previous script for each sample:
```bash
#!/bin/bash
#A wrapper that loops over the ATAC samples and calls batch_macs2callpeaks.sh for each sample.

samples="ATAC-seq_ETP_T_cells_rep1 ATAC-seq_ETP_T_cells_rep2 ATAC-seq_DN3_T_cells_rep1 ATAC-seq_DN3_T_cells_rep2"

for i in ${samples};do

        echo "sbatch --export=sample=${i} batch_macs2callpeaks.sh"
        sbatch --export=sample=${i} batch_macs2callpeaks.sh

done
```
Run this:
```bash
bash MACS2_wrapper.sh
```
**Discussion point:** It is quite common to add the --broad option to macs2 callpeaks. How would this affect the output?

**Note:** This is a more efficient way to setup the peak calling compared to how we called findPeaks with Homer. What is the difference?

**2.** **Subtract overlapping peaks between replicates using bedtools.**

Like in the ChIP-case we want to use reproducible peaks. Again, IDR is the best practice but we are gonna work with overlapping beaks between replicats.

Start by copying each .narrow peak-file into a common folder.

Load the bedtools module.
```bash
module load GCC/5.4.0-2.26  OpenMPI/1.10.3 icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 ifort/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 BEDTools/2.26.0
```
Like Homer, this is a suite that contain many usefool tools for processing of files in the bed-format and you can read about them here:
https://bedtools.readthedocs.io/en/latest/#

For now we are gonna use `bedtools intersect`.
Look at the options with:
```bash
bedtools intersect -h
```
```bash
bedtools intersect -a MACS2_ATAC-seq_ETP_T_cells_rep1_peaks.narrowPeak -b MACS2_ATAC-seq_ETP_T_cells_rep2_peaks.narrowPeak > ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed

bedtools intersect -a MACS2_ATAC-seq_DN3_T_cells_rep1_peaks.narrowPeak -b MACS2_ATAC-seq_DN3_T_cells_rep2_peaks.narrowPeak > ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed
```
How many peaks do the overlaps contain?
The number of lines in a file can be derived with the command `cat file | wc -l`.

> ETP: 56994
>
> DN3: 20616

This can of course be interpreted as there are much fewer open sites at the DN3 stage but it could also mean that the DN3 samples suffer from a lower signal-to-noise.



**3.** **Cell type specific peak subset analysis.**

We will now derive the peaks that are unique to ETPs and DN3, respectively, and investigate them further.

To only report entries in `-a` that is not in `-b` we will use the `-v` option with `bedtools intersect`.

```bash
bedtools intersect -a ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed -b ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed -v > ETP_unique_ATAC_peaks.bed

bedtools intersect -b ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed -a ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed -v > DN3_unique_ATAC_peaks.bed
```

An alternative way to do the analysis above would be to calculate statistical difference with, for instance, DESeq2 or edgeR. This would allow for detection of dynamic changes in chromatin accessibility whereas the analysis above gives a more black and white picture of peaks that are completely lost or gained between the two stages.

Load the Homer module and annotate these peak files.

```bash
module purge
module load GCC/7.3.0-2.30 homer/4.10

annotatePeaks.pl ETP_unique_ATAC_peaks.bed mm10 > ETP_unique_ATAC_peaks.mm10_annotated.txt

annotatePeaks.pl DN3_unique_ATAC_peaks.bed mm10 > DN3_unique_ATAC_peaks.mm10_annotated.txt
```
_Open and look at the files._


**4.** **Motif analysis on subset specific ATAC-peaks.**

**Task:** What motifs are enriched in the two different categories?
```bash
#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora projec
#SBATCH -n 4 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 01:00:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'findMotifsGenome_mm10' # name of the job
#SBATCH -o 'findMotifsGenome_mm10%j.out' # stdout log file
#SBATCH -e 'findMotifsGenome_mm10%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.

#Load Homer module
module load GCC/7.3.0-2.30 homer/4.10
#Basic peak-filea annotation.
findMotifsGenome.pl ETP_unique_ATAC_peaks.mm10_annotated.txt mm10 ./Motifs-ETP_unique_ATAC_peaks -size 200 -mask -p 4 -preparsedDir /home/jonun/NAS2/ChIP-ATAC-workshop-2019-course-material/homer_background_preparsed_dir

findMotifsGenome.pl DN3_unique_ATAC_peaks.mm10_annotated.txt mm10 ./Motifs-DN3_unique_ATAC_peaks -size 200 -mask -p 4 -preparsedDir /home/jonun/NAS2/ChIP-ATAC-workshop-2019-course-material/homer_background_preparsed_dir

module purge
```
Open the homerResults.html and find out :)


### Combined ChIP- and ATAC-seq analyses


**1. We will now try to answer a common question. Can we identify TF binding in accessible/active and inaccessible/inactive chromatin and how do they differ?**

In our case the questions would be: **How many PU.1 peaks are there in ATAC-accessible reagions? In H3K4me2 marked regions?**

These question can be answered with both Homer and bedtools intersect but for now let us use Homer's `mergePeaks`.

PU.1 overlapping with H3K4me2 marked histones.
```bash
mergePeaks 11307_DN1_PU1_peaks_15.txt DN1_H3K4me2_regions.overlapping.txt -prefix ETP_H3K4me2
	Max distance to merge: direct overlap required (-d given)
	Merging peaks...
	Comparing 11307_DN1_PU1_peaks_15.txt (37476 total) and 11307_DN1_PU1_peaks_15.txt (37476 total)
	Comparing 11307_DN1_PU1_peaks_15.txt (37476 total) and DN1_H3K4me2_regions.overlapping.txt (36998 total)
	Comparing DN1_H3K4me2_regions.overlapping.txt (36998 total) and 11307_DN1_PU1_peaks_15.txt (37476 total)
	Comparing DN1_H3K4me2_regions.overlapping.txt (36998 total) and DN1_H3K4me2_regions.overlapping.txt (36998 total)

11307_DN1_PU1_peaks_15.txt	DN1_H3K4me2_regions.overlapping.txt	Total	Name
	X	18031	DN1_H3K4me2_regions.overlapping.txt
X		13367	11307_DN1_PU1_peaks_15.txt
X	X	18967	11307_DN1_PU1_peaks_15.txt|DN1_H3K4me2_regions.overlapping.txt
```
```bash
mergePeaks 11735_DN2b_PU1_peaks_15.txt DN3_H3K4me2_regions.overlapping.txt -prefix DN3_H3K4me2
	Max distance to merge: direct overlap required (-d given)
	Merging peaks...
	Comparing 11735_DN2b_PU1_peaks_15.txt (6688 total) and 11735_DN2b_PU1_peaks_15.txt (6688 total)
	Comparing 11735_DN2b_PU1_peaks_15.txt (6688 total) and DN3_H3K4me2_regions.overlapping.txt (24031 total)
	Comparing DN3_H3K4me2_regions.overlapping.txt (24031 total) and 11735_DN2b_PU1_peaks_15.txt (6688 total)
	Comparing DN3_H3K4me2_regions.overlapping.txt (24031 total) and DN3_H3K4me2_regions.overlapping.txt (24031 total)

11735_DN2b_PU1_peaks_15.txt	DN3_H3K4me2_regions.overlapping.txt	Total	Name
	X	20729	DN3_H3K4me2_regions.overlapping.txt
X		3175	11735_DN2b_PU1_peaks_15.txt
X	X	3302	11735_DN2b_PU1_peaks_15.txt|DN3_H3K4me2_regions.overlapping.txt
```
PU.1 overlapping with ATAC-accessible chromatin.

```bash
mergePeaks 11307_DN1_PU1_peaks_15.txt ../MACS2_peak_files/ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed -prefix ETP_ATAC
	Max distance to merge: direct overlap required (-d given)
	Merging peaks...
	Duplicate peak name (MACS2_ATAC-seq_ETP_T_cells_rep1_peak_40) - this could potentially cause problems
		Sometimes unavoidable for BED/2DBED formats
		New name for this peak is MACS2_ATAC-seq_ETP_T_cells_rep1_peak_40--2
		Warning over 1000 peaks with duplicate names
	Comparing 11307_DN1_PU1_peaks_15.txt (37476 total) and 11307_DN1_PU1_peaks_15.txt (37476 total)
	Comparing 11307_DN1_PU1_peaks_15.txt (37476 total) and ../MACS2_peak_files/ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed (56994 total)
	Comparing ../MACS2_peak_files/ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed (56994 total) and 11307_DN1_PU1_peaks_15.txt (37476 total)
	Comparing ../MACS2_peak_files/ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed (56994 total) and ../MACS2_peak_files/ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed (56994 total)

11307_DN1_PU1_peaks_15.txt	../MACS2_peak_files/ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed	Total	Name
	X	35479	../MACS2_peak_files/ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed
X		15783	11307_DN1_PU1_peaks_15.txt
X	X	21458	11307_DN1_PU1_peaks_15.txt|../MACS2_peak_files/ATAC-seq_ETP_T_cells_rep1_2.overlapping.bed
```
```bash
mergePeaks 11735_DN2b_PU1_peaks_15.txt ../MACS2_peak_files/ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed -prefix DN3_ATAC
	Max distance to merge: direct overlap required (-d given)
	Merging peaks...
	Duplicate peak name (MACS2_ATAC-seq_DN3_T_cells_rep1_peak_10) - this could potentially cause problems
		Sometimes unavoidable for BED/2DBED formats
		New name for this peak is MACS2_ATAC-seq_DN3_T_cells_rep1_peak_10--2
	Comparing 11735_DN2b_PU1_peaks_15.txt (6688 total) and 11735_DN2b_PU1_peaks_15.txt (6688 total)
	Comparing 11735_DN2b_PU1_peaks_15.txt (6688 total) and ../MACS2_peak_files/ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed (20616 total)
	Comparing ../MACS2_peak_files/ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed (20616 total) and 11735_DN2b_PU1_peaks_15.txt (6688 total)
	Comparing ../MACS2_peak_files/ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed (20616 total) and ../MACS2_peak_files/ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed (20616 total)

11735_DN2b_PU1_peaks_15.txt	../MACS2_peak_files/ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed	Total	Name
	X	18329	../MACS2_peak_files/ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed
X		4398	11735_DN2b_PU1_peaks_15.txt
X	X	2284	11735_DN2b_PU1_peaks_15.txt|../MACS2_peak_files/ATAC-seq_DN3_T_cells_rep1_2.overlapping.bed
```
These output files can now be used for downstream analysis to answer, for instance, if PU.1 binding to promoters are enriched in open and active chromatin or can we identify co-enrichment of other motifs in a specific category etc?


**2. Use deepTools to plot the H3K4me2 and ATAC fragment distribution around PU.1 peaks.**

In this part we will plot the read distribution of PU.1, H3K4me2 and ATAC, around a combined list of DN1 + DN2b PU.1 peaks.

https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html

This is a multistep rocket:

* Start by combining the DN1 (ETP) and DN3 PU.1 peak files:

```bash
module load GCC/7.3.0-2.30 homer/4.10

mergePeaks 11307_DN1_PU1_peaks_15.txt 11735_DN2b_PU1_peaks_15.txt > PU1_peaks.txt
#Turn this into a bed-file.
pos2bed.pl PU1_peaks.txt > PU1_peaks.bed

module purge
```

* Then we will create the files (bigwigs) to be used as the input for the matrix to be visualized. **BONUS:** These bigwig tracks can also be visulized in a genome browser, like IGV or UCSC which is something I **ALWAYS** would recommend. This is done with the bamCoverage function in deepTools but first, load the deepTools module.


```bash
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 deepTools/2.5.4-Python-3.6.6
```
```bash
bamCoverage -h
```
Create a sbatch script for this task.
```bash
#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 4 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 00:20:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'bamCoverage' # name of the job
#SBATCH -o 'bamCoverage%j.out' # stdout log file
#SBATCH -e 'bamCoverage%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.

sample_path="/home/jonun/NAS2/ChIP-ATAC-workshop-2019-course-material"
#genome_size=2652783500 #mm10 genome size non-N bases
module load GCC/6.4.0-2.28  OpenMPI/2.1.2 SAMtools/1.7

#First we must index the bam-file with samtools.
samtools index ${sample_path}/${sample}/${sample}.filt.bam
module purge

#Load the deepTools module
module load GCC/7.3.0-2.30 OpenMPI/3.1.1
module load deepTools/2.5.4-Python-3.6.6


#Make bigwigs.
bamCoverage -b ${sample_path}/${sample}/${sample}.filt.bam  --normalizeUsingRPKM -p 4 --centerReads -o ${sample_path}/${sample}/${sample}.bw

module purge
```
Create a wrapper that loops through the files (to save time we select one of the replicates).
```bash
#!/bin/bash
#A wrapper that loops over the samples and calls batch_bamCoverage.sh for each sample.

samples="11307_DN1_PU1 11735_DN2b_PU1 10340_DN1_H3K4me2 10500_DN3_H3K4me2 ATAC-seq_ETP_T_cells_rep1 ATAC-seq_DN3_T_cells_rep1"

for i in ${samples};do

        echo "sbatch --export=sample=${i} batch_bamCoverage.sh"
        sbatch --export=sample=${i} batch_bamCoverage.sh

done
```
Run this:
```bash
bash bamCoverage_wrapper.sh
```

Great! Now we have created us some bigwig files to compute a matrix from with computeMatrix.
Write a sbatch script for this task.
```bash
#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 1 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 00:30:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'computeMatrix' # name of the job
#SBATCH -o 'computeMatrix%j.out' # stdout log file
#SBATCH -e 'computeMatrix%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.

#Load the deepTools module
module load GCC/7.3.0-2.30 OpenMPI/3.1.1
module load deepTools/2.5.4-Python-3.6.6

sample_path="/projects/fs3/jonun/ChIP-ATAC-workshop-2019-course-material"

#Make matrix.
computeMatrix reference-point -S ${sample_path}/11307_DN1_PU1/11307_DN1_PU1.bw ${sample_path}/11735_DN2b_PU1/11735_DN2b_PU1.bw ${sample_path}/10340_DN1_H3K4me2/10340_DN1_H3K4me2.bw ${sample_path}/10500_DN3_H3K4me2/10500_DN3_H3K4me2.bw ${sample_path}/ATAC-seq_ETP_T_cells_rep1/ATAC-seq_ETP_T_cells_rep1.bw ${sample_path}/ATAC-seq_DN3_T_cells_rep1/ATAC-seq_DN3_T_cells_rep1.bw -R ${sample_path}/merged_peak_files/PU1_peaks.bed -b 1250 -a 1250 --skipZeros --referencePoint center -o ${sample_path}/chip_atac_matrix.matrix.gz

module purge
```
This will create a matrix center in the middle of each peak and extended outwards with 1250 bp in both directions. Read densities will be calculated with a bin size of 10 bp but this can be changed with the `-bs` option.

Now are we finally ready to plot the heatmap and while doing so we will also cluster our data into two clusters.
```bash
#! /bin/bash
#SBATCH -A lsens2018-3-6 # the ID of our Aurora project
#SBATCH -n 1 # how many processor cores to use
#SBATCH -N 1 # how many processors to use (always use 1 here unless you know what you are doing)
#SBATCH -t 00:20:00 # kill the job after ths hh::mm::ss time
#SBATCH -J 'plotHeatmap' # name of the job
#SBATCH -o 'plotHeatmap%j.out' # stdout log file
#SBATCH -e 'plotHeatmap%j.err' # stderr log file
#SBATCH -p dell # which partition to use
# A script that preprocesses the ChIP samples.

#Load the deepTools module
module load GCC/7.3.0-2.30 OpenMPI/3.1.1
module load deepTools/2.5.4-Python-3.6.6

sample_path="/projects/fs3/jonun/ChIP-ATAC-workshop-2019-course-material"
#Plot heatmap
plotHeatmap -m ${sample_path}/chip_atac_matrix.matrix.gz --kmeans 2 --colorMap jet -out ${sample_path}/workshop_chip_atac_heatmap.png -T workshop_chip_atac_heatmap
module purge
```
Behold! Peaks separated between these clusters can furter be targeted for downstream analysis.


###Summary

> **1.** TF and histone peak finding with Homer.

> **2.** Use Homer's mergePeaks to derive unique peak subsets.

> **3.** Peak annotation with annotatePeaks.pl.

> **4.** Motif analysis with Homer and MEME.

> **5.** ATAC-seq peak calling with MACS2.

> **6.** Peak manipulation with bedtools.

> **7.** ChIP/ATAC peak overlap with mergePeaks.

> **8.** Create bigwigs and a basic heatmap with deepTools.

**BONUS:** **FRIP-score calculation with Homer.**

We will start by calculating the reads in peaks. We can do this with `annotatePeaks.pl` if we provide two extra arguments `-raw' and '-d /path-to-tag-directory'

```bash
module load GCC/7.3.0-2.30 homer/4.10

annotatePeaks.pl 11307_DN1_PU1_peaks.txt mm10 -raw -d ~/NAS2/ChIP-ATAC-workshop-2019-course-material/11307_DN1_PU1/11307_DN1_PU1_tagdir > 11307_DN1_PU1_peaks.counts.txt
```
This will create an annotate peak file with a new column (column 21) containing the read counts in all peaks. We can use awk to summarise this:

`awk -F "\t" '{sum += $20} END {print sum}' 11307_DN1_PU1_peaks.counts.txt`
> Reads in peaks = 1877944

The total tags/reads can be found in the tagInfo.txt file in the tag directory. Note that the counts will be halfed if the data is paired-end.

> Total tags/reads = 8848637

> FRIP = 1877944/8848637 = 0.21


There are many more ways to do this if Homer is not availiable. For instance featureCounts from the subread package can be used if providing a bam and a peak-file. Experiment away!
