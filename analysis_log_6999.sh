#Raam Sivakumar
#analysis log for all the steps in the pipeline

#upload fastq.tar.gz file to scratch directory on cluster with SCP
scp bov_sperm_concat_trimmed.fastq.tar.gz rsivakum@graham.computecanada.ca:~/scratch

#log in to cluster and change directories to scratch (containing fastq.tar.gz)
ssh rsivakum@graham.computecanda.ca
cd scratch/

#make new directory and move fastq.tar.gz into that
mkdir 6999
mv bov_sperm_concat_trimmed.fastq.tar.gz 6999/

#list all files within tar archive
tar -tf bov_sperm_concat_trimmed.fastq.tar.gz

#count files
tar -tf bov_sperm_concat_trimmed.fastq.tar.gz | wc -l
#there are 44 files

#extract first 2 files from archive
tar -zxvf bov_sperm_concat_trimmed.fastq.tar.gz 01Bb_S6_R1_001.trim.cat
tar -zxvf bov_sperm_concat_trimmed.fastq.tar.gz 01Ba_S13_R1_001.trim.cat

#read counts
grep "@" 01Bb_S6_R1_001.trim.cat | wc -l
#there are 3419851 reads in this file
grep "@" 01Ba_S13_R1_001.trim.cat | wc -l
#there are 9835594 reads in this file

#check and load seqtk module for randomly downsampling 10k reads
module spider seqtk
module load seqtk

#downsampling both files and checking for downsampling
seqtk sample 01Ba_S13_R1_001.trim.cat 10000 > 01Ba_S13_R1_001.trim.cat.downsampled.fastq
grep "@" 01Ba_S13_R1_001.trim.cat.downsampled.fastq | wc -l
#10000 reads found
seqtk sample 01Bb_S6_R1_001.trim.cat 10000 > 01Bb_S6_R1_001.trim.cat.downsampled.fastq
grep "@" 01Bb_S6_R1_001.trim.cat.downsampled.fastq | wc -l
#10000 reads found

#check and load fastqc module
module spider fastqc
module load fastqc

#perform preliminary statistics with fastqc and move output into new directory
mkdir fastqc_out
fastqc 01Ba_S13_R1_001.trim.cat.downsampled.fastq -o fastqc_out/
#####
#Picked up JAVA_TOOL_OPTIONS: -Xmx2g
#Started analysis of 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 10% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 20% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 30% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 40% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 50% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 60% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 70% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 80% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 90% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Approx 100% complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#Analysis complete for 01Ba_S13_R1_001.trim.cat.downsampled.fastq
#####
fastqc 01Bb_S6_R1_001.trim.cat.downsampled.fastq -o fastqc_out/
###
#Picked up JAVA_TOOL_OPTIONS: -Xmx2g
#Started analysis of 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 10% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 20% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 30% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 40% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 50% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 60% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 70% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 80% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 90% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Approx 100% complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
#Analysis complete for 01Bb_S6_R1_001.trim.cat.downsampled.fastq
###

#download html results from output directory onto current local machine directory for viewing
#this is done from a new terminal on the local machine
scp rsivakum@graham.computecanada.ca:/home/rsivakum/scratch/6999/fastqc_out/*.html .