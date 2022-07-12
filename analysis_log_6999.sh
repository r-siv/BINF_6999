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
#there are 44 sample files total (sample index indicated by number after letter S)
#files starting with 01 to 25 (bull id) are all bull sRNA files to be used
#bulls 01-14 are low fertility
#bulls 15-24 are high fertility
#capitalized letter after bull id represents collection type
#lower case letter after collection type represents technical replicate (may not always be there)
#ex: 01Ba_S13_R1_001.trim.cat -> bull 1, collection B, technical replicate a, sample 13
#ex: 04A_S12_R1_001.trim.cat -> bull 4, collection A, sample 12
#files starting with I are IVF-produced embryos
#files starting with P are parthenogenesis-produced embryos
#DNA files are other extra files run on same flow cell

#extract 3 files from archive
tar -xvf bov_sperm_concat_trimmed.fastq.tar.gz 01Bb_S6_R1_001.trim.cat 01Ba_S13_R1_001.trim.cat 24A_S44_R1_001.trim.cat
#read counts
grep "@" 01Bb_S6_R1_001.trim.cat | wc -l
#there are 3419851 reads in this file
grep "@" 01Ba_S13_R1_001.trim.cat | wc -l
#there are 9835594 reads in this file
grep "@" 24A_S44_R1_001.trim.cat | wc -l
#there are 2306411 reads in this file

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
scp rsivakum@graham.computecanada.ca:~/scratch/6999/fastqc_out/*.html .

#check and load python module
module spider python
module load python

#set up and activate virtual environment for installing multiqc through pip
virtualenv --no-download ENV
source ENV/bin/activate
pip install --no-index --upgrade pip
pip install multiqc

#generate multiqc reports from downsampled fastq files in current directory
multiqc .

#download html results from output directory onto current local machine directory for viewing
#this is done from a new terminal on the local machine
scp rsivakum@graham.computecanada.ca:~/scratch/6999/multiqc_report.html .

#upload unitas perl script (obtained from https://www.smallrnagroup.uni-mainz.de/software.html) to 6999 directory on cluster
#this is done from a new terminal on the local machine
scp unitas_1.7.0.zip rsivakum@graham.computecanada.ca:~/scratch/6999

#unzip unitas script into current 6999 directory and change directories
unzip unitas_1.7.0.zip
cd unitas

#check and load perl module
module spider perl
module load perl

#trim adaptors, collapse and filter low complexity reads (default filter at 75% repeats for a sequence) automatically on 01Ba file
perl unitas_1.7.0.pl -i ~/scratch/6999/01Ba_S13_R1_001.trim.cat.downsampled.fastq -s bos_taurus -trim ?
#####
#NOTE: Fraction of reads containing putative adapter sequence is low!
#Continue nevertheless.
#Clip most abundant sequence (0.85%): TCCCCGACGGGG-(N)x-3'... OK.
#Total processed reads: 10000
#Clipped reads: 126
#10-50 nt: 126
#NOTE: Fraction of clipped reads is low!
#####
#trim adaptors, collapse and filter low complexity reads automatically on 01Bb file
perl unitas_1.7.0.pl -i ~/scratch/6999/01Bb_S6_R1_001.trim.cat.downsampled.fastq -s bos_taurus -trim ?
#####
#NOTE: Fraction of reads containing putative adapter sequence is low!
#Continue nevertheless.
#Clip most abundant sequence (2.15%): TAGATCCGAACT-(N)x-3'... OK.
#Total processed reads: 10000
#Clipped reads: 0
#10-50 nt: 0
#NOTE: Fraction of clipped reads is low!
#####

#note when using -trim ? option, different files are created so unsure if it should be performed (lines 117 & 128)

#count number of annotated piRNA sequences from 01Bb file (#2 is untrimmed file version)
grep ">" UNITAS_10-04-2022_01Bb_S13_R1_001.trim.cat.downsampled.fastq_#2/fasta/unitas.piRcandidates.fas | wc -l
#there are 1600 annotated reads
#count number of unannotated reads from this file
grep ">" UNITAS_10-04-2022_01Bb_S13_R1_001.trim.cat.downsampled.fastq_#2/fasta/unitas.no-annotation.fas | wc -l
#there are 3742 unannotated reads
#count number of annotated piRNA sequences from 01Ba file (#2 is untrimmed file version)
grep ">" UNITAS_10-04-2022_01Ba_S13_R1_001.trim.cat.downsampled.fastq_#2/fasta/unitas.piRcandidates.fas | wc -l
#there are 1021 annotated reads
#count number of unannotated reads from this file
grep ">" UNITAS_10-04-2022_01Ba_S13_R1_001.trim.cat.downsampled.fastq_#2/fasta/unitas.no-annotation.fas | wc -l
#there are 5473 unannotated reads

#download whole folder to display html results from both 01Ba and 01Bb directories onto current local machine directory
#this is done from a new terminal on the local machine
scp -r rsivakum@graham.computecanada.ca:~/scratch/6999/unitas/UNITAS_10-04-2022_01Bb_S6_R1_001.trim.cat.downsampled.fastq_#2 .
scp -r rsivakum@graham.computecanada.ca:~/scratch/6999/unitas/UNITAS_10-04-2022_01Ba_S13_R1_001.trim.cat.downsampled.fastq_#2 .
