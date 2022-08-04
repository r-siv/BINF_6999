#BINF*6999 Project 
#BASH script for Small RNA-Seq analysis in bovine sperm
#This analysis log contains all the steps taken in the pipeline

#upload fastq.tar.gz file to scratch directory on cluster with SCP protocol
scp bov_sperm_concat_trimmed.fastq.tar.gz rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/

#log in to cluster and change directories to projects (containing fastq.tar.gz)
ssh rsivakum@graham.computecanda.ca
cd projects/def-jlamarre/rsivakum/

#if files were placed in scratch directory and are to be moved into projects directory, follow instructions at below link: 
#https://docs.alliancecan.ca/wiki/Frequently_Asked_Questions#Moving_files_across_the_project.2C_scratch_and_home_filesystems

#make new directory and move fastq.tar.gz into that
mkdir 6999
mv bov_sperm_concat_trimmed.fastq.tar.gz 6999/

#change directories and list all files within tar archive
cd 6999/
tar -tf bov_sperm_concat_trimmed.fastq.tar.gz

#count files
tar -tf bov_sperm_concat_trimmed.fastq.tar.gz | wc -l
#there are 44 sample files total (sample index indicated by number after letter S)
#files starting with 01 to 25 (bull id) are all bull sRNA files to be used (all in fastq format)
#bulls 01-14 are low fertility
#bulls 15-24 are high fertility
#capitalized letter after bull id represents collection type
#lower case letter after collection type represents technical replicate (may not always be there)
#ex: 01Ba_S13_R1_001.trim.cat -> bull 1, collection B, technical replicate a, sample 13
#ex: 04A_S12_R1_001.trim.cat -> bull 4, collection A, sample 12
#files starting with I are IVF-produced embryos
#files starting with P are parthenogenesis-produced embryos
#DNA files are other extra files run on same flow cell

#extract 3 files from archive for testing
tar -xvf bov_sperm_concat_trimmed.fastq.tar.gz 01Ba_S13_R1_001.trim.cat 01Bb_S6_R1_001.trim.cat 24A_S44_R1_001.trim.cat
#make new testing directory and move the newly extracted files into that
mkdir testing_data
mv *.cat testing_data

#check read counts for all test files
for file in testing_data/*; 
do 
    echo "$file"; grep "@" $file | wc -l; 
done
#there are 9835594 reads in 01Ba_S13_R1_001.trim.cat
#there are 3419851 reads in 01Bb_S6_R1_001.trim.cat
#there are 2306411 reads in 24A_S44_R1_001.trim.cat

#check and load seqtk module for randomly downsampling 10k reads
module spider seqtk
module load seqtk

#setup interactive node (to be done for less computationally intensive steps)
salloc --time=4:0:0 --ntasks=2 --account=def-jlamarre

#randomly downsample all files
#NOTE: Skip this step when working with complete data (downsampling done only for testing pipeline)
for file in testing_data/*; 
do 
    seqtk sample $file 10000 > $file.down.fastq; 
done

#check for correct read numbers
for file in testing_data/*.fastq; 
do 
    echo "$file"; grep "@" $file | wc -l; 
done
#10000 reads found in all newly downsampled files

#check and load fastqc module
module spider fastqc
module load fastqc

#make new fastqc directory to perform preliminary fastqc statistics on all downsampled files. Move output into new directory
mkdir fastqc_out
for file in testing_data/*.fastq; 
do 
    fastqc $file -o fastqc_out/; 
done

#download html results from output directory onto current local machine directory for viewing
#this is done from a new terminal on the local machine
scp rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/fastqc_out/*.html .

#check and load python module
module spider python
module load python

#set up and activate virtual environment for installing multiqc through pip
virtualenv --no-download ENV
source ENV/bin/activate #run this line everytime multiqc is to be performed (similar to module load) and exit virtual environment using deactivate
pip install --no-index --upgrade pip
pip install multiqc

#generate multiqc reports from fastqc files into current directory
multiqc fastqc_out/*
deactivate

#download html results from output directory onto current local machine directory for viewing
#this is done from a new terminal on the local machine
scp rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/multiqc_report.html .

#make new directory called scripts for storing scripts to be used downstream
mkdir scripts

#upload unitas, length filter, TBr2 basic, sRNAmapper, reallocate, proTrac, pingpong and repeatmasker files (obtained from https://www.smallrnagroup.uni-mainz.de/software.html) to scripts directory on cluster
#this is done from a new terminal on the local machine
scp unitas_1.7.0.zip TBr2_length-filter.pl TBr2_basic-analyses.pl sRNAmapper.pl reallocate.pl proTRAC_2.4.2.pl TBr2_pingpong.pl RMvsMAP.pl rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/scripts

#unzip unitas script
unzip scripts/unitas_1.7.0.zip

#check and load perl module
module spider perl
module load perl

#copy fastq files into temporary directory for easier parsing with Unitas
mkdir for_unitas
cp testing_data/*.fastq for_unitas

#record all user input and printed output to terminal using script (end recording with CTRL-D or exit)
script unitaslog.txt

#annotate, collapse and filter low complexity reads (default filter at 75% repeats for a sequence) automatically on all files (4 threads used)
perl unitas/unitas_1.7.0.pl -i for_unitas -s bos_taurus -threads 4
#data is already trimmed so no need for -trim

#for 01Ba_S13_R1_001.trim.cat.down.fastq 9893 sequences are ok, 107 sequences are dusty (low-complexity)
#  Monomer:             17
#  Biased for 2 bases:  61
#  2 nt motif enriched: 3
#  3 nt motif enriched: 8
#  4 nt motif enriched: 1
#  5 nt motif enriched: 3

#for 01Bb_S6_R1_001.trim.cat.down.fastq 9930 sequences are ok, 70 sequences are dusty (low-complexity)
#  Monomer:             8
#  Biased for 2 bases:  50
#  2 nt motif enriched: 1
#  3 nt motif enriched: 5
#  4 nt motif enriched: 2
#  5 nt motif enriched: 1

#for 24A_S44_R1_001.trim.cat.down.fastq 9789 sequences are ok, 211 sequences are dusty (low-complexity)
#  Monomer:             43
#  Biased for 2 bases:  113
#  2 nt motif enriched: 1
#  3 nt motif enriched: 7
#  4 nt motif enriched: 1
#  5 nt motif enriched: 4

#count number of annotated piRNA sequences from each file
for group in *\#1; 
do 
    echo $group; grep ">" $group/fasta/unitas.piRcandidates.fas | wc -l; 
done
#there are 1021 piRNA reads in 01Ba_S13_R1_001.trim.cat
#there are 1600 piRNA reads in 01Bb_S6_R1_001.trim.cat
#there are 1806 piRNA reads in 24A_S44_R1_001.trim.cat

#count number of unannotated reads from each file (should contain known and unknown piRNA reads)
for group in *\#1; 
do 
    echo $group; grep ">" $group/fasta/unitas.no-annotation.fas | wc -l; 
done
#there are 5473 unannotated reads in 01Ba_S13_R1_001.trim.cat
#there are 3742 unannotated reads in 01Bb_S6_R1_001.trim.cat
#there are 4560 unannotated reads in 24A_S44_R1_001.trim.cat

#download all html files from one Unitas directory to display results onto current local machine directory (UNITAS_12-07-2022_24A_S44_R1_001.trim.cat.down.fastq_#1 was used as an example)
#this is done from a new terminal on the local machine
#scp should not be used for transferring files from multiple remote directories to multiple local directories
scp -r rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/UNITAS_12-07-2022_24A_S44_R1_001.trim.cat.down.fastq_#1/html .

#copy unannotated files (will contain piRNA reads) from Unitas directories into testing directory
for file in *\#1;
do
    cp $file/fasta/unitas.no-annotation.fas ./testing_data/$file.fas;
done

#make directory for filtered files
mkdir len_filtered

#record terminal output
script length_filter_log.txt

#run piRNA length filter on unannotated files and move them to filtered directory
for file in testing_data/U*;
do
    perl scripts/TBr2_length-filter.pl -i $file -o len_filtered/$file"_len" -min 24 -max 32;
done
#for 01Ba_S13_R1_001.trim.cat.down.fastq_#1.fas
#Found 1560 sequences < 24 nt.
#Found 1124 sequences > 32 nt.
#Found 2789 sequences from 24-32 nt.

#for 01Bb_S6_R1_001.trim.cat.down.fastq_#1.fas
#Found 936 sequences < 24 nt.
#Found 603 sequences > 32 nt.
#Found 2203 sequences from 24-32 nt.

#for 24A_S44_R1_001.trim.cat.down.fastq_#1.fas
#Found 1670 sequences < 24 nt.
#Found 204 sequences > 32 nt.
#Found 2686 sequences from 24-32 nt.

#rename new files with suffix _len and shortened name (remove first 18 characters including "UNITAS", underscores and date)
for file in len_filtered*;
do
    mv "$file" "${file:18}";
done

#check and load fastx toolkit module
module spider fastx-toolkit
module load fastx-toolkit

#convert collapsed _len files into fasta format (uncollapsed) for obtaining length and nucleotide distributions
for file in len_filtered/*_len;
do 
    fasta_formatter -t -i $file -o $file"_tab";
done

for file in len_filtered/*_tab;
do
    awk '{c=$1; while (c-->0) print}' $file > ${file%.*}""_expand.txt;
done

for file in len_filtered/*_expand.txt;
do
    awk '{print ">"$1"\n"$2}' $file > ${file%.*}".fa";
done

#basic analysis script requires uncollapsed fasta/fastq files for obtaining length and nucleotide distributions
for file in len_filtered/*_expand.fa;
do
    perl scripts/TBr2_basic-analyses.pl -i $file -o ${file%.*}""_len_expand_analysis.txt"";
done

#download all .txt files from len_filtered directory for plotting onto current local machine directory
#this is done from a new terminal on the local machine
scp -r rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/len_filtered/*.txt .

#upload twoBitToFa tool (from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/), Bos taurus genome, repeatmasker annotation and ENSEMBL geneset files to main directory on cluster
#this is done from a new terminal on the local machine
scp twoBitToFa bosTau9.2bit bosTau9.fa.out.gz Bos_taurus.ARS-UCD1.2.106.chr.gtf.gz rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999

#unzip all the new Bos taurus files
for file in *.gz;
do
    gunzip $file;
done

#unpack twoBitToFa tool and extract genome from 2bit file
chmod 744 twoBitToFa
./twoBitToFa bosTau9.2bit bosTau9.fa

#append "chr" to the start of all numbers in the first column for the gtf file to match repeatmasker and genome formats
for file in *.gtf;
do
    awk 'OFS="\t" {if (NR > 5) $1="chr"$1}' $file;
done

#make map file directory
mkdir maps

#run sRNAmapper on len files (best alignments in terms of mismatch counts written for each sequence)
#if run on interactive node, follow commands below:

for file in len_filtered/*_len;
do
    perl scripts/sRNAmapper.pl -i $file -g bosTau9.fa -a best -o maps/$file.map;
done

#mapping should be submitted as an array job through SLURM scheduler using SBATCH
touch job.sh
nano job.sh
#the contents of this shell script should be as follows:

#!/bin/bash
#SBATCH --array=0-29
#SBATCH --job-name=" mapping"
#SBATCH --time=03:00:00
#SBATCH --mem=1000
#SBATCH --account=def-jlamarre
files=(len_filtered/*_len)
file=${files[$SLURM_ARRAY_TASK_ID]}

echo $file
output=$(basename ${file}).map

module load perl
perl scripts/sRNAmapper.pl -i $file -g bosTau9.fa -a best -o maps/$output
sleep 30

#submit file as array job
sbatch job.sh

#when all computations are complete, check SLURM output files if any errors occurred
grep -r  "No such" *.out
#increase allocated memory (#SBATCH --mem=) if oom-kill events occur in slurm output

#run ping-pong check on non-weighted map files and rename new files with suffix _pp.txt (remove map suffix)
#if run on interactive node, follow commands below:

for file in maps/*.map;
do
    perl scripts/TBr2_pingpong.pl -i $file -o ${file%.*}"_pp.txt";
done

#ping-pong checks should be submitted as an array job through SLURM scheduler using SBATCH
touch job4.sh
nano job4.sh
#the contents of this shell script should be as follows:

#!/bin/bash
#SBATCH --array=0-29
#SBATCH --job-name=" pingpong"
#SBATCH --time=01:00:00
#SBATCH --mem=5000
#SBATCH --account=def-jlamarre
files=(maps/*.map)
file=${files[$SLURM_ARRAY_TASK_ID]}

echo $file

module load perl
perl scripts/TBr2_pingpong.pl -i $file -o ${file%.*}"_pp.txt"
sleep 30

#submit file as array job
sbatch job4.sh

#record terminal output
script pingpong_log.txt

#Analyzing map file maps/01Ba_S13_R1_001.trim.cat.down.fastq_#1.fas_len.map
#This is a standard map file: +hits: 6401 / -hits: 6056

#Analyzing map file maps/01Bb_S6_R1_001.trim.cat.down.fastq_#1.fas_len.map
#This is a standard map file: +hits: 1790 / -hits: 2074

#Analyzing map file maps/24A_S44_R1_001.trim.cat.down.fastq_#1.fas_len.map
#This is a standard map file: +hits: 19447 / -hits: 19061

#when all computations are complete, check SLURM output files if any errors occurred
grep "error" *.out

#make ping-pong directory
mkdir pingpong

#move ping-pong files to pingpong directory
mv *_pp.txt pingpong/

#download all .txt files from pingpong directory for plotting onto current local machine directory
#this is done from a new terminal on the local machine
scp -r rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/pingpong/* .

#setup interactive node for repeatmasker
salloc --time=2:30:0 --ntasks=2 --account=def-jlamarre

#run repeatmasker annotation on non-weighted map files
for file in maps/*.map;
do
    perl scripts/RMvsMAP.pl bosTau9.fa.out $file
    name = $(basename ${file})
    mv TE-fasta TE-fasta"_${name:0:4}"
    mv repeat_counts.table "${name:0:4}"repeat_counts
done

#download all repeatmasker output files from current directory for plotting onto current local machine directory
#this is done from a new terminal on the local machine
scp -r rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/TE-fasta* .
scp -r rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/*.table .

#run reallocation on map files (perimeter of 5000, resolution of 1000, bell shape function and 0 to reject loci with no allocated reads)
#if run on interactive node, follow commands below:

for file in maps/*;
do
    perl scripts/reallocate.pl -i $file 5000 1000 b 0;
done

#reallocation should be submitted as an array job through SLURM scheduler using SBATCH
touch job2.sh
nano job2.sh
#the contents of this shell script should be as follows:

#!/bin/bash
#SBATCH --array=0-29
#SBATCH --job-name=" reallocating"
#SBATCH --time=01:00:00
#SBATCH --mem=1000
#SBATCH --account=def-jlamarre
files=(analysis/maps/*)
file=${files[$SLURM_ARRAY_TASK_ID]}

echo $file

module load perl
perl scripts/reallocate.pl -i $file 5000 1000 b 0
sleep 30

#submit file as array job
sbatch job2.sh

#move reallocated files to maps directory
mv *.weighted-5000-1000-b-0 maps/

#run proTrac clustering on new weighted files using maximum piRNA length of 32 and output mapped reads as tab-delimited table
#if run on interactive node, follow commands below:

for file in maps/*.weighted-5000-1000-b-0;
do
    perl scripts/proTRAC_2.4.2.pl -map $file -genome bosTau9.fa -repeatmasker bosTau9.fa.out -geneset Bos_taurus.ARS-UCD1.2.106.chr.gtf -pimax 32 -pti;
done

#clustering should be submitted as an array job through SLURM scheduler using SBATCH
touch job3.sh
nano job3.sh
#the contents of this shell script should be as follows:

#!/bin/bash
#SBATCH --array=0-29
#SBATCH --job-name=" clustering"
#SBATCH --time=03:00:00
#SBATCH --mem=5000
#SBATCH --account=def-jlamarre
files=(maps/*.weighted-5000-1000-b-0)
file=${files[$SLURM_ARRAY_TASK_ID]}

echo $file

module load perl
perl scripts/proTRAC_2.4.2.pl -map $file -genome bosTau9.fa -repeatmasker bosTau9.fa.out -geneset Bos_taurus.ARS-UCD1.2.106.chr.gtf -pimax 32 -pti
sleep 30

#submit file as array job
sbatch job3.sh

#record terminal output
script protrac_log.txt

#minimum size of piRNA cluster set to 1000 bp (default)
#genome size (without gaps): 2715825630 bp
#gaps (N/X/-): 939 bp
#Number of Chromosomes/Scaffolds: 2211

#for 01Ba_S13_R1_001.trim.cat.down.fastq_#1.fas_len.map.weighted-5000-1000-b-0
#mapped reads: 1630
#non-identical sequences: 1388
#Genomic hits: 10270 (+:5234 -:5036)
#Total size of 69 predicted piRNA clusters: 866403 bp (0.032%)
#Non identical sequences that can be assigned to clusters: 756 (55.022%)
#Sequence reads that can be assigned to clusters: 849 (52.537%)
#Total repeat-masked bases in clusters: 388451 (44.83%)
#Total repeat-masked bases in genome: 1348220998 (49.64%)

#for 01Bb_S6_R1_001.trim.cat.down.fastq_#1.fas_len.map.weighted-5000-1000-b-0
#mapped reads: 1843
#non-identical sequences: 1707
#Genomic hits: 3337 (+:1507 -:1830)
#Total size of 68 predicted piRNA clusters: 994670 bp (0.037%)
#Non identical sequences that can be assigned to clusters: 1261 (73.959%)
#Sequence reads that can be assigned to clusters: 1349 (73.275%)
#Total repeat-masked bases in clusters: 349019 (35.09%)
#Total repeat-masked bases in genome: 1348220998 (49.64%)

#for 24A_S44_R1_001.trim.cat.down.fastq_#1.fas_len.map.weighted-5000-1000-b-0
#mapped reads: 2437
#non-identical sequences: 2253
#Genomic hits: 7601 (+:3879 -:3722)
#Total size of 72 predicted piRNA clusters: 1315428 bp (0.048%)
#Non identical sequences that can be assigned to clusters: 1292 (57.346%)
#Sequence reads that can be assigned to clusters: 1394 (57.201%)
#Total repeat-masked bases in clusters: 470769 (35.79%)
#Total repeat-masked bases in genome: 1348220998 (49.64%)

#make clustering directory and move proTrac files to clustering directory
mkdir clustering
mv pro* clustering/

#download all protrac files from clustering directory for plotting onto current local machine directory
#this is done from a new terminal on the local machine
scp -r rsivakum@graham.computecanada.ca:~/projects/def-jlamarre/rsivakum/6999/clustering/* .

#all following steps are done on new terminal from local machine for preparing files for data visualization
#assuming all downloaded files are in same folder

#for plotting read length distributions, remove first 5 lines, everything after first 10 lines and remove 3rd column
for file in *_analysis.txt; 
do
    sed -i 1,5d $file
    head -n 10 "$file" > "$file.10" && mv "$file.10" "$file"
    cut -f1,2 "$file" > "$file.cut" && mv "$file.cut" "$file"
done

#for plotting ping-pong results, remove first 2 lines and everything after first 32 lines
for file in *_pp.txt;
do 
    sed -i 1,2d $file
    head -n 32 "$file" > "$file.32" && mv "$file.32" "$file"
done

#for plotting cluster expressions, remove file name prefix, first 72 lines and last 4 lines. Also extract columns 1, 2 and 6
for file in *.table; 
do
    mv $file ${file#proTRAC_}
    sed -i 1,72d $file
    head -n -4 $file > "$file.-4" && mv "$file.-4" "$file"
    cut -f1,2,6 "$file" > "$file.cut" && mv "$file.cut" "$file"
done

#manually move high fertility files (15 - 25) to temporary folder and add fertility labels
mkdir temporary
for file in temporary*; 
do
    sed -i "s/$/\t"high"/" $file;
done

#add low fertility labels to low fertility files in current directory (01 - 14)
cd ..
for file in *; 
do
    sed -i "s/$/\t"low"/" $file;
done

#add headers to each column for both sets of files (current and temporary directory)
for file in *; 
do
    echo -e 'Cluster\tChromosome\tNormalized_Hits\tFertility' | cat - $file > temp && mv temp $file;
done

for file in temporary*; 
do
    echo -e 'Cluster\tChromosome\tNormalized_Hits\tFertility' | cat - $file > temp && mv temp $file;
done

#for plotting TE targets remove lines 2 - 14072 and last 3 lines
for file in *repeat.table; 
do
    sed -i 2,14072d $file
    head -n -3 $file > "$file.-3" && mv "$file.-3" "$file"
done