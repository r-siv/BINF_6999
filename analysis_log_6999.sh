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

#change directories and list all files within tar archive
cd 6999/
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

#randomly downsample all files and check for correct read numbers
for file in testing_data/*; 
do 
    seqtk sample $file 10000 > $file.down.fastq; 
done
for file in testing_data/*.fastq; 
do 
    echo "$file"; grep "@" $file | wc -l; 
done
#10000 reads found in all newly downsampled files

#check and load fastqc module
module spider fastqc
module load fastqc

#change directories and make new fastqc directory to perform preliminary fastqc statistics on all downsampled files. Move output into new directory
cd ..
mkdir fastqc_out
for file in testing_data/*.fastq; 
do 
    fastqc $file -o fastqc_out/; 
done

#download html results from output directory onto current local machine directory for viewing
#this is done from a new terminal on the local machine
scp rsivakum@graham.computecanada.ca:~/scratch/6999/fastqc_out/*.html .

#check and load python module
module spider python
module load python

#set up and activate virtual environment for installing multiqc through pip
virtualenv --no-download ENV
source ENV/bin/activate #run this line everytime multiqc is to be performed (similar to module load) and leave virtual environment using deactivate
pip install --no-index --upgrade pip
pip install multiqc

#generate multiqc reports from fastqc files into current directory
multiqc fastqc_out/*
deactivate

#download html results from output directory onto current local machine directory for viewing
#this is done from a new terminal on the local machine
scp rsivakum@graham.computecanada.ca:~/scratch/6999/multiqc_report.html .

#make new directory called scripts for storing downstream scripts to be used
mkdir scripts

#upload unitas, length filter, sRNAmapper, reallocate, proTrac, pingpong and repeatmasker files (obtained from https://www.smallrnagroup.uni-mainz.de/software.html) to scripts directory on cluster
#this is done from a new terminal on the local machine
scp unitas_1.7.0.zip TBr2_length-filter.pl sRNAmapper.pl reallocate.pl proTRAC_2.4.2.pl TBr2_pingpong.pl RMvsMAP.pl rsivakum@graham.computecanada.ca:~/scratch/6999/scripts

#unzip unitas script
unzip scripts/unitas_1.7.0.zip

#check and load perl module
module spider perl
module load perl

#collapse and filter low complexity reads (default filter at 75% repeats for a sequence) automatically on all files
for file in testing_data/*.fastq; 
do 
    perl unitas/unitas_1.7.0.pl -i $file -s bos_taurus; 
done
#note when using -trim ? option, different files are created so unsure if it should be performed (lines 117 & 128)

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

#this still could be simplified

#count number of unannotated reads from each file
for group in *\#1; 
do 
    echo $group; grep ">" $group/fasta/unitas.no-annotation.fas | wc -l; 
done
#there are 5473 unannotated reads in 01Ba_S13_R1_001.trim.cat
#there are 3742 unannotated reads in 01Bb_S6_R1_001.trim.cat
#there are 4560 unannotated reads in 24A_S44_R1_001.trim.cat

#download all Unitas folders to display html results onto current local machine directory
#this is done from a new terminal on the local machine
scp -r rsivakum@graham.computecanada.ca:~/scratch/6999/*\#1 .

#copy unannotated files (will contain piRNA reads) from Unitas directories into testing directory
for file in *\#1;
do
    cp $file/fasta/unitas.no-annotation.fas ./testing_data/$file.fas;
done

#run piRNA length filter on unannotated files and rename new files with suffix _len
for file in testing_data/U*;
do
    perl scripts/TBr2_length-filter.pl -i $file -o $file""_len"" -min 24 -max 32;
done
#for UNITAS_12-07-2022_01Ba_S13_R1_001.trim.cat.down.fastq_#1.fas
#Found 1560 sequences < 24 nt.
#Found 1124 sequences > 32 nt.
#Found 2789 sequences from 24-32 nt.

#for UNITAS_12-07-2022_01Bb_S6_R1_001.trim.cat.down.fastq_#1.fas
#Found 936 sequences < 24 nt.
#Found 603 sequences > 32 nt.
#Found 2203 sequences from 24-32 nt.

#for UNITAS_12-07-2022_24A_S44_R1_001.trim.cat.down.fastq_#1.fas
#Found 1670 sequences < 24 nt.
#Found 204 sequences > 32 nt.
#Found 2686 sequences from 24-32 nt.

#upload Bos taurus genome, repeatmasker annotation and ENSEMBL geneset files to main directory on cluster
#this is done from a new terminal on the local machine
scp GCF_002263795.2_ARS-UCD1.3_genomic.fna.gz bosTau7.fa.out.gz Bos_taurus.ARS-UCD1.2.106.chr.gtf.gz rsivakum@graham.computecanada.ca:~/scratch/6999/

#unzip all the new Bos taurus files
for file in *.gz;
do
    gunzip $file;
done

#setup interactive node
salloc --time=4:0:0 --ntasks=2 --account=def-jlamarre

#run sRNAmapper on len files (best alignments in terms of mismatch counts written for each sequence)
for file in testing_data/*_len;
do
    perl scripts/sRNAmapper.pl -input $file -genome GCF_002263795.2_ARS-UCD1.3_genomic.fna -alignments best;
done

#run reallocation on new map files (perimeter of 5000, resolution of 1000, bell shape function and 0 to reject loci with no allocated reads)
for file in testing_data/*.map;
do
    perl scripts/reallocate.pl $file 5000 1000 b 0;
done

#run proTrac clustering on new weighted files using using maximum piRNA length of 32 and output mapped reads as tab-delimited table
for file in testing_data/*.weighted-5000-1000-b-0;
do
    perl scripts/proTRAC_2.1.2.pl -map $file -genome GCF_002263795.2_ARS-UCD1.3_genomic.fna -repeatmasker bosTau7.fa.out -geneset Bos_taurus.ARS-UCD1.2.106.chr.gtf -pimax 32 -pti;
done
#minimum size of piRNA cluster set to 1000 bp (default)
#genome size (without gaps): 2711181669 bp
#gaps (N/X/-): 730 bp
#Number of Chromosomes/Scaffolds: 1957

#for UNITAS_12-07-2022_01Ba_S13_R1_001.trim.cat.down.fastq_#1.fas_len.map.weighted-5000-1000-b-0
#mapped reads: 1616
#non-identical sequences: 1374
#Genomic hits: 10258 (+:5250 -:5008)
#Total size of 69 predicted piRNA clusters: 866403 bp (0.032%)
#Non identical sequences that can be assigned to clusters: 756 (55.022%)
#Sequence reads that can be assigned to clusters: 849 (52.537%)
#Total repeat-masked bases in clusters: 0 (0%)
#Total repeat-masked bases in genome: 1394750942 (51.44%)

#for UNITAS_12-07-2022_01Bb_S6_R1_001.trim.cat.down.fastq_#1.fas_len.map.weighted-5000-1000-b-0
#mapped reads: 1841
#non-identical sequences: 1705
#Genomic hits: 3331 (+:1507 -:1824)
#Total size of 68 predicted piRNA clusters: 994670 bp (0.037%)
#Non identical sequences that can be assigned to clusters: 1261 (73.959%)
#Sequence reads that can be assigned to clusters: 1349 (73.275%)
#Total repeat-masked bases in clusters: 0 (0%)
#Total repeat-masked bases in genome: 1394750942 (51.44%)

#for UNITAS_12-07-2022_24A_S44_R1_001.trim.cat.down.fastq_#1.fas_len.map.weighted-5000-1000-b-0
#mapped reads: 2435
#non-identical sequences: 2251
#Genomic hits: 7599 (+:3879 -:3720)
#Total size of 72 predicted piRNA clusters: 1315428 bp (0.049%)
#Non identical sequences that can be assigned to clusters: 1292 (57.397%)
#Sequence reads that can be assigned to clusters: 1394 (57.248%)
#Total repeat-masked bases in clusters: 0 (0%)
#Total repeat-masked bases in genome: 1394750942 (51.44%)

#run ping-pong check on non-weighted map files and rename new files with suffix _pp.txt (remove map suffix)
for file in testing_data/*.map;
do
    perl scripts/TBr2_pingpong.pl -i $file -o  ${file%.*}"_pp.txt";
done

#run repeatmasker annotation on non-weighted map files
for file in testing_data/*.map;
do
    perl scripts/RMvsMAP.pl scripts/bosTau7.fa.out $file;
done