library("ggplot2")
library("ggExtra")
library("tidyverse")
library("dplyr")
library("tidyr")
library("ggrepel")
library("directlabels")

##length distributions
filelist = list.files(pattern="*_analysis.txt")
files = lapply(filelist,read.delim)

table = reduce(files, full_join, by="length")
table = table %>%
  mutate(table[,-1]/100000)
table_df = as.data.frame(table)
col_names = list("length","01Ba","01Bb","02A","03A","04A","04Ba","04Bb","05B",
                 "06A","07D","07Da","07Db","08B","09B","10B","11A","12A","13A",
                 "14B","15A","17A","18A","19A","20A","21A","22A","23A","24A",
                 "25A")
names(table_df)=c(col_names)

newtable = gather(table_df,sample,distribution,"01Ba":"25A",factor_key=TRUE)

newtable %>%
  mutate(across(length, as.numeric))

high_fertility = c("15A","17A","18A","19A","20A","21A","22A","23A","24A","25A")

newtable$Fertility = as.factor(ifelse(newtable$sample %in% high_fertility,
                                      "High", "Low"))

ggplot(newtable,aes(x=length,y=distribution, group=sample, color=Fertility,
                    label=sample)) +
  geom_line(size = 0.9) +
  geom_point(size = 1.9 ,alpha = 0.7) +
  removeGridX() +
  scale_x_continuous(breaks=seq(24,32,1)) +
  scale_y_continuous(breaks=seq(0,40,2)) +
  geom_dl(aes(label = sample), method = 
            list(dl.trans(x = x + 0.2),"last.points"),cex = 0.8) +
  ggtitle("Read Length Distribution between Low and High Fertility Bulls") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Reads (per 100,000 reads)") +
  xlab("Read Length")

##pingpong analysis
filelist2 = list.files(pattern="*_pp.txt")
files2 = lapply(filelist2,read.delim)
table2 = reduce(files2, full_join, by="overlap..nt.")
table2 = table2 %>%
  mutate(table2[,-1]/1000)
table2[is.na(table2)] = 0
table_df2 = as.data.frame(table2)
col_names2 = list("overlap","01Ba","01Bb","02A","03A","04Ba","04Bb","05B",
                 "06A","07D","07Da","07Db","08B","09B","10B","11A","12A","13A",
                 "14B","15A","17A","18A","19A","20A","21A","22A","23A","24A",
                 "25A")
names(table_df2)=c(col_names2)

newtable2 = gather(table_df2,sample,read_pairs,"01Ba":"25A",factor_key=TRUE)

newtable2 %>%
  mutate(across(overlap, as.numeric))

high_fertility = c("15A","17A","18A","19A","20A","21A","22A","23A","24A","25A")

newtable2$Fertility = as.factor(ifelse(newtable2$sample %in% high_fertility,
                                      "High", "Low"))

ggplot(newtable2,aes(x=overlap,y=read_pairs,fill=Fertility)) +
  geom_bar(position = "dodge",stat="identity") +
  scale_x_continuous(breaks=seq(1,31,1)) +
  scale_y_continuous(breaks=seq(0,2800,500)) +
  ggtitle("Ping-Pong Signature Comparison between Low and High Fertility Bulls") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Read Pair Counts (per 1000 reads)") +
  xlab("Overlap Position")
##Cluster expression
filelist3 = list.files(pattern="*results.table")
files3 = lapply(filelist3,read.delim)

table3 = reduce(files3, full_join)
table3 = table3[,2:4]

table3 = data.frame(lapply(table3,function(x){
  gsub("Location: ","",x)}))
table3 = data.frame(lapply(table3,function(x){
  gsub("Hits \\(normalized): ","",x)}))

table3$Chromosome = strtrim(table3$Chromosome,5)

table3$Chromosome = as.character(table3$Chromosome)
table3$Chromosome = factor(table3$Chromosome, levels = 
                             unique(table3$Chromosome))
ggplot(data=table3,aes(x=Chromosome,fill=Fertility)) +
  geom_bar() +
  ggtitle("piRNA Cluster Expression Across Chromosomes") +
  removeGridX() +
  scale_y_continuous(breaks=seq(0,500,50)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Normalized Reads") +
  xlab("Chromosome")+
  theme(axis.text.x = element_text(angle = 90))
##TE representation
filelist4 = list.files(pattern="*repeat.table")
files4 = lapply(filelist4,read.delim)

files4[11] = NULL
table4 = reduce(files4, full_join)

table4[table4==0] = NA

table4 = subset(table4,Sense>1)
table4 = subset(table4,Antisense>1)

table4 = table4 %>%
  filter(!if_all(c(Sense,Antisense), is.na)) %>%
  mutate(table4[2]/1000) %>%
  mutate(table4[3]/1000)

ggplot(data=table4,aes(x=TE,y=Sense,fill=Fertility)) +
  geom_bar(position = "dodge",stat="identity") +
  ggtitle("piRNA-TE Representation (Sense)") +
  removeGridX() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Reads (per 1000 reads)") +
  xlab("TE")
ggplot(data=table4,aes(x=TE,y=Antisense,fill=Fertility)) +
  geom_bar(position = "dodge",stat="identity") +
  ggtitle("piRNA-TE Representation (Antisense)") +
  removeGridX() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Reads (per 1000 reads)") +
  xlab("TE")

