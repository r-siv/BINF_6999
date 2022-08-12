library("ggplot2")
library("ggExtra")
library("tidyverse")
library("dplyr")
library("tidyr")
library("rstatix")

#length distributions----

#download files
filelist = list.files(pattern="*_analysis.txt")
files = lapply(filelist,read.delim)

#merge all data tables
table = reduce(files, full_join, by="length")
table_df = as.data.frame(table)

#add sample names to columns
col_names = list("length","02A","04A","07Da","11A","13A","14B","15A","20A",
                 "21A","22A","23A","24A","25A")
names(table_df)=c(col_names)

#divide values in each column by total read count to get proportions
total_reads_per_file = c(2762499,794227,5349077,5514557,3341907,2322465,
                         3976507,2396285,1180511,2134892,2062643,654471,2890194)
table_df[,-1] = sweep(table_df[,-1],2,total_reads_per_file,FUN = "/")
table_df[,-1] = table_df[,-1]*100

#convert to long format
newtable = gather(table_df,sample,distribution,"02A":"25A",factor_key=TRUE)

newtable %>%
  mutate(across(length, as.numeric))

#add fertility as categorical variable
high_fertility = c("15A","17A","18A","19A","20A","21A","22A","23A","24A","25A")
newtable$Fertility = as.factor(ifelse(newtable$sample %in% high_fertility,
                                      "High", "Low"))

#plot read lengths
ggplot(newtable,aes(x=length,y=distribution, group=sample, color=Fertility,
                    label=sample)) +
  geom_line(size = 0.9) +
  geom_point(size = 1.9 ,alpha = 0.7) +
  theme_bw() +
  removeGridX() +
  scale_x_continuous(breaks=seq(24,32,1)) +
  scale_y_continuous(breaks=seq(0,40,2)) +
  ggtitle("Read Length Distribution between Low and High Fertility Bulls") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 20)) +
  ylab("Read Count Percentages (Proportion per sample)") +
  xlab("Read Length") +
  theme(legend.position = "bottom")

#ping-pong analysis----

#download files
filelist2 = list.files(pattern="*_pp.txt")
files2 = lapply(filelist2,read.delim)

#merge all data tables
table2 = reduce(files2, full_join, by="overlap..nt.")

#scale down all values by dividing by 1000 and remove all NA values
table2 = table2 %>%
  mutate(table2[,-1]/1000)
table2[is.na(table2)] = 0
table_df2 = as.data.frame(table2)

#add sample names to columns
col_names2 = list("overlap","02A","04A","07Da","11A","13A","14B","15A","20A",
                  "21A","22A","23A","24A","25A")
names(table_df2)=c(col_names2)

#apply normalization (scale between 0 - 1) to values 
table_df2.scal = apply(table_df2[,-1], MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
table_df2.scal = cbind(table_df2[,1],table_df2.scal)
colnames(table_df2.scal)[1] = "overlap"

#convert to long format
newtable2 = gather(as.data.frame(table_df2.scal),sample,read_pairs,"02A":"25A",factor_key=TRUE)

#remove overlaps over 15 nt in length
newtable2 = newtable2 %>%
  mutate(across(overlap, as.numeric)) %>%
  filter(overlap < 15)

#add fertility as categorical variable
high_fertility = c("15A","17A","18A","19A","20A","21A","22A","23A","24A","25A")
newtable2$Fertility = as.factor(ifelse(newtable2$sample %in% high_fertility,
                                       "High", "Low"))

#plot ping pong overlap distribution
ggplot(newtable2,aes(x=overlap,y=read_pairs,fill=Fertility)) +
  geom_bar(position = "dodge",stat="identity") +
  theme_bw() +
  removeGridX() +
  scale_x_continuous(breaks=seq(1,14,1)) +
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  ggtitle("Ping-Pong Signature Comparison between Low and High Fertility Bulls") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 18)) +
  theme(plot.title = element_text(size=17)) +
  ylab("Normalized Read Pair Counts") +
  xlab("Overlap Length")+
  theme(legend.position = "bottom")

#pingpong z-scores
#generating data frame from Z-score values
zscoresdataframe = data.frame(Z_Score = c(5.18500831671655,22.4375016494233,
                                          30.02752602621,-0.604257960862786,
                                          9.26720930799129,33.7005133929326,
                                          0.198500168812854,2.75038541104538,
                                          15.3455393682021,38.7255160652387,
                                          2.93468257373282,7.05967034064722,
                                          2.84804109723292),
                              Samples = c("02A","04A","07Da","11A","13A","14B",
                                          "15A","20A","21A","22A","23A","24A",
                                          "25A"))

zscoresdataframe$Fertility = as.factor(ifelse(zscoresdataframe$Samples %in% high_fertility,
                                              "High", "Low"))
str(zscoresdataframe)
summary(zscoresdataframe)

#plotting Z-scores
ggplot(zscoresdataframe,aes(x = Samples,y = Z_Score, fill = Fertility)) +
  geom_col() +
  theme_bw() +
  removeGridX() +
  ggtitle("Z-Scores for Ping-Pong Overlaps between Low and High Fertility Samples") +
  scale_y_continuous(breaks=seq(-1,40,6)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 18)) +
  theme(plot.title = element_text(size=17)) +
  ylab("Z Scores") +
  xlab("Samples")+
  theme(legend.position = "bottom")

#normality testing
shapiro.test(zscoresdataframe$Z_Score)

#testing for differences in Z-Scores between fertility states
t.test(zscoresdataframe$Z_Score ~ zscoresdataframe$Fertility, 
       data = zscoresdataframe)

#Cluster expression----

#download files
filelist3 = list.files(pattern="*results.table")
files3 = lapply(filelist3,read.delim)

#merge all data tables and remove first column of clusters
table3 = reduce(files3, full_join)
table3 = table3[,2:4]

#edit column values to remove all prefixes
table3 = data.frame(lapply(table3,function(x){
  gsub("Location: ","",x)}))
table3 = data.frame(lapply(table3,function(x){
  gsub("Hits \\(normalized): ","",x)}))
table3$Chromosome = strtrim(table3$Chromosome,5)

#make chromosome column a factor
table3$Chromosome = as.character(table3$Chromosome)
table3$Chromosome = factor(table3$Chromosome, levels = 
                             unique(table3$Chromosome))

#remove all rows containing chrUn as a chromosome
table3 = table3[!grepl("chrUn", table3$Chromosome),]

#split table by fertility state
table3LHX<-split(table3, table3$Fertility)

#plot cluster expression for high fertility samples
ggplot(data=table3LHX$high,aes(x=Chromosome)) +
  geom_bar(fill = "#F8766D") +
  ggtitle("piRNA Cluster Expression Across Chromosomes (high fertility)") +
  theme_bw() +
  removeGridX() +
  scale_y_continuous(breaks=seq(0,500,25)) +
  scale_x_discrete() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 21)) +
  ylab("Normalized Read Counts") +
  xlab("Chromosome")+
  theme(axis.text.x = element_text(size = 21, angle = 90, 
                                   hjust = 0.95,vjust = 0.2))

#plot cluster expression for low fertility samples
ggplot(data=table3LHX$low,aes(x=Chromosome)) +
  geom_bar(fill = "#00BFC4") +
  ggtitle("piRNA Cluster Expression Across Chromosomes (low fertility)") +
  theme_bw() +
  removeGridX() +
  scale_y_continuous(breaks=seq(0,500,25)) +
  scale_x_discrete() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 21)) +
  ylab("Normalized Read Counts") +
  xlab("Chromosome")+
  theme(axis.text.x = element_text(size = 21, angle = 90, 
                                   hjust = 0.95,vjust = 0.2))

#plot cluster expression for all fertility samples
ggplot(data=table3,aes(x=Chromosome,fill=Fertility)) +
  geom_bar(position = "dodge") +
  ggtitle("piRNA Cluster Expression Across Chromosomes") +
  theme_bw() +
  removeGridX() +
  scale_y_continuous(breaks=seq(0,80,10)) +
  scale_x_discrete() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 21)) +
  ylab("Number of Clusters") +
  xlab("Chromosome")+
  theme(axis.text.x = element_text(size = 21, angle = 90, 
                                   hjust = 0.95,vjust = 0.2))+
  theme(legend.position = "bottom")

#TE representation----

#download files
filelist4 = list.files(pattern="*repeat.table")
files4 = lapply(filelist4,read.delim)

#merge data tables
table4 = reduce(files4, full_join)

#rename zero values as NA
table4[table4==0] = NA

#remove rows if NA values exist in both Sense,Antisense columns
#remove rows if any values are less than 10
table4 = table4 %>%
  filter(!if_all(c(Sense,Antisense), is.na)) %>%
  filter_at(vars(2:3), any_vars(. > 10))

#scale all values by dividing all values by 1000
table4 = table4 %>%
  mutate(table4[2]/1000) %>%
  mutate(table4[3]/1000) %>%
  filter_at(vars(2:3), any_vars(. > 10))

#remove all rows with rRNA targets
table4 = table4[!grepl("rRNA", table4$TE),]

#make dataframe with very low values under opposite fertility state for 
# specific TE targets
zerodataframe = data.frame(TE = c("LINE/RTE-BovB","LINE/RTE-BovB",
                                  "LTR/ERVL-MaLR","LTR/ERVL-MaLR"),
                           Sense = c(0.01,0.01,0.01,0.01),
                           Antisense = c(0.01,0.01,0.01,0.01),
                           Fertility = c("high","high","high","high"))

#merge above data frame with main table to remove any wide bars (occurs when data
# only exists under 1 fertility state)
table4 = rbind(table4,zerodataframe)

#plot TE targets for sense direction
ggplot(data=table4,aes(x=TE,y=Sense,fill=Fertility)) +
  geom_bar(position = "dodge",stat="identity") +
  coord_flip() +
  ggtitle("piRNA-TE Representation (Sense)") +
  theme_bw() +
  removeGridX() +
  removeGridY() +
  scale_y_continuous(breaks=seq(0,80,10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 18)) +
  ylab("Read Pair Counts (thousands)") +
  xlab("Predicted TE targets")+
  theme(legend.position = "bottom")

#plot TE targets for antisense direction
ggplot(data=table4,aes(x=TE,y=Antisense,fill=Fertility)) +
  geom_bar(position = "dodge",stat="identity") +
  coord_flip() +
  ggtitle("piRNA-TE Representation (Antisense)") +
  theme_bw() +
  removeGridX() +
  removeGridY() +
  scale_y_continuous(breaks=seq(0,80,10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 18)) +
  ylab("Read Pair Counts (thousands)") +
  xlab("Predicted TE targets")+
  theme(legend.position = "bottom")