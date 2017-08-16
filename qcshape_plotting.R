#SHAPE-like data plotting
#Operates on output of qcshape

#The simple bare necessities
library(ggplot2)
library(readr)
library(reshape2)
library(viridis)
library(scales)

#Import data
biotypes <- read_delim("biotypes.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

stopbases <- read_delim("stopbases.txt",
                        "\t", escape_double = FALSE, trim_ws = TRUE)
stopbases_perc <- read_delim("stopbases_perc.txt",
                             "\t", escape_double = FALSE, trim_ws = TRUE)

rpkmcorr <- read_delim("rpkmcorr.txt",
                       "\t", escape_double = FALSE, col_names = FALSE,
                       trim_ws = TRUE)
colnames(rpkmcorr) <- c("Sample1","Sample2","r","p")

rtcorr <- read_delim("rtcorr.txt",
                     "\t", escape_double = FALSE, col_names = FALSE,
                     trim_ws = TRUE)
colnames(rtcorr) <- c("Sample1","Sample2","r","p")

reads_mapped <- read_delim("reads_mapped.txt",
                           "\t", escape_double = FALSE, trim_ws = TRUE)

reads_mapped_long <- melt(reads_mapped)

shaperegions <- read_delim("shaperegions.txt",
                           "\t", escape_double = FALSE, col_names = FALSE,
                           col_types = cols(X7 = col_double()),
                           na = "NULL", trim_ws = TRUE)

colnames(shaperegions) <- c("Transcript","Feature","Start","Stop","Strand","LocalPos","Shape")

enrichbases <- read_delim("enrichbases.txt",
                          "\t", escape_double = FALSE, trim_ws = TRUE)
enrichbases_perc <- read_delim("enrichbases_perc.txt",
                               "\t", escape_double = FALSE, trim_ws = TRUE)

enrichbases_list <- read_delim("enrichbases_list.txt",
                               "\t", escape_double = FALSE, trim_ws = TRUE)

x18S <- read_delim("18S_scores.txt",
                   "\t", escape_double = FALSE, trim_ws = TRUE, na = "NULL")
x28S <- read_delim("28S_scores.txt",
                   "\t", escape_double = FALSE, trim_ws = TRUE, na = "NULL")

#####Plotting#####

###Reads mapped

#absolute counts, by sample
ggplot(subset(reads_mapped_long,variable != "Percmapped")) +
  geom_bar(aes(Sample,value,fill=variable),stat="identity") +
  labs(y="Reads",fill=NULL,title="Reads Mapped per Sample") +
  scale_fill_viridis(discrete = TRUE, begin=0.1,end=0.6) +
  scale_y_continuous(labels=comma)
ggsave(filename="mapped_read_counts.png",device="png", dpi=100, width = 5.5 ,height = 4.5)

#proportion mapped, by sample
ggplot(subset(reads_mapped_long,variable != "Percmapped")) +
  geom_bar(aes(Sample,value,fill=variable),stat="identity",position="fill") +
  labs(y="Prop. Reads",fill=NULL,title="Reads Mapped to Transcriptome") +
  scale_fill_viridis(discrete = TRUE, begin=0.1,end=0.6)
ggsave(filename="mapped_read_props.png",device="png", dpi=100, width = 5.5 ,height = 4.5)

#absolute counts, by group
ggplot(subset(reads_mapped_long,variable != "Percmapped")) +
  geom_bar(aes(Group,value,fill=variable),stat="identity") +
  labs(y="Reads",fill=NULL,title="Reads Mapped per Sample") +
  scale_fill_viridis(discrete = TRUE, begin=0.1,end=0.6) +
  scale_y_continuous(labels=comma)
ggsave(filename="mapped_read_counts_grouped.png",device="png", dpi=100, width = 5.5 ,height = 4.5)

#proportion mapped, by group
ggplot(subset(reads_mapped_long,variable != "Percmapped")) +
  geom_bar(aes(Group,value,fill=variable),stat="identity",position="fill") +
  labs(y="Prop. Reads",fill=NULL,title="Reads Mapped to Transcriptome") +
  scale_fill_viridis(discrete = TRUE, begin=0.1,end=0.6)
ggsave(filename="mapped_read_props_grouped.png",device="png", dpi=100, width = 5.5 ,height = 4.5)

###Reads mapped to biotypes

#by sample
ggplot(biotypes) +
  geom_bar(aes(Biotype,ReadsMapped,fill=Sample),stat="identity",position="dodge") +
  coord_flip() +
  labs(title="Reads Mapping to RNA Biotypes",y="Reads Mapped") +
  scale_y_continuous(labels=comma)
ggsave(filename="mapped_reads_biotypes.png",device="png", dpi=100, width = 7, height = 6.5)

#by group
ggplot(biotypes) +
  geom_bar(aes(Biotype,ReadsMapped,fill=Group),stat="identity",position="dodge") +
  coord_flip() +
  labs(title="Reads Mapping to RNA Biotypes",y="Reads Mapped") +
  scale_y_continuous(labels=comma)
ggsave(filename="mapped_reads_biotypes_grouped.png",device="png", dpi=100, width = 7, height = 6.5)

###Correlations

#RPKM
ggplot(rpkmcorr) +
  geom_tile(aes(Sample1,Sample2,fill=r)) +
  geom_label(aes(Sample1,Sample2,label=round(r,2))) +
  labs(x=NULL,y=NULL,title="RPKM Correlation (Pearson's r)")
ggsave(filename="rpkm_corr_heat.png",device="png", dpi=100, width = 5.5, height = 4.5)

#RT stops
ggplot(rtcorr) +
  geom_tile(aes(Sample1,Sample2,fill=r)) + 
  geom_label(aes(Sample1,Sample2,label=round(r,2))) +
  labs(x=NULL,y=NULL,title="RT Stop Position Correlation (Pearson's r)")
ggsave(filename="rt_corr_heat.png",device="png", dpi=100, width = 5.5, height = 4.5)

###RT stop base composition

#By sample
ggplot(stopbases) +
  geom_bar(aes(Sample,Count,fill=Base),stat="identity",position="fill") +
  ggtitle("RT Stop Base Composition") +
  scale_fill_viridis(discrete = TRUE )
ggsave(filename="basecomp_bysample.png",device="png", dpi=100, width = 5.5, height = 4.5)

#By group
ggplot(stopbases) +
  stat_summary(aes(Base,Count,fill=Group),position="dodge",fun.y=sum,geom="bar") +
  ggtitle("RT Stop Base Composition") +
  scale_y_continuous(labels=comma) +
ggsave(filename="basecomp_bygroup.png",device="png", dpi=100, width = 5.5, height = 4.5)

###Enrichment scores by base

ggplot(subset(enrichbases_list)) +
  geom_histogram(aes(Score,fill=Base),position="fill",binwidth = 0.01)+
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("Base Count by Score")
ggsave(filename="enrichscores_bybase.png",device="png", dpi=100, width = 5.5, height = 4.5)

###Specific transcripts

#rRNA
ggplot(x18S) +
  geom_point(aes(Position,Score,color=Base)) +
  scale_color_viridis(discrete = TRUE) +
  ggtitle("18S rRNA")
ggsave(filename="enrich_18S.png",device="png", dpi=100, width = 15, height = 3)

ggplot(x28S) +
  geom_point(aes(Position,Score,color=Base)) +
  scale_color_viridis(discrete = TRUE) +
  ggtitle("28S rRNA")
ggsave(filename="enrich_28S.png",device="png", dpi=100, width = 15, height = 3)
