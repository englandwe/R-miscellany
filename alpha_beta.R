library(vegan)
library(ggplot2)
library(reshape2)
library(viridis)
source("/path/to/multiplot.R")
#web source for multiplot script: http://peterhaschke.com/Code/multiplot.R

###DATA###
#Data should be a matrix-style dataframe.  Samples as rows, OTUs as columns.  Names should be row/column labels.
#This works for a biom table converted from .biom to .txt using biom convert
indata <- t(read.csv("/path/to/OTUtable.txt",sep = "\t",header = TRUE,row.names = 1,skip = 1))

###METADATA###
#aka your mapping file
inmeta <- read.csv("/path/to/metadata.txt",sep = "\t",header = TRUE)

###RAREFACTION###
#Min num samples
raremax <- min(rowSums(indata))
#Generate rarefaction data.  Set step to something reasonable for your data.
rc_data <- rarecurve(indata, raremax, step = 1000, label=FALSE)
#Make a df for plotting with ggplot
#probably not the most R-esque way to do this since the word "apply" is not involved, but it works
rare_df=data.frame("Sample" = character(0), "Subsample" = integer(0), "TaxObs" = double(0), stringsAsFactors=FALSE)
for (i in seq_along(rc_data)) {
  for (j in seq_along(rc_data[[i]])) {
    rare_df[nrow(rare_df) + 1, ] <- c(row.names(indata)[[i]],attr(rc_data[[i]], "Subsample")[[j]],rc_data[[i]][[j]])
  }
}

###RAREFACTION - ADD METADATA###
rare_df_meta <- merge(rare_df,inmeta,by.x="Sample",by.y="X.SampleID",sort=FALSE)

###RAREFACTION PLOTS###
#Rarefaction curves
ggplot(rare_df_meta) +
  geom_line(aes(as.numeric(Subsample),as.numeric(TaxObs),group=Sample,color=D.T,linetype=Material),size=1) +
  scale_color_brewer(type="qualitative",palette = "Paired") +
  labs(x="Subsample Size",y="Contigs Observed") +
  theme_bw()

###ALPHA DIVERSITY###
alpha_list <- list()
alpha_list[[1]] <- diversity(indata, index="shannon")
alpha_list[[2]] <- diversity(indata, index="simpson")
alpha_list[[3]] <- fisher.alpha(indata ,MARGIN = 1)

#getting the data in a format ggplot likes
alpha_df <- as.data.frame(alpha_list)
alpha_df[,4] <- row.names(alpha_df)
colnames(alpha_df) <- c("Shannon","Simpson","Fisher","Sample")
alpha_df_long <- melt(alpha_df,id.vars = "Sample")

###ALPHA DIVERSITY - ADD METADATA###
alpha_df_long_meta <- merge(alpha_df_long,inmeta,by.x="Sample",by.y="X.SampleID",sort=FALSE)

###ALPHA DIVERSITY PLOTS###
ggplot(alpha_df_long_meta) +
  geom_point(aes(Sample,value,color=D.T,shape=Material),size=3) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_wrap(~variable,scales = "free") +
  scale_color_brewer(type="qualitative",palette = "Paired") +
  ggtitle("Alpha Diversity")

###BETA DIVERSITY###
#Calculate Bray-Curtis
#binary = FALSE gives you Bray-Curtis, TRUE would be Sorenson
beta_bc <- vegdist(indata, method="bray", binary=FALSE)
beta_bc_mds <- metaMDS(beta_bc, distance = "bray", k = 2, trymax = 20, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
beta_bc_mds_df <- as.data.frame(beta_bc_mds$points)
beta_bc_mds_stress <- beta_bc_mds$stress
beta_bc_mds_df[,3] <- row.names(beta_bc_mds_df)
colnames(beta_bc_mds_df)[[3]] <- "Sample"

###BETA DIVERSITY - ADD METADATA###
beta_bc_mds_df_meta <- merge(beta_bc_mds_df,inmeta,by.x="Sample",by.y="X.SampleID",sort=FALSE)

###BETA DIVERSITY - PLOTTING###
ggplot(beta_bc_mds_df_meta) +
  geom_point(aes(MDS1,MDS2,color=D.T,shape=Material),size=3) +
  theme_bw() + 
  annotate("text",label=paste("stress=",round(beta_bc_mds_stress,digits = 2),sep=""),x=Inf,y=Inf,hjust=1.1,vjust=1.5) +
  scale_color_brewer(type="qualitative",palette = "Paired") +
  labs(x="NMDS1",y="NMDS2")
