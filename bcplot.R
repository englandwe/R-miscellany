library(readr)
library(ggplot2)
library(scales)

#barcodes.txt: output of splitFastq.pl run through bcCount.py

barcodes <- read_delim("barcodes.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(barcodes) <- c("Barcode","Count","BarcodeID")

bcplot <- ggplot(barcodes) +
  geom_bar(aes(1,Count,fill=BarcodeID),stat="identity") +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_y_continuous(labels=comma) +
  labs(x=NULL,y="Reads",fill="Barcode") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(filename="barcode_counts.png",plot=bcplot,device="png", dpi=100, width = 3.5 ,height = 5.5)

