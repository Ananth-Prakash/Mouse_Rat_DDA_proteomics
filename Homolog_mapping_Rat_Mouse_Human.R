# Heatmap showing homolog bin values in tissues across mouse human and rat

library(ggplot2)
library(reshape2)
library(viridis)
library(stringr)
library(dplyr)
library(plyr)
library(tidyr)
library(pagenum)
library(ggpubr)
library(gridExtra)

homolog_bins <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/homologs_mouse_rat_human-bins-one_to_one_mappings_only_common_organs.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# to check how many one-to-one mappings by counting number of times gene names have occurred in the dataset
# mapping count 1 in all three genes are one-to-one mappings and anything else are either one-to-many or many-to-many mappings
mappings_human <- ddply(homolog_bins,.(GeneName_human),nrow)
mappings_mouse <- ddply(homolog_bins,.(GeneName_mouse),nrow)
mappings_rat <- ddply(homolog_bins,.(GeneName_rat),nrow)

colnames(mappings_human) <- c("GeneName_human", "GeneName_counts_human")
colnames(mappings_mouse) <- c("GeneName_mouse", "GeneName_counts_mouse")
colnames(mappings_rat) <- c("GeneName_rat", "GeneName_counts_rat")

homolog_bins <- merge(homolog_bins, mappings_human,
                      by.x=c("GeneName_human"), by.y=c("GeneName_human"),
                      all.x=FALSE, all.y=FALSE)
homolog_bins <- merge(homolog_bins, mappings_mouse,
                      by.x=c("GeneName_mouse"), by.y=c("GeneName_mouse"),
                      all.x=FALSE, all.y=FALSE)
homolog_bins <- merge(homolog_bins, mappings_rat,
                      by.x=c("GeneName_rat"), by.y=c("GeneName_rat"),
                      all.x=FALSE, all.y=FALSE)

# Filter down to include only one-to-one mappings
homolog_bins <- homolog_bins[(homolog_bins$GeneName_counts_human ==1 & homolog_bins$GeneName_counts_mouse == 1 & homolog_bins$GeneName_counts_rat == 1),]

homolog_bins <- homolog_bins[order(homolog_bins$GeneName_human),]
homolog_bins$median <- apply(homolog_bins[,c(8:25)], 1, median, na.rm = T)

homolog_bins <- homolog_bins[homolog_bins$Present_in_samples != 0,]
data <- homolog_bins[,-c(1:7,26:30)]

data$Homologs <- factor(paste(homolog_bins$GeneName_human, homolog_bins$UniProt_ID_human, homolog_bins$mouse_ensg, homolog_bins$human_ensg, homolog_bins$rat_ensg, sep=" / "))


data_long <- gather(data, Organs, Bins, Brain_mouse:Testis_rat, factor_key=TRUE)
data_long$Species <- gsub(".*_","", data_long$Organs, perl=TRUE)
data_long$Organs <- gsub("_.*","", data_long$Organs, perl=TRUE)

data_long$Species <- gsub("mouse","Mouse",data_long$Species)
data_long$Species <- gsub("human","Human",data_long$Species)
data_long$Species <- gsub("rat","Rat",data_long$Species)


data_summary <- data
data_summary$Mouse_all_median <- apply(data_summary[,grepl("mouse", colnames(data_summary))],1,median, na.rm =TRUE)
data_summary$Human_all_median <- apply(data_summary[,grepl("human", colnames(data_summary))],1,median, na.rm =TRUE)
data_summary$Rat_all_median <- apply(data_summary[,grepl("rat", colnames(data_summary))],1,median, na.rm =TRUE)

data_summary$Mouse_count <- apply(data_summary[,grepl("_mouse", colnames(data_summary))],1,function(x) (length(which(!is.na(x)))))
data_summary$Human_count <- apply(data_summary[,grepl("_human", colnames(data_summary))],1,function(x) (length(which(!is.na(x)))))
data_summary$Rat_count <- apply(data_summary[,grepl("_rat", colnames(data_summary))],1,function(x) (length(which(!is.na(x)))))

# A subdata where orthologs are expressed is expressed in human, rat and mouse 
subdata <- data_summary[(data_summary$Human_count >= 3 & data_summary$Mouse_count >= 3  & data_summary$Rat_count >= 3 ),]

#take random 10 entries
subdata <- subdata[sample(nrow(subdata), 10), ]
subdata_long <- gather(subdata, Organs, Bins, Brain_human:Testis_rat, factor_key=TRUE)
subdata_long$Species <- gsub(".*_","", subdata_long$Organs, perl=TRUE)
subdata_long$Organs <- gsub("_.*","", subdata_long$Organs, perl=TRUE)

subdata_long$Species <- gsub("mouse","Mouse",subdata_long$Species)
subdata_long$Species <- gsub("human","Human",subdata_long$Species)
subdata_long$Species <- gsub("rat","Rat",subdata_long$Species)

subdata_long$Homologs <- gsub("\\/.*","", subdata_long$Homologs, perl=TRUE)
ggplot(subdata_long, aes(x=Organs,y=Homologs))+
  geom_tile(alpha=0, color = "gray",
            lwd = 0.3,
            linetype = 1)+
  geom_point(aes(colour = Species, 
                 size =Bins), position=position_dodge(width=0.5))+
  scale_y_discrete(limits=rev)+
  theme_bw()

x=0
#setPagenum(1)
input <- data.frame()
plot_list = list()

no_of_rows_per_plotpage <- 50
no_of_plot_pages <- ceiling(nrow(homolog_bins)/no_of_rows_per_plotpage)
  
for (i in 1:no_of_plot_pages) {
  
  subdata<-tail(head(data, n=50+x), n=50)
  
  input <- rbind(input,subdata)
  x=nrow(input)
  
  input_long <- gather(subdata, Organs, Bins, Brain_human:Testis_rat, factor_key=TRUE)
  input_long$Species <- gsub(".*_","", input_long$Organs, perl=TRUE)
  input_long$Organs <- gsub("_.*","", input_long$Organs, perl=TRUE)
  
  input_long$Species <- gsub("mouse","Mouse",input_long$Species)
  input_long$Species <- gsub("human","Human",input_long$Species)
  input_long$Species <- gsub("rat","Rat",input_long$Species)
  
  p = ggplot(input_long,aes(x=Organs,y=Homologs))+
    geom_tile(alpha=0, color = "gray",
              lwd = 0.3,
              linetype = 1)+
    geom_point(aes(colour = Species, 
                   size =Bins), position=position_dodge(width=0.5))+
    scale_y_discrete(limits=rev)+
    theme_bw()
    #pagenum(text="ABC Corp - ",x=.95, y=.2, just=c('right','bottom'))
  plot_list[[i]] = p
}

# Another option: create pdf where each page is a separate plot.
pdf(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/Plots/Homolog_binned_expression_comparison.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 15) # The height of the plot in inches
for (i in 1:no_of_plot_pages) {
  print(plot_list[[i]])
}

# Step 3: Run dev.off() to create the file!
dev.off()


#### To plot pairwise scatterplot
homolog_ppb <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/homologs_mouse_rat_human-ppb-one_to_one_mappings_only_common_organs.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


data_ppb <- homolog_ppb[,-c(1:6)]

#data$mouse_ensg <- factor(data$mouse_ensg)
data_ppb$Homologs <- factor(paste(homolog_ppb$GeneName_human, homolog_ppb$mouse_ensg, homolog_ppb$human_ensg, homolog_ppb$rat_ensg, sep=" / "))

mouse_long <- gather(data_ppb, Organs, Mouse, Brain_median_mouse:Testis_median_mouse, factor_key=TRUE)
mouse_long$Organs <- gsub("_.*","", mouse_long$Organs, perl=TRUE)
mouse_long <- mouse_long[,c("Homologs","Organs","Mouse")]


human_long <- gather(data_ppb, Organs, Human, Brain_median_human:Testis_median_human, factor_key=TRUE)
human_long$Organs <- gsub("_.*","", human_long$Organs, perl=TRUE)
human_long <- human_long[,c("Homologs","Organs","Human")]

rat_long <- gather(data_ppb, Organs, Rat, Brain_median_rat:Testis_median_rat, factor_key=TRUE)
rat_long$Organs <- gsub("_.*","", rat_long$Organs, perl=TRUE)
rat_long <- rat_long[,c("Homologs","Organs","Rat")]

#1
Scatter_mouse_human <- merge(x=mouse_long, y=human_long,
                             by.x=c("Homologs","Organs"), by.y=c("Homologs","Organs"))

ggplot(Scatter_mouse_human , aes(x=Mouse, y=Human)) + 
  geom_bin2d(bins = 60) +
  scale_fill_continuous(type = "viridis")  + 
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.6,
           label.y.npc = 0.2)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Homologs_ppb_expression_Mouse_Human")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  facet_wrap(~Organs)

#2
Scatter_rat_human <- merge(x=rat_long, y=human_long,
                             by.x=c("Homologs","Organs"), by.y=c("Homologs","Organs"))

ggplot(Scatter_rat_human , aes(x=Rat, y=Human)) + 
  geom_bin2d(bins = 60) +
  scale_fill_continuous(type = "viridis")  + 
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.6,
           label.y.npc = 0.2)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Homologs_ppb_expression_Rat_Human")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  facet_wrap(~Organs)

#2
Scatter_rat_mouse <- merge(x=rat_long, y=mouse_long,
                           by.x=c("Homologs","Organs"), by.y=c("Homologs","Organs"))

ggplot(Scatter_rat_mouse , aes(x=Rat, y=Mouse)) + 
  geom_bin2d(bins = 60) +
  scale_fill_continuous(type = "viridis")  + 
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.6,
           label.y.npc = 0.2)+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Homologs_ppb_expression_Rat_Mouse")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  facet_wrap(~Organs)


#### To compare correlation within species across organs

human_ppb <- data_ppb[,c(1:6)]
mouse_ppb <- data_ppb[,c(7:12)]
rat_ppb <- data_ppb[,c(13:18)]

colnames(human_ppb) <- gsub("_median*","", colnames(human_ppb), perl=TRUE)
colnames(mouse_ppb) <- gsub("_median*","", colnames(mouse_ppb), perl=TRUE)
colnames(rat_ppb) <- gsub("_median*","", colnames(rat_ppb), perl=TRUE)

human_var_list = combn(names(human_ppb)[1:6], 2, simplify=FALSE)
mouse_var_list = combn(names(mouse_ppb)[1:6], 2, simplify=FALSE)
rat_var_list = combn(names(rat_ppb)[1:6], 2, simplify=FALSE)

# Make plots.
plot_list_human = list()
for (i in 1:15) {
  p = ggplot(human_ppb, aes_string(x=human_var_list[[i]][1], y=human_var_list[[i]][2])) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis")  + 
    theme_bw()+
    scale_x_log10()+
    scale_y_log10()+
    stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.6,
             label.y.npc = 0.2)+
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.position = "none")+
    # ggtitle("Human")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  plot_list_human[[i]] = p
}
do.call(grid.arrange,plot_list_human)

#2
plot_list_mouse = list()
for (i in 1:15) {
  p = ggplot(mouse_ppb, aes_string(x=mouse_var_list[[i]][1], y=mouse_var_list[[i]][2])) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis")  + 
    theme_bw()+
    scale_x_log10()+
    scale_y_log10()+
    stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.6,
             label.y.npc = 0.2)+
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.position = "none")+
    # ggtitle("Mouse")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  plot_list_mouse[[i]] = p
}
do.call(grid.arrange,plot_list_mouse)

#3
plot_list_rat = list()
for (i in 1:15) {
  p = ggplot(rat_ppb, aes_string(x=rat_var_list[[i]][1], y=rat_var_list[[i]][2])) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis")  + 
    theme_bw()+
    scale_x_log10()+
    scale_y_log10()+
    stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.6,
             label.y.npc = 0.2)+
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.position = "none")+
    # ggtitle("Rat")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  plot_list_rat[[i]] = p
}
do.call(grid.arrange,plot_list_rat)
