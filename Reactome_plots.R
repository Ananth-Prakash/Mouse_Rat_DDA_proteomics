library(plotly)
library(dplyr)
library(stringi)

setwd('/Users/shengbo_wang/Documents/Projects/DDA/Mouse/Bin_intensity_analysis/')
#Load table:
Reactome_analysis <- read.csv("/Users/shengbo_wang/Documents/Projects/DDA/Mouse/Bin_intensity_analysis/Thawfeek_reactome-pathway-rat.csv",
                                sep = ",",
                              header = TRUE,
                                check.names = FALSE)
Reactome_analysis<-Reactome_analysis %>%
  filter(!(Organ==""))

Reactome_analysis<-Reactome_analysis %>%
  filter(!(Reactome_analysis$`Entities pValue`>=0.05))

Reactome_analysis<-Reactome_analysis %>%
  filter(!(Reactome_analysis$`Pathway name`=="Disease"))


new_df <- data.frame("Pathway_identifier" = rep(NA, length(unique(Reactome_analysis$`Pathway identifier`))),
                     "Pathway_name"= rep(NA, length(unique(Reactome_analysis$`Pathway identifier`))))

for (i in unique(Reactome_analysis$Organ)){
  new_df[[i]]<-NA
}

new_df$Pathway_identifier <- unique(Reactome_analysis$`Pathway identifier`)

for (i in seq(new_df$Pathway_identifier)){
  new_df$Pathway_name[i] <- Reactome_analysis$`Pathway name`[which(new_df$Pathway_identifier[i]==Reactome_analysis$`Pathway identifier`)[1]]
  
}

for (j in seq(dim(Reactome_analysis)[1])){
  
  new_df[which(new_df$Pathway_identifier==Reactome_analysis$`Pathway identifier`[j]),Reactome_analysis$Organ[j]] <- Reactome_analysis$`Entities pValue`[j]
}


#Remove Disease rows since they are not very informative


new_df<-new_df %>%
  filter(!(Pathway_name=="Disease"))

new_df <- new_df[new_df[,1]!="",]

new_df$Pathway_identifier <- na.omit(as.numeric(unlist(strsplit(new_df$Pathway_identifier, "R-HSA-"))))
new_df<- new_df %>% arrange(Pathway_identifier) %>% arrange(Pathway_name)
new_df$Pathway_identifier <- stri_sort(as.character(seq(new_df$Pathway_identifier)))



library(ggplot2)
library(tidyr)
library(heatmaply)
mine.long <- pivot_longer(data = new_df,
                          col = -c(1:2),
                          names_to = "Organ",
                          values_to = "p_value")

mine.heatmap <- ggplot(data = mine.long,
                       mapping = aes(x = Pathway_identifier,
                                     y = Organ,
                                     fill = p_value)) +
                geom_tile() +
                xlab(label = "Pathway id")
mine.heatmap

mine.heatmap <- ggplot(data = mine.long,
                       mapping = aes(x = Pathway_identifier,
                                     y = Organ,
                                     fill = p_value)) +
  geom_tile(colour = "white") +
  xlab("Pathway identifier") +
  scale_y_discrete(limits=rev)+
  theme(axis.ticks.x = element_line(colour =  rainbow(36)[as.numeric(as.factor(new_df$Pathway_name))],
                                    size = 1.5,
                                    lineend =  "square"),
        axis.ticks.length.x.bottom = unit(8, "pt"),
        axis.text.x = element_blank(),
        legend.position = "right")
mine.heatmap


#Print a legend without plot
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =unique(new_df$Pathway_name), pch=16, pt.cex=1.3, cex=0.7, bty='n',
       col = rainbow(36))
mtext("Pathway name", at=0.1, cex=1)


mine.heatmap <- ggplot(data = mine.long,
                       mapping = aes(x = Pathway_identifier,
                                     y = Organ,
                                     fill = p_value)) +
  geom_tile(colour = "lightgrey") +
  scale_fill_continuous(low="darkred", high="thistle2", 
                        guide="colorbar",na.value="white")+
  xlab("Pathway identifier") +
  scale_y_discrete(limits=rev)+
  labs(fill="p value")+
  
  theme(axis.ticks.x = element_line(colour =  rainbow(24)[as.numeric(as.factor(new_df$Pathway_name))],
                                    size = 1.5,
                                    lineend =  "square"),
        axis.ticks.length.x.bottom = unit(8, "pt"),
        axis.text.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = "right")
mine.heatmap


#Print a legend without plot
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend = unique(new_df$Pathway_name), pch=36, pt.cex=1, cex=0.7, bty='n',
       col = rainbow(36))
mtext("Pathway name", at=0.1, cex=1)


# Dendrogram figure for manuscript
# https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html

library(ggdendro)
library(ggplot2)
library(reshape2)
library(grid)

dendro_df <- new_df
#dendro_df[, c(3:34)] <- scale(dendro_df[, 3:34])
dendro_df.matrix <- as.matrix(dendro_df[, -c(1,2)])
rownames(dendro_df.matrix) <- dendro_df$Pathway_identifier
dendro_df.matrix[is.na(dendro_df.matrix)] <- 0
dendro_df.matrix_transpose <- t(dendro_df.matrix)
df.dendro <- as.dendrogram(hclust(d = dist(x = dendro_df.matrix_transpose)))
dendroplot <- ggdendrogram(data = df.dendro, rotate = TRUE)

new_long <- mine.long

reordered <- order.dendrogram(df.dendro)


new_long$accession <- factor(x = new_long$Organ,
                             levels = new_long$Organ[reordered], 
                             ordered = TRUE)

heatmap_plot <- ggplot(data = new_long, 
                       aes(x = Pathway_identifier, 
                           y = accession,
                           fill = p_value)) +
  geom_tile(colour = "lightgrey") +
  labs(y="")+
  scale_fill_continuous(low="darkred", high="thistle2", 
                        guide="colorbar",na.value="white")+
  xlab("Pathway identifier") +
  labs(fill="p value")+
  theme(axis.ticks.x = element_line(colour =  rainbow(36)[as.numeric(as.factor(new_df$Pathway_name))],
                                    size = 1.5,
                                    lineend =  "square"),
        axis.ticks.length.x.bottom = unit(8, "pt"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.position = "right")


grid.newpage()
print(heatmap_plot, vp = viewport(x = 0.6, y = 0.5, width = 0.8, height = 0.9995))
print(dendroplot+scale_y_reverse(), vp = viewport(x = 0.1, y = 0.5, width = 0.2, height = 1.05))
