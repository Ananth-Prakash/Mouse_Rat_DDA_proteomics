library(ggplot2)
library(dplyr)
library(stats)
library(tidyr)
library(viridis)
library(grid)
library(RColorBrewer)
library(heatmap3)
library(stringr)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
library(sva)

BiocManager::install("limma")
library(limma)


# Input: ppb (FOT) normalised intensity values from MaxQuant,
#        Normal tissue samples only
#        Rename coloumns to add PXDID and Sample number
# 

setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/')

PXDID <- "PXD016958"

tmp  <- read.table(file = paste(PXDID, "/proteinGroups_ppb_final-tissue_names_", PXDID, ".txt", sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")


#Some of the gene entries (ex. IGHA2 has two Ensembl gene ids ENSG00000211890 & ENSG00000276173
#                          ex. IGHV2-70 has two Ensembl gene ids ENSG00000274576 & ENSG00000282453)
# because of this there are duplicate entries of such genes. These are aggregated by taking the median of them

colnames(tmp)[2] <- "Gene.Symbol"
dataset <- aggregate(tmp[,-c(1)], list(tmp$Gene.Symbol), median, na.rm =TRUE)
dataset <- subset(dataset, select = -c(Gene.Symbol))
colnames(dataset)[1] <- "Gene" 
dataset[dataset==0] <- NA

colnames(dataset)

colnames(dataset) <- gsub(".*_", "", colnames(dataset), perl=TRUE)
gene_names <- dataset[,c(1)] 


convert_to_scaled_data <- function(data, bins) {
  apply(data, 1, function(sample) {
    ceiling((sample / max(sample, na.rm = T)) * bins)
  })
}

convert_to_ranked_data <- function(data) {
  apply(data, 2, function(sample_data){
    unique_vals = unique(sample_data)
    unique_vals = unique_vals[order(unique_vals)]
    sapply(sample_data, function(s){which(s == unique_vals)}) - 1
  })
}

Binning <- function(data_input, gene_names) {
  # Make a data frame of sample + protein/peptide abundances. 
  data = data.frame(data_input)
  
  # Rank data.
  data_ranked = as.data.frame(t(convert_to_ranked_data(data)))
  colnames(data_ranked) <- gene_names
  
  # Different size of bins to try.
  binning_size_approaches = list(
    bins_2 = 1,
    bins_5 = 4,
    bins_10 = 9,
    bins_100 = 99,
    bins_1000 = 999,
    bins_10000 = 9999
  )
  
  # Transform data using different bin sizes.
  results = lapply(binning_size_approaches, function(binning_size_approach) {
    # 'Bin' data by scaling the ranked data.
    binned_data = as.data.frame(t(convert_to_scaled_data(data_ranked, binning_size_approach)))
  })
}

# If a dataset has different tissues or tissue regions, do not bin all together,
# these tissues/regions must be separated and binned separately.

colnames(dataset)

#binning_input <- dataset[,-c(1)]
binning_input <- dataset[,c(31:33), drop=FALSE]

#backup before changing NAs to 0
binning_input_with_NA <- binning_input

# non-detected genes (NA or 0 values) are assigned 0 values here to help in binning.
binning_input[is.na(binning_input)] <- 0


Dataset_all_bins <- Binning(binning_input, gene_names)


# Consider only data grouped into 5 bins
Get_Bin_5 <- function(binned_data_list){
  
  Bins_5 <- as.data.frame(t(binned_data_list$bins_5))
  # rename bins 0 to 4 --> 1 to 5
  # lowest intensities are in bin1 and highest in bin5
  Bins_5[Bins_5 == 4] <- 5
  Bins_5[Bins_5 == 3] <- 4
  Bins_5[Bins_5 == 2] <- 3
  Bins_5[Bins_5 == 1] <- 2
  Bins_5[Bins_5 == 0] <- 1
  
  Bins_5 <- tibble::rownames_to_column(Bins_5, "Gene")
}

Dataset_5_bins <- Get_Bin_5(Dataset_all_bins)


# Now reassign NA values back to those genes which were undetected or 0 in ALL samples 
# NOTE: Only those genes are reassigned back to NA that have missing values in
# all samples (ie., if at least one sample had a value then it is left as 0)
# this particular reassigning is done after binning (later in the code).


Dataset_5_bins[which(apply(binning_input_with_NA, 1, function(x) all(is.na(x)) )),-c(1)] <- NA

write.table(Dataset_5_bins, file = paste("Binned_expression_", PXDID ,"_Third-segment-of-proximal-tubule_Kidney.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

# Count the number of occurrences of each bin
Bins_count <- function(dataset){
  Bin_1 <- apply(dataset, 1, function(x) length(which(x == 1)) )
  Bin_2 <- apply(dataset, 1, function(x) length(which(x == 2)) )
  Bin_3 <- apply(dataset, 1, function(x) length(which(x == 3)) )
  Bin_4 <- apply(dataset, 1, function(x) length(which(x == 4)) )
  Bin_5 <- apply(dataset, 1, function(x) length(which(x == 5)) )
  
  result <- as.data.frame(cbind(Bin_1, Bin_2, Bin_3, Bin_4, Bin_5))
}

Bincounts <- cbind( Gene=Dataset_5_bins[,c(1)], as.data.frame( Bins_count ( as.data.frame(Dataset_5_bins[,-c(1)])  ) ) )
write.table(Bincounts, file = paste("Bin_counts_", PXDID, "_Third-segment-of-proximal-tubule_Kidney", sep=""), sep = "\t", row.names = FALSE, quote = FALSE )

# After this when comparing different datasets of the same tissue type, merge them
# this way there will be same genes between datasets when comparing

##### Comparing datasets

dataset1 <- read.table("Binned_expression_PXD001839_LeftVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset2 <- read.table("Binned_expression_PXD003375_Spinalcord-CaudalSegment.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset3 <- read.table("Binned_expression_PXD003375_Spinalcord-RostralSegment.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset4 <- read.table("Binned_expression_PXD004364_Testis.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset5 <- read.table("Binned_expression_PXD006692_Lung.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset6 <- read.table("Binned_expression_PXD012677_Amygdala.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset7 <- read.table("Binned_expression_PXD013543_Heart-LeftVentricle.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset8 <- read.table("Binned_expression_PXD015928_Tendon.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset9 <- read.table("Binned_expression_PXD016793_Liver.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset10 <- read.table("Binned_expression_PXD016958_Connecting-tubule_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset11 <- read.table("Binned_expression_PXD016958_Cortical-collecting-duct_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset12 <- read.table("Binned_expression_PXD016958_Cortical-thick-ascending-limb_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset13 <- read.table("Binned_expression_PXD016958_Distal-convoluted-tubule_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset14 <- read.table("Binned_expression_PXD016958_First-segment-of-proximal-tubule_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset15 <- read.table("Binned_expression_PXD016958_Inner-medullary-collecting-duct_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset16 <- read.table("Binned_expression_PXD016958_Medullary-thick-ascending-limb_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset17 <- read.table("Binned_expression_PXD016958_Outer-medullary-collecting-duct_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset18 <- read.table("Binned_expression_PXD016958_Second-segment-of-proximal-tubule_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
dataset19 <- read.table("Binned_expression_PXD016958_Third-segment-of-proximal-tubule_Kidney.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# Keep the minimum number of genes that are common in all datsets
Merge_data <- function(dataset1, dataset2){
  merged <- merge(x=dataset1, y=dataset2,
                  by.x=c("Gene"), by.y=c("Gene"), all.x=TRUE, all.y=TRUE)
}

merged_data <- Merge_data(dataset1, dataset2)
merged_data <- Merge_data(merged_data, dataset3)
merged_data <- Merge_data(merged_data, dataset4)
merged_data <- Merge_data(merged_data, dataset5)
merged_data <- Merge_data(merged_data, dataset6)
merged_data <- Merge_data(merged_data, dataset7)
merged_data <- Merge_data(merged_data, dataset8)
merged_data <- Merge_data(merged_data, dataset9)
merged_data <- Merge_data(merged_data, dataset10)
merged_data <- Merge_data(merged_data, dataset11)
merged_data <- Merge_data(merged_data, dataset12)
merged_data <- Merge_data(merged_data, dataset13)
merged_data <- Merge_data(merged_data, dataset14)
merged_data <- Merge_data(merged_data, dataset15)
merged_data <- Merge_data(merged_data, dataset16)
merged_data <- Merge_data(merged_data, dataset17)
merged_data <- Merge_data(merged_data, dataset18)
merged_data <- Merge_data(merged_data, dataset19)

filtered_data <- merged_data

# count number of samples for each gene that have non NA values (i.e., detected)
non_missing_sample_count_percentage <- data.frame(non_missing_sample_count_percentage= apply(filtered_data[2:ncol(filtered_data)], 1, function(x) (length(which(!is.na(x))))/(ncol(filtered_data)-1)*100 ))

filtered_data <- cbind(filtered_data, non_missing_sample_count_percentage)

# Filter genes that are present in at least 33%, 50% or 75% of samples
filtered_data <- filtered_data[filtered_data$non_missing_sample_count_percentage >= 50,]

### Since only common genes are considered, this condition is not required
# The undetected genes are also put into bin1 
# (After discussing with Andy Jones [20/01/2021]: since undetected genes could be below detection threshold due to their low expression or abundance)
filtered_data[is.na(filtered_data)] <- 1

write.table(filtered_data, file = "Binned_filtered_data_non_missing_sample_count_percentage_greaterthanequalto_50.txt", sep = "\t", row.names = FALSE, quote = FALSE )

pca_input <- filtered_data[,-c(1, ncol(filtered_data))]


#### PCA plot
# 
pca_bins <- prcomp(t(pca_input), scale = FALSE)
pca_plot_data <- data.frame(pca_bins$x[,1:2]) # Take components 1 and 2
pca_plot_data <- tibble::rownames_to_column(pca_plot_data, "Samples")

pca_plot_data$Organs <- gsub(".*\\.", "", pca_plot_data$Samples, perl=TRUE)
pca_plot_data$Organs <- gsub("Amygdala","Brain", pca_plot_data$Organs, perl=TRUE)
pca_plot_data$Organs <- gsub("LeftVentricle","Heart", pca_plot_data$Organs, perl=TRUE)
pca_plot_data$Organs <- gsub("CaudalSegment1|CaudalSegment2|CaudalSegment3|RostralSegment1|RostralSegment2|RostralSegment3","SpinalCord", pca_plot_data$Organs, ignore.case=T, perl=TRUE)
pca_plot_data$Organs <- gsub("Connectingtubule|Corticalcollectingduct|Corticalthickascendinglimb|Distalconvolutedtubule|Firstsegmentofproximaltubule|Innermedullarycollectingduct|Medullarythickascendinglimb|Outermedullarycollectingduct|Secondsegmentofproximaltubule|Thirdsegmentofproximaltubule","Kidney", pca_plot_data$Organs, ignore.case=T, perl=TRUE)

pca_plot_data <- pca_plot_data %>% group_by(Organs) %>% mutate( OrganID = cur_group_id() )

pca_plot_data$Datasets <- gsub("\\..*", "", pca_plot_data$Samples, perl=TRUE)
pca_plot_data <- pca_plot_data %>% group_by(Datasets) %>% mutate( DatasetID = cur_group_id() )


#Adding number of datasets next to each tissue sample
pca_plot_data <- pca_plot_data %>% group_by(Organs) %>% mutate( Organs = paste(Organs, "(", length(unique(Datasets)), ")") )

# To add lables on the legend
# From here: https://stackoverflow.com/questions/49965758/change-geom-texts-default-a-legend-to-label-string-itself

oldK <- GeomText$draw_key
# define new key
GeomText$draw_key <- function (data, params, size, 
                               var=unique(pca_plot_data$DatasetID), 
                               cols=scales::hue_pal()(length(var))) {
  
  # sort as ggplot sorts these alphanumerically / or levels of factor
  txt <- if(is.factor(var)) levels(var) else sort(var)
  txt <- txt[match(data$colour, cols)]
  
  textGrob(txt, 0.5, 0.5,  
           just="center", 
           gp = gpar(col = alpha(data$colour, data$alpha), 
                     fontfamily = data$family, 
                     fontface = data$fontface, 
                     fontsize = data$size * .pt))
}

ggplot(pca_plot_data, aes(x=PC1, y = PC2, colour = Datasets))+
  geom_text(aes(label=DatasetID))+
  labs(x="PC1", y="PC2")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(0.3,"line"))+
  guides(col = guide_legend(ncol = 3))+
  ggtitle("Rat-Samples bin_values_batch-per-dataset-regions\n[filter: genes detected in at least 50% of samples]")

# IMPORTANT reset key
GeomText$draw_key <- oldK

