# To get the median bin values of organs as Summary after binning tissues separately
library(RColorBrewer)
library(heatmap3)
library(pheatmap)
library(Hmisc)
setwd('/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/')

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

write.table(merged_data, file = "Binned_intensities_all_samples-RAT.txt", sep = "\t", row.names = FALSE, quote = FALSE )


binned_samples <- merged_data

colnames(binned_samples) <- gsub("PXD013543.Sample10.LeftVentricle", "PXD013543.Sample10a.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD013543.Sample12.LeftVentricle", "PXD013543.Sample12a.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD013543.Sample13.LeftVentricle", "PXD013543.Sample13a.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD013543.Sample14.LeftVentricle", "PXD013543.Sample14a.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD013543.Sample16.LeftVentricle", "PXD013543.Sample16a.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD013543.Sample21.LeftVentricle", "PXD013543.Sample21a.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD013543.Sample3.LeftVentricle", "PXD013543.Sample3a.LeftVentricle", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD013543.Sample8.LeftVentricle", "PXD013543.Sample8a.LeftVentricle", colnames(binned_samples))

colnames(binned_samples) <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                                  "Brain", colnames(binned_samples))

colnames(binned_samples) <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                                 "Heart", colnames(binned_samples))

colnames(binned_samples) <- gsub("PXD003375.Control1.CaudalSegment2", "PXD003375.Control1a.CaudalSegment1", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control1.CaudalSegment3", "PXD003375.Control2a.CaudalSegment1", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control2.CaudalSegment1", "PXD003375.Control3a.CaudalSegment1", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control2.CaudalSegment2", "PXD003375.Control4a.CaudalSegment2", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control2.CaudalSegment3", "PXD003375.Control5a.CaudalSegment3", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control3.CaudalSegment1", "PXD003375.Control6a.CaudalSegment1", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control3.CaudalSegment2", "PXD003375.Control7a.CaudalSegment2", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control3.CaudalSegment3", "PXD003375.Control8a.CaudalSegment3", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control1.RostralSegment1", "PXD003375.Control9a.RostralSegment1", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control1.RostralSegment2", "PXD003375.Control10a.RostralSegment2", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control1.RostralSegment3", "PXD003375.Control11a.RostralSegment3", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control2.RostralSegment1", "PXD003375.Control12a.RostralSegment1", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control2.RostralSegment2", "PXD003375.Control13a.RostralSegment2", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control2.RostralSegment3", "PXD003375.Control14a.RostralSegment3", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control3.RostralSegment1", "PXD003375.Control15a.RostralSegment1", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control3.RostralSegment2", "PXD003375.Control16a.RostralSegment2", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD003375.Control3.RostralSegment3", "PXD003375.Control17a.RostralSegment3", colnames(binned_samples))

colnames(binned_samples) <- gsub("PXD016958.Sample1.CorticalCollectingDuct", "PXD016958.Sample13a.CorticalCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.CorticalCollectingDuct", "PXD016958.Sample14a.CorticalCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.CorticalCollectingDuct", "PXD016958.Sample15a.CorticalCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample4.CorticalCollectingDuct", "PXD016958.Sample16a.CorticalCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample1.CorticalThickAscendingLimb", "PXD016958.Sample17a.CorticalThickAscendingLimb", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.CorticalThickAscendingLimb", "PXD016958.Sample18a.CorticalThickAscendingLimb", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.CorticalThickAscendingLimb", "PXD016958.Sample19a.CorticalThickAscendingLimb", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample1.DistalConvolutedTubule", "PXD016958.Sample20a.DistalConvolutedTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.DistalConvolutedTubule", "PXD016958.Sample21a.DistalConvolutedTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.DistalConvolutedTubule", "PXD016958.Sample22a.DistalConvolutedTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample1.FirstSegmentOfProximalTubule", "PXD016958.Sample23a.FirstSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.FirstSegmentOfProximalTubule", "PXD016958.Sample24a.FirstSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.FirstSegmentOfProximalTubule", "PXD016958.Sample25a.FirstSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample1.SecondSegmentOfProximalTubule", "PXD016958.Sample26a.SecondSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.SecondSegmentOfProximalTubule", "PXD016958.Sample27a.SecondSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.SecondSegmentOfProximalTubule", "PXD016958.Sample28a.SecondSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample1.ThirdSegmentOfProximalTubule", "PXD016958.Sample29a.ThirdSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.ThirdSegmentOfProximalTubule", "PXD016958.Sample30a.ThirdSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.ThirdSegmentOfProximalTubule", "PXD016958.Sample31a.ThirdSegmentOfProximalTubule", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample1.InnerMedullaryCollectingDuct", "PXD016958.Sample32a.InnerMedullaryCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.InnerMedullaryCollectingDuct", "PXD016958.Sample33a.InnerMedullaryCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.InnerMedullaryCollectingDuct", "PXD016958.Sample34a.InnerMedullaryCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample1.OuterMedullaryCollectingDuct", "PXD016958.Sample35a.OuterMedullaryCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.OuterMedullaryCollectingDuct", "PXD016958.Sample36a.OuterMedullaryCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.OuterMedullaryCollectingDuct", "PXD016958.Sample37a.OuterMedullaryCollectingDuct", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample1.MedullaryThickAscendingLimb", "PXD016958.Sample38a.MedullaryThickAscendingLimb", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample2.MedullaryThickAscendingLimb", "PXD016958.Sample39a.MedullaryThickAscendingLimb", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample3.MedullaryThickAscendingLimb", "PXD016958.Sample40a.MedullaryThickAscendingLimb", colnames(binned_samples))
colnames(binned_samples) <- gsub("PXD016958.Sample4.MedullaryThickAscendingLimb", "PXD016958.Sample41a.MedullaryThickAscendingLimb", colnames(binned_samples))


colnames(binned_samples) <- gsub("CaudalSegment1|CaudalSegment2|CaudalSegment3|RostralSegment1|RostralSegment2|RostralSegment3","SpinalCord", colnames(binned_samples), ignore.case=T, perl=TRUE)
colnames(binned_samples) <- gsub("Connectingtubule|Corticalcollectingduct|Corticalthickascendinglimb|Distalconvolutedtubule|Firstsegmentofproximaltubule|Innermedullarycollectingduct|Medullarythickascendinglimb|Outermedullarycollectingduct|Secondsegmentofproximaltubule|Thirdsegmentofproximaltubule","Kidney", colnames(binned_samples), ignore.case=T, perl=TRUE)


cor_results1 <- rcorr(as.matrix(binned_samples[,-c(1)]), type = c("pearson"))
correl_matrix <- cor_results1$r
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correl_summary <- flattenCorrMatrix(cor_results1$r, cor_results1$P) 
  
Pairwise_only_Brain <- with(correl_summary, correl_summary[ grepl("Brain", row) & grepl("Brain", column),])
Pairwise_only_Heart <- with(correl_summary, correl_summary[ grepl("Heart", row) & grepl("Heart", column),])
Pairwise_only_Kidney <- with(correl_summary, correl_summary[ grepl("Kidney", row) & grepl("Kidney", column),])
Pairwise_only_Lung <- with(correl_summary, correl_summary[ grepl("Lung", row) & grepl("Lung", column),])
Pairwise_only_SpinalCord <- with(correl_summary, correl_summary[ grepl("SpinalCord", row) & grepl("SpinalCord", column),])
Pairwise_only_Liver <- with(correl_summary, correl_summary[ grepl("Liver", row) & grepl("Liver", column),])

Significant <- correl_summary[correl_summary$p < 0.05,]

Organs <- data.frame(ID= factor(gsub(".*\\.", "", colnames(binned_samples)[-1], perl=TRUE)))

annotation <- data.frame(Organs = gsub(".*\\.", "", colnames(correl_matrix), perl=TRUE))
rownames(annotation) <- colnames(correl_matrix)
annotation$Datasets <- gsub("\\..*", "", rownames(annotation), perl=TRUE)


pheatmap(correl_matrix, show_rownames = F, show_colnames = F,
         annotation_col = annotation,
         border_color = NA,
         fontsize = 6)


merged_data_long <- gather(merged_data, Datasets, Bins, colnames(merged_data[2]):colnames(merged_data[ncol(merged_data)]), factor_key=TRUE)
merged_data_long$Tissues <- gsub(".*\\.", "", merged_data_long$Datasets)

merged_data_long$Organs <- merged_data_long$Tissues
merged_data_long$Organs <- gsub("Brain|DorsoLateralPreFrontalCortex|Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|Precuneus|CaudateNucleus|CerebellarHemisphericCortex|Cerebellum|CorpusCallosum|EntorhinalCortex|FrontalGyrus|InferiorParietalLobule|TemporalCortex|MiddleTemporalLobe|FrontalPole|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|Thalamus|VisualCortex",
                                "Brain", merged_data_long$Organs)
merged_data_long$Organs <- gsub("Heart|Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|PulmonaryValve|PulmonaryVein|RightAtrium|RightVentricle|TricuspidValve|VentricularSeptum",
                                "Heart", merged_data_long$Organs)
merged_data_long$Organs <- gsub("CaudalSegmentOfSpinalCord|RostralSegmentOfSpinalCord",
                                "SpinalCord", merged_data_long$Organs)
merged_data_long$Organs <- gsub("CaudalSegment1|CaudalSegment2|CaudalSegment3|RostralSegment1|RostralSegment2|RostralSegment3","SpinalCord", merged_data_long$Organs, ignore.case=T, perl=TRUE)
merged_data_long$Organs <- gsub("Connectingtubule|Corticalcollectingduct|Corticalthickascendinglimb|Distalconvolutedtubule|Firstsegmentofproximaltubule|Innermedullarycollectingduct|Medullarythickascendinglimb|Outermedullarycollectingduct|Secondsegmentofproximaltubule|Thirdsegmentofproximaltubule","Kidney", merged_data_long$Organs, ignore.case=T, perl=TRUE)

merged_data_long_aggregate <- aggregate(merged_data_long[,-c(2,4)], list("Gene" = merged_data_long$Gene, "Organs" = merged_data_long$Organs), median, na.rm =TRUE)

merged_data_long_aggregate <- merged_data_long_aggregate[,-c(3,5)]
colnames(merged_data_long_aggregate) <- c("Gene", "Tissues", "Median_bins")

All_tissues_median_bins<- spread(merged_data_long_aggregate, Tissues, Median_bins)

colnames(All_tissues_median_bins[,-c(1)]) <- gsub("^(.*)", "Median_bins_", colnames(All_tissues_median_bins[,-c(1)]), perl=TRUE)

Gene_info <- read.table("Gene_distribution_in_organs-GeneNames-Median_intensities-RAT.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
Gene_info <- Gene_info[,c(1,2)]

All_tissues_median_bins <- merge(x=Gene_info, y=All_tissues_median_bins,
                                 by.x=c("GeneName"), by.y=c("Gene"))

write.table(All_tissues_median_bins, file = "Gene_distribution_in_organs-GeneNames-Median_bin_values-RAT.txt", sep = "\t", row.names = FALSE, quote = FALSE )
