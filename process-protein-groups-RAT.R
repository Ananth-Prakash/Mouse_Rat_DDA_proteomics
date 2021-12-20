library(matrixStats)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(mygene)

setwd("/Users/ananth/Documents/MaxQuant_Bechmarking/Rat/PXD016958/")


if ( !exists("EXP_TYPE")) warning("Please specify experiment type variable: EXP_TYE")
## Specify experiment type, i.e. what type of quantification is used, MS1-based or MS2-based. As a rule of thumb, label free and SILAC is MS1-based, and iTRAQ and TMT is MS2-based.
# EXP_TYPE <- "MS2-quant"
EXP_TYPE <- "MS1-quant"


##############
# define FOT normalisation function
# FOT stands for Fraction Of Total. 
# In this normalisation method each protein iBAQ intensity value is scaled to the total amount of signal in a given MS run (column) and transformed to parts per billion (ppb)
fot.normalise <- function(x){
  data.sum <-   apply(x, 2, function(y){sum(y, na.rm=TRUE)})
  # barplot((data.sum), log = "y")
  #
  ##### do ppm normalisation
  x.mat <- as.matrix(x)
  x.mat.ppb <- apply(x.mat, 2, function(i) i/sum(i, na.rm = T) * 1000000000 )
  x.mat.ppb <- as.data.frame(x.mat.ppb)
  colnames(x.mat.ppb) <- paste("ppb.", colnames(x.mat.ppb), sep = "")
  return(x.mat.ppb)
}
##############


##### read  protein expression matrix produced by MaxQuant
tmp  <- read.table( "proteinGroups.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

## clean up 
tmp <- tmp[ tmp[, "Reverse"] != "+", ]
tmp <- tmp[ tmp[, "Potential.contaminant"] != "+", ]
## consider those protein groups with more than 1 peptides matched
tmp <- tmp[ tmp[, "Peptides"] > 1, ]

## 
# if experiment is label free use iBAQ quant:
if(EXP_TYPE == "MS1-quant"){
  message("Collecting iBAQ quantification")
  tmp <- tmp[ ,c(2, grep("iBAQ.", colnames(tmp))) ]
}


# if experiment is TMT, use reporter intensities, if experiment iTRAQ use Intensities. Note, MaxQuant might change how it reports iTRAQ and TMT intensities.
if(EXP_TYPE == "MS2-quant"){
  message("Collecting MS2 intensities")
  if( any(grepl("Reporter.intensity.corrected", colnames(tmp))) ){
    
    # Filter to consider only protein groups for which the at least 50 % of the sample replicates 
    # have non-zero intensity values
    
    ##### Read sample replicate information
    ##### Note: User defined annotation file
    sample_replicates <- read.table("sample_replicates.txt" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
    
    for(i in 1:nrow(sample_replicates))
    {
      sample_replicates[i,] <- gsub("^","Reporter.intensity.corrected.", sample_replicates[i,], perl=TRUE)
      sample_replicates[i,] <- gsub(" ",".", sample_replicates[i,], perl=TRUE)
      sample_replicates[i,] <- gsub("-",".", sample_replicates[i,], perl=TRUE)
      
      subdata <- tmp[, c(unname(unlist(sample_replicates[i,-1]))), drop=FALSE]
      tmp$replicate_nonzero_instensity_count <- apply(subdata[1:ncol(subdata)], 1, function(x) length(which(x != 0)) )
      
      tmp <- tmp[tmp$replicate_nonzero_instensity_count >= ncol(subdata)/2,]
      
    }

    tmp <- tmp[ ,c(2, grep("Reporter.intensity.corrected.[0-9].{1,}.", colnames(tmp))), drop=FALSE ]
  } else {
    tmp <- tmp[ ,c(2, grep("Intensity.", colnames(tmp))), drop=FALSE ]
  }
}


#####
Majority.protein.IDs <- tmp$Majority.protein.IDs
tmp <- tmp[ , -1]

if(EXP_TYPE == "MS1-quant"){
  tmp <- fot.normalise(tmp)  
}
#
tmp <- data.frame( cbind(Majority.protein.IDs, tmp, stringsAsFactors = FALSE) )
##############
tmp[tmp == 0] <- NA
tmp[ tmp == "NaN"] <- NA


data.to.map <- tmp
data.to.map$"ENSG" <- "NA"
data.to.map$"Gene.Symbol" <- "NA"
data.to.map$"unique.gene.count" <- "NA"

for(i in 1:nrow(data.to.map)){
  
  x <- data.frame(strsplit(data.to.map[ i, "Majority.protein.IDs"], split = ";"), stringsAsFactors = FALSE)
  x_temp <- regmatches(x[,1],regexpr("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", x[,1]))

  f = file()
  sink(file=f)
  res <- tryCatch(queryMany(x_temp, scopes="uniprot", fields=c("ensembl.gene", "symbol"), species=c("rat")),
                  error = function(e) {print(0)})
  sink()
  close(f)
  
  if (class(res)=="DFrame" | class(res) == "DataFrame"){
    data.to.map[ i, "ENSG"] <- paste( unique(unlist(res$ensembl.gene[!is.na(res$ensembl.gene)])), collapse = ";")
    data.to.map[ i, "Gene.Symbol"] <- paste( unique(unlist(res$symbol[!is.na(res$symbol)])), collapse = ";")}
  temp_symb <- data.to.map[i,"Gene.Symbol"]
  data.to.map[ i , "unique.gene.count"] <- str_count(unique(temp_symb), ";")+1
  
  print(paste0("Processing protein groups... ", as.character(round(i*100/nrow(data.to.map),1)),"%"))
}

# data backup before filtering
data.to.map_before_filtering <- data.to.map

data.to.map <- data.to.map[ data.to.map$ENSG != "" , ]
data.to.map <- data.to.map[ data.to.map$ENSG != "NA" , ]
data.to.map <- data.to.map[ data.to.map[, "unique.gene.count"] == 1, ]

# remove all protein groups that map to multiple ENSG gene IDs (this removes a lot of proteins) - the reasoning to remove these cases is that we cannot establish for sure which gene is contributing the signal to the protein abundance; all genes contribute equally or one gene is a majority? 
data.to.map <- data.to.map[ grep(";", data.to.map$ENSG, invert = TRUE) , ]

colnames(data.to.map)
##aggregating columns
data.to.map <- aggregate(data.to.map[ , 2:(ncol(data.to.map)-3)], list("Gene ID" = data.to.map$ENSG, "Gene.Symbol" = data.to.map$Gene.Symbol), median, na.rm =TRUE)

colnames(data.to.map) <- gsub("Reporter.intensity.corrected.[0-9].", "", colnames(data.to.map), perl=TRUE)

if(EXP_TYPE == "MS2-quant"){
  write.table(data.to.map, "proteinGroups_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )} else {
    write.table(data.to.map, "proteinGroups_ppb_final.txt", sep = "\t", row.names = FALSE, quote = FALSE )   
  }

#Write to file the list of protein groups which are mapped to more than 1 gene ID
ambiguous_gene_mapped_protein_groups <- data.to.map_before_filtering[ data.to.map_before_filtering[, "unique.gene.count"] > 1 | data.to.map_before_filtering[, "ENSG"] == "" | is.na(data.to.map_before_filtering$ENSG), ]
ambiguous_gene_mapped_protein_groups <- data.to.map_before_filtering[ data.to.map_before_filtering[, "unique.gene.count"] > 1 ,]
ambiguous_gene_mapped_protein_groups$ENSG <- unlist(lapply(ambiguous_gene_mapped_protein_groups$ENSG,function(x) paste(sort(strsplit(x, split=";")[[1]]), collapse=";")))
ambiguous_gene_mapped_protein_groups$Gene.Symbol <- unlist(lapply(ambiguous_gene_mapped_protein_groups$Gene.Symbol,function(x) paste(sort(strsplit(x, split=";")[[1]]), collapse=";")))

ambiguous_gene_mapped_protein_groups <- ambiguous_gene_mapped_protein_groups[,c(1, ncol(ambiguous_gene_mapped_protein_groups)-2, ncol(ambiguous_gene_mapped_protein_groups)-1, ncol(ambiguous_gene_mapped_protein_groups), 2:(ncol(ambiguous_gene_mapped_protein_groups)-3) )]

write.table(ambiguous_gene_mapped_protein_groups, "ambiguous_gene_mapped_protein_groups.txt", sep = "\t", row.names = FALSE, quote = FALSE )
