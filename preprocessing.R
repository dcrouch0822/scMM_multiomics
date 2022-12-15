###########################################################
#   Script to run scMultiOMIC-seq preprocessing workflow  #
#                       D.Croucher                        #
#                        Feb 2022                         #
###########################################################

###########################################################


### ------------ GENERAL OVERVIEW OF THIS SCRIPT
### Load packages 
### Parse Options
### Pipeline Variables
### Collect sequencing metrics outputted by cellranger multi -- omit for now
### Read count data and merge into single seurat object 
### Add metadata
### Add cell-level quality control metrics 
### Add QC flag
### Read VDJ data and add to seurat object


### ------------ THIS PIPELINE OUTPUTS THE FOLLOWING:
### Summary of sequencing metrics (csv) -- omit for now
### Sample-specific, unfiltered seurat object with QC flags for nGenes and percentMito. 
### VDJ info for BCR and TCR data has been added to meta.data slot
### Sample-specific VDJ metadata (csv) 


### ------------ PENDING
### Sequencing metrics from cellranger
### Doublet calling
### Background removal
### Plots


###########################################################

###########################################################
### EXAMPLE EXECUTION ON H4H
##
## #!/bin/bash
## #SBATCH -t 72:00:00
## #SBATCH --mem=60G
## #SBATCH -p himem
## #SBATCH -c 15
## #SBATCH -N 1
##
## module load R/4.1.0
##
## Rscript /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/scripts/preprocessing.R 
##    --cellrangerOutputDir /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/cr_outs/ \
##    --rawDGEinputDir /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/cr_outs/ \
##    --outputDir /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/analysis/preprocess/ \
##    --samples /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/analysis/ALQ_MCRN007_BCMA_sampleNames.txt \
##    --sampleMetaData /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/analysis/ALQ_MCRN007_BCMA_sampleMetaData.csv \
##    --minGenesPerCell 200 \
##    --maxMitoPerCell 20 \
##    --runDoubletFinder TRUE \
##    --outputPlots TRUE 
##
###########################################################

### start clocks
StopWatchStart <- list()
StopWatchEnd   <- list()
StopWatchStart$Overall  <- Sys.time()
StopWatchEnd$Overall  <- Sys.time()

#########################################
# Load packages
#########################################

print("")
print("********************")
print("Load packages")
print(Sys.time())
print("********************")
print("")
StopWatchStart$LoadPackages  <- Sys.time()


suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
#suppressMessages(library(doubletFinder))
suppressMessages(library(Matrix))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(reshape))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gplots))
suppressMessages(library(readxl))
suppressMessages(library(stringr))
suppressMessages(library(Hmisc))
suppressMessages(library(devtools))
suppressMessages(library(plyr))
suppressMessages(library(optparse))


StopWatchEnd$LoadPackages  <- Sys.time()

#########################################
#  Parse Options
#########################################

option_list <- list(make_option("--cellrangerOutputDir",
                                type = "character",
                                help = "path to dir with cr count outs to read in",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--rawDGEinputDir",
                                type = "character",
                                help = "same as cellrangerOutputDir",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--outputDir",
                                type = "character",
                                help = "path to dir where all pipeline outputs will go",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--samples",
                                type = "character",
                                default = NULL,
                                help = "path to txt file with column SAMPLEID; each sample on own line",
                                metavar= "character"
                               ),
                    make_option("--sampleMetaData",
                                type = "character",
                                help = "path to dir where sample level metadata file is located",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--minGenesPerCell",
                                type = "integer",
                                help = "minimum number of genes per cell - will set whether QC flag is added but won't filter out cells",
                                default = NULL,
                                metavar= "integer"
                               ),
                    make_option("--maxMitoPerCell",
                                type = "integer",
                                help = "maximum percentage of mitochondrial UMIs per cell - will set whether QC flag is added but won't filter out cells",
                                default = NULL,
                                metavar= "integer"
                               ),     
                    make_option("--runDoubletFinder",
                                type = "logical",
                                help = "TRUE/FALSE to run DoubletFinder - will add doublet status in meta.data but won't filter out cells",
                                default = NULL,
                                metavar= "logical"
                               ),
                    make_option("--outputPlots",
                                type = "logical",
                                help = "TRUE/FALSE to output plots",
                                default = NULL,
                                metavar= "logical"
                               )                    
                   )


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cellrangerOutputDir <- opt$cellrangerOutputDir
rawDGEinputDir <- opt$rawDGEinputDir
outputDir <- opt$outputDir
samples <- opt$samples
sampleMetaData <- opt$sampleMetaData
minGenesPerCell <- opt$minGenesPerCell
maxMitoPerCell <- opt$maxMitoPerCell
runDoubletFinder <- opt$runDoubletFinder
outputPlots <- opt$outputPlots

#########################################
# Pipeline Variables
#########################################

print("")
print("********************")
print("Pipeline Variables")
print(Sys.time())
print("********************")
print("")
StopWatchStart$SetUpVaribles  <- Sys.time()

sampleNames <- c(as.character(read.table(samples, header = T)$SAMPLEID))


if(length(sampleNames) == 1){
    print(paste0("Sample IDs: ", sampleNames))
} else if(length(sampleNames > 1)){
    print(paste0("Sample IDs: ", paste(sampleNames, collapse = ", ")))
}

  

dataDirs <- c(paste(cellrangerOutputDir,
                    sampleNames, 
                    sep = ""))

dataDirCounts <- c(paste(dataDirs,
                         "/count/sample_feature_bc_matrix/", 
                         sep = ""))

dataDirTCR <- c(paste(dataDirs, 
                      "/vdj_t/", 
                      sep = ""))

dataDirBCR <- c(paste(dataDirs, 
                      "/vdj_b/", 
                      sep = ""))

outputPreprocessDirs <- c(paste(sampleNames, 
                      "/vdj_b/", 
                      sep = ""))

StopWatchEnd$SetUpVaribles  <- Sys.time()


for (name in sampleNames) {

 setwd(outputDir)
 dir.create(name)
 setwd(name)
 dir.create("data")
 if (outputPlots == TRUE){dir.create("figures")}
  
}



#########################################
# Collect sequencing metrics outputted by cellranger multi -- omit for now
#########################################

          
                       
#########################################
# Read count data and merge into single seurat object 
#########################################

print("")
print("*********************")
print("Read data and merge")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ReadDataMerge <- Sys.time()

#Load and merge counts into a single seurat object

for (i in seq_along(sampleNames)) {
  
  data <- Read10X(dataDirCounts[i], strip.suffix = TRUE) #returns the a list with counts matrices for GEX and ADT
  colnames(data$`Gene Expression`) <- paste(sampleNames[i], colnames(data$`Gene Expression`), sep="_") 
  colnames(data$`Antibody Capture`) <- paste(sampleNames[i], colnames(data$`Antibody Capture`), sep="_") 
  
  seuratTmp <- CreateSeuratObject(counts = data$`Gene Expression`)
  seuratTmp[['ADT']] <- CreateAssayObject(counts = data$`Antibody Capture`)
  
    if (i == as.integer(1)){

      seurat <- seuratTmp
    
    }
  
    if (i > as.integer(1)){
      
      seurat <- merge(x=seurat, y=seuratTmp)
    }

}


adt.proteins <- rownames(data$`Antibody Capture`)
length(adt.proteins)

# Validate that the object now contains multiple assays
Assays(seurat)

# List the current default assay
DefaultAssay(seurat)

# Clean up 
rm(seuratTmp)


StopWatchEnd$ReadDataMerge <- Sys.time()                     
                       
                       
#########################################
# Add metadata
#########################################

print("")
print("*********************")
print("Add metadata")
print(Sys.time())
print("********************")
print("")
StopWatchStart$AddMetaData <- Sys.time()


seurat@meta.data$Subject_ID <- gsub("-T.*", "", seurat@meta.data$orig.ident)
seurat@meta.data$Timepoint <- paste0("T", gsub(".*-T", "", seurat@meta.data$orig.ident)) #consider gsubing the -O part of (says T0-O) 
seurat@meta.data$Sample_ID <- seurat@meta.data$orig.ident


# Add clinical meta data 

sample.meta.data <- read.csv(sampleMetaData)

samples_for_meta <- unique(seurat@meta.data$Sample_ID) #vector of sample names which is added to "orig.ident" in previous steps
meta.name <- colnames(sample.meta.data) #name of sample-level metadata (columns)

#set sample.meta.data row order to match samples

sample.meta.data$Sample_ID <- factor(sample.meta.data$Sample_ID, levels = samples_for_meta)
sample.meta.data <- sample.meta.data[order(sample.meta.data$Sample_ID), ]

meta <- seurat@meta.data

#lift sample-level metadata over to each cell (row in meta.data)
for (i in seq_along(meta.name)) {
  
  for (j in seq_along(samples_for_meta)) {
    
    meta[meta$Sample_ID %in% samples_for_meta[j], meta.name[i]] <- sample.meta.data[j,i]
    
  }
  
  #set factor levels based on order or sample-level metadata entries 
  meta[,meta.name[i]] <- factor(meta[,meta.name[i]], levels = unique(sample.meta.data[,i]))
  
}


seurat <- AddMetaData(seurat, meta)

StopWatchEnd$AddMetaData <- Sys.time()


                       
#########################################
# Add cell-level quality control metrics 
#########################################

print("")
print("*********************")
print("Add QC metrics")
print(Sys.time())
print("********************")
print("")
StopWatchStart$AddQC <- Sys.time()


#Add mito metadata
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

#Fix nUMI/nGene
seurat[["nGene_RNA"]] <- seurat[["nFeature_RNA"]]
seurat[["nUMI_RNA"]] <- seurat[["nCount_RNA"]]              
seurat[["nProtein_ADT"]] <- seurat[["nFeature_ADT"]]
seurat[["nUMI_ADT"]] <- seurat[["nCount_ADT"]]   

StopWatchEnd$AddQC <- Sys.time()

                       
                       
#########################################
# Add QC flag
#########################################

print("")
print("*********************")
print("Add QC flag")
print(Sys.time())
print("********************")
print("")
StopWatchStart$QCFlag <- Sys.time()
                       
                       
###Flag outliers of percent.mt 

mitoQC <- WhichCells(seurat, expression = percent.mt > maxMitoPerCell)
seurat@meta.data[rownames(seurat@meta.data) %in% mitoQC, "mitoQC"] <- "FAIL"
seurat@meta.data[is.na(seurat@meta.data$mitoQC), "mitoQC"] <- "PASS"


###Flag outliers of nGene

geneQC <- WhichCells(seurat, expression = nGene_RNA < minGenesPerCell)
seurat@meta.data[rownames(seurat@meta.data) %in% geneQC, "geneQC"] <- "FAIL"
seurat@meta.data[is.na(seurat@meta.data$geneQC), "geneQC"] <- "PASS"


###Flag outliers of nUMI

#umiQC <- WhichCells(seurat, expression = nUMI < 1000)
#seurat@meta.data[rownames(seurat@meta.data) %in% umiQC, "umiQC"] <- "FAIL"
#seurat@meta.data[is.na(seurat@meta.data$umiQC), "umiQC"] <- "PASS"

###Create Cell-level QC flag (Cell_QC)

a <- paste0(seurat@meta.data$nGene_QC,
            seurat@meta.data$mitoQC)
seurat@meta.data$Cell_QC <- ifelse(grepl("FAIL",a),
                                       "FAIL",
                                       "PASS"
                                      )                       

StopWatchEnd$QCFlag <- Sys.time()


#########################################
# Read VDJ data and add to seurat object
#########################################


clono_metaTMP_TCR <- list()
clono_metaTMP_BCR <- list()

## ---------------------------- Load and append TCR data

for (i in seq_along(sampleNames)) {
  
  matrix_dir = dataDirTCR[i]
  
  #read in all_contig_annotations.csv from cellranger output folder
  tmp <- read.csv(file = paste0(matrix_dir, "/filtered_contig_annotations.csv"))

  #remove barcoes with non-productive TCRs - working off filtered_contig_annotations so all entries are production 
  #productive <- tmp[tmp$productive == "True", ]
  
  #subset csv to only include barcode name and raw_clonotype_id
  productive <- tmp[, c("barcode", "raw_clonotype_id")]
  
  #create another column that combines barcode and raw_clonotype_id so you can ultimately have a unique clonotype for each barcode
  productive$unique_barcode_clono <- paste(productive$barcode, 
                                           productive$raw_clonotype_id,
                                           sep = "_")
  
  unique <- unique(productive$unique_barcode_clono)
  
  productive_unique <- data.frame(t(data.frame(strsplit(unique, "_")  )))
  colnames(productive_unique) <- c("barcodes", "raw_clonotype_id")
  rownames(productive_unique) <- productive_unique[,"barcodes"]
  #productive_unique$barcodes <- NULL
  
  #paste sample name before barcode and remove -1 suffix so it matches seurat dge
  
  productive_unique$barcodes <- paste(sampleNames[i], 
                                       rownames(productive_unique), 
                                       sep = "_")
  productive_unique$barcodes <- gsub("-1$", "", productive_unique$barcodes)
  
  #subset productive_unique dataframe so it only includes barcodes in seurat dge
  barcodesWithSCdataTMP <- productive_unique[productive_unique$barcodes %in%
                           rownames(seurat@meta.data), ]
  
  rownames(barcodesWithSCdataTMP) <- barcodesWithSCdataTMP$barcodes
  barcodesWithSCdataTMP$raw_clonotype_id <- paste(sampleNames[i], 
                                       barcodesWithSCdataTMP$raw_clonotype_id, 
                                       sep = "_")
  
  #create a clonotype metadata list

  clono_metaTMP_TCR[[i]] <- read.csv(file = paste0(matrix_dir, "/clonotypes.csv"))
  clono_metaTMP_TCR[[i]]$clonotype_id <- paste(sampleNames[i], 
                                           clono_metaTMP_TCR[[i]]$clonotype_id, 
                                        sep = "_")
  
  
  clonotype.file <- paste0("./data/", sampleNames[i], "_clono_meta_tcr.csv")
  write.csv(clono_metaTMP_TCR,  file = clonotype.file)
  
  
  #Add CDR3_aa data to seurat object meta data 
  
  subset <- subset(seurat, subset = Sample_ID == sampleNames[i])

  barcodesWithSCdata_TCR <- barcodesWithSCdataTMP
  barcodesWithSCdata_TCR$barcodes <- NULL
  colnames(barcodesWithSCdata_TCR) <- "raw_tcr_clonotype_id"

  subset <- AddMetaData(subset, barcodesWithSCdata_TCR, col.name = "raw_tcr_clonotype_id")
  
  #Add CDR3_aa data to seurat object meta data 
  
  for (b in unique(clono_metaTMP_TCR[[1]]$clonotype_id)) {
  
  cdr3 <- as.character(clono_metaTMP_TCR[[1]][clono_metaTMP_TCR[[1]]$clonotype_id == b, "cdr3s_aa"])
  subset@meta.data[subset@meta.data$raw_tcr_clonotype_id %in% b, "cdr3_aa_tcr"] <- cdr3

  }

 
  ## ---------------------------- Load and append BCR data
  
  matrix_dir = dataDirBCR[i]
  
  #read in all_contig_annotations.csv from cellranger output folder
  tmp <- read.csv(file = paste0(matrix_dir, "/filtered_contig_annotations.csv"))

  #remove barcoes with non-productive BCRs - working off filtered_contig_annotations so all entries are production 
  #productive <- tmp[tmp$productive == "True", ]
  
  #subset csv to only include barcode name and raw_clonotype_id
  productive <- tmp[, c("barcode", "raw_clonotype_id")]
  
  #create another column that combines barcode and raw_clonotype_id so you can ultimately have a unique clonotype for each barcode
  productive$unique_barcode_clono <- paste(productive$barcode, 
                                           productive$raw_clonotype_id,
                                           sep = "_")
  
  unique <- unique(productive$unique_barcode_clono)
  
  productive_unique <- data.frame(t(data.frame(strsplit(unique, "_")  )))
  colnames(productive_unique) <- c("barcodes", "raw_clonotype_id")
  rownames(productive_unique) <- productive_unique[,"barcodes"]
  #productive_unique$barcodes <- NULL
  
  #paste sample name before barcode and remove -1 suffix so it matches seurat dge
  
  productive_unique$barcodes <- paste(sampleNames[i], 
                                       rownames(productive_unique), 
                                       sep = "_")
  productive_unique$barcodes <- gsub("-1$", "", productive_unique$barcodes)

  
  #subset productive_unique dataframe so it only includes barcodes in seurat dge
  barcodesWithSCdataTMP <- productive_unique[productive_unique$barcodes %in%
                           rownames(seurat@meta.data), ]
  
  rownames(barcodesWithSCdataTMP) <- barcodesWithSCdataTMP$barcodes
  barcodesWithSCdataTMP$raw_clonotype_id <- paste(sampleNames[i], 
                                       barcodesWithSCdataTMP$raw_clonotype_id, 
                                       sep = "_")
  
  #create a clonotype metadata list

  clono_metaTMP_BCR[[i]] <- read.csv(file = paste0(matrix_dir, "/clonotypes.csv"))
  clono_metaTMP_BCR[[i]]$clonotype_id <- paste(sampleNames[i], 
                                           clono_metaTMP_BCR[[i]]$clonotype_id, 
                                        sep = "_")
  
  
  clonotype.file <- paste0("./data/", sampleNames[i], "_clono_meta_bcr.csv")
  write.csv(clono_metaTMP_BCR,  file = clonotype.file)
  
  
  #Add CDR3_aa data to seurat object meta data 
  
  barcodesWithSCdata_BCR <- barcodesWithSCdataTMP
  barcodesWithSCdata_BCR$barcodes <- NULL
  colnames(barcodesWithSCdata_BCR) <- "raw_ig_clonotype_id"

  subset <- AddMetaData(subset, barcodesWithSCdata_BCR, col.name = "raw_ig_clonotype_id")
  
  #Add CDR3_aa data to seurat object meta data 
  
  for (c in unique(clono_metaTMP_BCR[[1]]$clonotype_id)) {
  
  cdr3 <- as.character(clono_metaTMP_BCR[[1]][clono_metaTMP_BCR[[1]]$clonotype_id == c, "cdr3s_aa"])
  subset@meta.data[subset@meta.data$raw_ig_clonotype_id %in% c, "cdr3_aa_bcr"] <- cdr3

  }

   
  #Save sample-specific seurat object

  seurat.file <- paste0("./data/", sampleNames[i], "_preprocessed_seurat.rds")
  saveRDS(subset, file = seurat.file)

}




StopWatchEnd$SaveData <- Sys.time()
StopWatchEnd$Overall <- Sys.time()


#########################################
# Stop watches
#########################################

print("")
print("********************")
print("Run Time")
print(Sys.time())
print("********************")
print("")

clocks <- mapply('-', StopWatchEnd, StopWatchStart, SIMPLIFY = FALSE)
print(clocks)
#bb <- as.numeric(bb, units = "mins")
#lapply(bb, function(x){as.numeric(x, units = "mins")})


#########################################
# Session Information
#########################################

print("")
print("********************")
print("Session Info")
print(Sys.time())
print("********************")
print("")

print(sessionInfo())

print("***************************")
print("*******END OF SCRIPT*******")
print("***************************")
