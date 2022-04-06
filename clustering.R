############################################################################
#   Script to run multi-resolution clustering pipeline on multiOMICs data  #
#     L.Richards (multi-res clustering) D.Croucher (mods/qc/multiomics)    #
#                             April 2022                                   #
############################################################################

###########################################################
### GENERAL OVERVIEW OF THIS SCRIPT
### Load preprocessed scRNA-seq data (seurat objects)
### Remove flagged cells (based on QC cutoffs from pre-processing script)
### Normalize gene and protein UMIs 
### Scale gene expression matrix
### Identify variable genes (vst)
### Cell cyle scoring
### Run PCA and determine number of significant PCs (can auto-detect scree)
### Non-linear dimensionality reduction (tSNE, UMAP)
### Cluster cells over range of 6 resolutions
### Run differential gene expression (DGE) over all clustering resolutions
### Calculate silhouette width for clusters
### Select optimal clustering resolution based on DE-gene cutoff and max silhouette width
### Save data
### Output plots
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
## Rscript /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/scripts/multiresclustering.R 
##     --seuratObj /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/analysis/preprocess/ \
##     --outputDir /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/analysis/clustered/ \
##     --samples /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/analysis/ALQ_MCRN007_BCMA_sampleNames.txt \
##     --applyQCfilter TRUE \
##     --gene.filtering 0.001 \
##     --normalizationMethod LogNormalize \
##     --varsRegress /cluster/projects/pughlab/projects/scVKMYC/analysis/cohort/varsRegress.txt \
##     --numVarGenes 3000 \
##     --pcaScalingGenes var \
##     --computePC 75 \
##     --numPC 0 \
##     --addPC 0 \
##     --outputPlots TRUE \
##     --minResolution 0.4 \
##     --maxResolution 1.4 \
##     --optimalCluster FALSE \
##     --selectedResolution 0.4 \
##     --deGeneCutoff 10 \
##     --fdrCutoff 0.05 \
##     --markerGenes /cluster/projects/pughlab/projects/ALQ_MCRN007_BCMA/scomics/analysis/genelists/CellTypeMarkerGenes_B_PC_Cells.csv \
##     --species human
###########################################################

### start clocks
StopWatchStart <- list()
StopWatchEnd   <- list()
StopWatchStart$Overall  <- Sys.time()
StopWatchEnd$Overall  <- Sys.time()

#########################################
# 1) Load packages
#########################################

print("")
print("********************")
print("Load packages")
print(Sys.time())
print("********************")
print("")
StopWatchStart$LoadPackages  <- Sys.time()

#suppressMessages(library(scran))
suppressMessages(library(Matrix))
suppressMessages(library(scales))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(DropletUtils))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggExtra))
suppressMessages(library(ggpubr))
suppressMessages(library(sctransform))
suppressMessages(library(future))
suppressMessages(library(KneeArrower))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(ggrepel))
suppressMessages(library(patchwork))
suppressMessages(library(devtools))
suppressMessages(library(plyr))
suppressMessages(library(readxl))
suppressMessages(library(stringr))
suppressMessages(library(Hmisc))
suppressMessages(library(grid))
suppressMessages(library(reshape))
suppressMessages(library(reshape2))
suppressMessages(library(gplots))

StopWatchEnd$LoadPackages  <- Sys.time()


#########################################
#  Parse Options
#########################################

option_list <- list(make_option("--seuratObj",
                                type = "character",
                                default = NULL,
                                help = "path to seurat object to use as input (post-QC doublets)",
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
                    make_option("--applyQCfilter",
                                type = "logical",
                                help = "TRUE/FALSE to remove cells with FAIL QC flag",
                                default = NULL,
                                metavar= "logical"
                               ),
                    make_option("--gene.filtering",
                                type = "double",
                                help = "percent to be removed fro mavg library size ie. 0.005",
                                default = NULL,
                                metavar= "double"
                               ),
                    make_option("--normalizationMethod",
                                type = "character",
                                help = "LogNormalize, scran or scTransform",
                                default = NULL,
                                metavar= "character"
                               ),
                     make_option("--varsRegress",
                                type = "character",
                                default = NULL,
                                help = "path to csv file with header REGRESS; or list NONE",
                                metavar= "character"
                               ),
                    make_option("--numVarGenes",
                                type = "integer",
                                default = NULL,
                                help = "number of varibale genes to identify",
                                metavar= "integer"
                               ),
                    make_option("--pcaScalingGenes",
                                type = "character",
                                help = "Use all or var for genes used in data scaling and PCA",
                                default = NULL,
                                metavar= "character"
                               ),
                    make_option("--computePC",
                                type = "integer",
                                help = "how many PCs to compute in PCA",
                                default = NULL,
                                metavar= "integer"
                               ),
                    make_option("--numPC",
                                type = "integer",
                                help = "specify how many PCs to use for dim rediction, or list 0 for autodetection",
                                default = NULL,
                                metavar= "integer"
                               ),
                    make_option("--addPC",
                                type = "integer",
                                help = "how many PCs to add to numPC to capture additional biological variation",
                                default = NULL,
                                metavar= "integer"
                               ),
                    make_option("--outputPlots",
                                type = "logical",
                                help = "TRUE/FALSE to output plots",
                                default = NULL,
                                metavar= "logical"
                               ),
                    make_option("--minResolution",
                                type = "double",
                                help = "minimum clustering resolution; ie. 0.5",
                                default = NULL,
                                metavar= "double"
                               ),
                    make_option("--maxResolution",
                                type = "double",
                                default = NULL,
                                help = "maximum clustering resolution; ie. 1.0",
                                metavar= "double"
                               ),
                    make_option("--deGeneCutoff",
                                type = "integer",
                                default = NULL,
                                help = "cutoff for number DE genes per cluster; ie. 50",
                                metavar= "integer"
                               ),
                    make_option("--fdrCutoff",
                                type = "double",
                                default = NULL,
                                help = "Filter DE genes based on this FDR threshold, ir. 0.01",
                                metavar= "double"
                               ),
                   make_option("--optimalCluster",
                                type = "logical",
                                default = TRUE,
                                help = "TRUE/FALSE to select optimal clustering solution?",
                                metavar= "charater"
                                ),
                    make_option("--selectedResolution",
                                type = "integer",
                                default = NULL,
                                help = "selected clustering resolution; ie. 1.0",
                                metavar= "integer"
                               ),
                   make_option("--markerGenes",
                                type = "character",
                                default = NULL,
                                help = "Path to csv with marker genes, or use NONE",
                                metavar= "charater"
                               ),
                    make_option("--species",
                                type = "character",
                                default = NULL,
                                help = "species used in expeirment, current options include mouse or human)",
                                metavar= "charater"
                               )
)



opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

seuratObj <- opt$seuratObj
outputDir <- opt$outputDir
samples <- opt$samples
applyQCfilter <- opt$applyQCfilter
gene.filtering <- opt$gene.filtering
normalizationMethod <- opt$normalizationMethod
varsRegress <- opt$varsRegress
numVarGenes <- opt$numVarGenes
pcaScalingGenes <- opt$pcaScalingGenes
computePC <- opt$computePC
numPC  <- opt$numPC
addPC <- opt$addPC
outputPlots <- opt$outputPlots
minResolution <- opt$minResolution
maxResolution <- opt$maxResolution
deGeneCutoff <- opt$deGeneCutoff
fdrCutoff <- opt$fdrCutoff
optimalCluster <- opt$optimalCluster
selectedResolution <- opt$selectedResolution
markerGenes <- opt$markerGenes
species <- opt$species



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


seuratDirs <- c(paste0(seuratObj,
                       sampleNames, 
                       "/data/",
                       sampleNames,
                       "_preprocessed_seurat.rds"))

print(paste0("Path to seurat object: ", seuratDirs))


if (varsRegress == "NONE"){

    vars.regress <- NULL
    print(paste0("Variables to Regress: ", vars.regress))

} else {

    vars.regress <-  c(as.character(read.table(varsRegress, header = T)$REGRESS))
    print(paste0("Variables to Regress: ", paste(vars.regress, collapse = ", ")))

}


print(paste0("Output directory: ", outputDir))
print(paste0("Lowly Expressed Gene Cutoff: ", gene.filtering))
print(paste0("Normalization Method: ", normalizationMethod))
print(paste0("Number Variable Genes: ", numVarGenes))
print(paste0("Genes used for PCA and scaling: ", pcaScalingGenes))
print(paste0("Number of PCs to compute: ", computePC))
print(paste0("Number of PCs to use: ", ifelse(numPC == 0, "Auto-detect", numPC)))
print(paste0("Number of PCs to add: ", ifelse(normalizationMethod == "scTransform", "1.5x", addPC)))
print(paste0("Output Plots?: ", outputPlots))
print(paste0("Minimum Clustering Resolution: ", minResolution))
print(paste0("Maximum Clustering Resolution: ", maxResolution))
print(paste0("FDR Cutoff for Marker Genes: ", fdrCutoff))
print(paste0("Clusters must have a minimum of: ", deGeneCutoff, " DE genes"))
print(paste0("Run optimal clustering solution?: ", optimalCluster))

if (markerGenes == "NONE"){

    marker.genes <- NULL
    print(paste0("Assessing marker genes: ", markerGenes))

} else {

    marker.genes <- read.csv(markerGenes)
    print(paste0("Assessing marker genes: ", paste(marker.genes$GENE, collapse = ",")))

}

StopWatchEnd$SetUpVaribles  <- Sys.time()

for (name in sampleNames) {

 setwd(outputDir)
 dir.create(name)
 setwd(name)
 dir.create("data")
 if (outputPlots == TRUE){dir.create("figures")}
  
}


#########################################
# Load preprocessed scRNA-seq data (seurat objects)
#########################################

print("")
print("*********************")
print("Read in Data")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ReadData <- Sys.time()


seurat.obj <- readRDS(seuratDirs)

print("Total Gene Expression Dataset Size")
print(dim(seurat.obj@assays$RNA@data)) ## print out final size of object

print("Total Protein Expression Dataset Size")
print(dim(seurat.obj@assays$ADT@data)) ## print out final size of object

StopWatchEnd$ReadData <- Sys.time()


#########################################
# Remove flagged cells (based on QC cutoffs from pre-processing script)
#########################################

print("")
print("********************")
print("Remove Failed QC Cells")
print(Sys.time())
print("********************")
print("")
StopWatchStart$RemoveFlaggedCells <- Sys.time()

print(table(seurat.obj@meta.data$Cell_QC))
seurat.obj <- subset(seurat.obj , subset = Cell_QC == "PASS")
print(dim(seurat.obj@assays$RNA@data))

StopWatchEnd$RemoveFlaggedCells <- Sys.time()



#########################################
# Remove lowly expressed genes
#########################################

print("")
print("********************")
print("Filter out lowly expressed genes")
print(Sys.time())
print("********************")
print("")
StopWatchStart$FilterGenes <- Sys.time()

####remove gene.filtering % of lowest library size

nCells <- mean((table(seurat.obj@meta.data$Sample_ID))) * gene.filtering
print(paste0("Cutoff: ", gene.filtering, " of average cell size (", round(nCells, 2), " cells)"))

tmp <- seurat.obj
seurat.obj <- CreateSeuratObject(counts = tmp@assays$RNA@counts,
                                 meta.data = tmp@meta.data,
                                 min.cells = nCells
                                )
seurat.obj[['ADT']] <- CreateAssayObject(counts = tmp@assays$ADT@counts)

print(dim(seurat.obj@assays$RNA@counts))
print(dim(seurat.obj@assays$RNA@data))

print(dim(seurat.obj@assays$ADT@counts))
print(dim(seurat.obj@assays$ADT@data))

StopWatchEnd$FilterGenes <- Sys.time()


#########################################
# Normalization
#########################################

StopWatchStart$Normalization <- Sys.time()


print("")
print("********************")
print("Normalize Data")
print(Sys.time())
print("********************")
print("")

seurat.obj <- NormalizeData(object = seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000
                           )
seurat.obj <- NormalizeData(seurat.obj, assay = "ADT", normalization.method = 'CLR', margin = 2)


StopWatchEnd$Normalization <- Sys.time()


#########################################
# VarGenes and Data Scaling
#########################################

if (normalizationMethod == "LogNormalize" | normalizationMethod == "scran"){

    print("")
    print("********************")
    print("Identify Variable Genes (vst)")
    print(Sys.time())
    print("********************")
    print("")
    StopWatchStart$VarGenes <- Sys.time()
    #use default Seurat Method, vst
    print(paste0("Calculating fixed number of variable genes... ", numVarGenes))
    seurat.obj <- FindVariableFeatures(seurat.obj,
                                       selection.method = "vst",
                                       nfeatures = numVarGenes
                                      )
    StopWatchEnd$VarGenes <- Sys.time()

    print("")
    print("********************")
    print("Scale Data")
    print(Sys.time())
    print("********************")
    print("")
    StopWatchStart$ScaleData <- Sys.time()
    print(paste0("Regressing out....", paste(vars.regress, collapse = ", ")))
    print(paste0("Genes included in scaling... ", pcaScalingGenes))


    if (pcaScalingGenes == "all"){

        all.genes <- rownames(seurat.obj)
        seurat.obj <- ScaleData(seurat.obj,
                                vars.to.regress = vars.regress,
                                features = all.genes
                               )

    } else if (pcaScalingGenes == "var"){

            seurat.obj <- ScaleData(seurat.obj,
                                vars.to.regress = vars.regress,
                                features =  VariableFeatures(object = seurat.obj)
                               )
    }


}

StopWatchEnd$ScaleData <- Sys.time()


#########################################
# Cell Cycle Scoring
#########################################

print("")
print("********************")
print("Score Cells for Cell Cycle Markers")
print(Sys.time())
print("********************")
print("")
StopWatchStart$ScoreCellCycle <- Sys.time()
#As of now, there is no option to regress out cell cycle in this pipeline

# Also read in a list of cell cycle markers, from Tirosh et al, 2015

if (species == "human") {
 
  cc.genes <- readLines(con = "/cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/regev_lab_cell_cycle_genes.txt")
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]

}

if (species == "mouse") {
 
  cc.genes.human <- readLines(con = "/cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/regev_lab_cell_cycle_genes.txt")
  s.genes.human <- cc.genes.human[1:43]
  g2m.genes.human <- cc.genes.human[44:97]
  
  cc.genes <- read.table("/cluster/projects/pughlab/projects/scVKMYC/analysis/genelists/regev_lab_cell_cycle_genes_both.txt")
  
  s.genes <- cc.genes[cc.genes$HGNC.symbol %in% s.genes.human, ]
  s.genes <- s.genes$MGI.symbol

  g2m.genes <- cc.genes[cc.genes$HGNC.symbol %in% g2m.genes.human, ]
  g2m.genes <- g2m.genes$MGI.symbol

  
}

seurat.obj <- CellCycleScoring(seurat.obj,
                               s.features = s.genes,
                               g2m.features = g2m.genes,
                               set.ident = FALSE
                              )
#calculate CC.Differennc
#just in case you want to regress with alt method in future
seurat.obj$CC.Difference <- seurat.obj$S.Score - seurat.obj$G2M.Score

StopWatchEnd$ScoreCellCycle <- Sys.time()




#########################################
# Run PCA
#########################################

print("")
print("********************")
print("Run PCA")
print(Sys.time())
print("********************")
print("")
StopWatchStart$PCA <- Sys.time()

#print(paste0("Calculating... ", computePC, " PCs"))
print(paste0("Calculating... ", 100, " PCs"))
print(paste0("Genes included in PCA... ", pcaScalingGenes))

if (pcaScalingGenes == "all"){

    all.genes <- rownames(seurat.obj)
    seurat.obj <- RunPCA(seurat.obj,
                         features = all.genes,
                         #npcs = computePC,
                         npcs = 100,
                         verbose = FALSE
                        )

} else if (pcaScalingGenes == "var"){

    seurat.obj <- RunPCA(seurat.obj,
                         features = VariableFeatures(object = seurat.obj),
                         npcs = computePC,
                         #npcs = 100,
                         verbose = FALSE
                        )
}


#########################################
# Determine significant PCs
#########################################


### isolate PCA data for output
pca <- data.frame(seurat.obj@reductions$pca@stdev)
pca$PC <- seq(1:nrow(pca))
pca <- pca[ ,c(2,1)]
colnames(pca)[2] <- "st.dev"
eigs <- pca$st.dev**2
pca$prop.var <- eigs / sum(eigs)


if(numPC == 0){

    print(paste0("Determining significant PCs...Auto (KneeArrower)"))


    cutoff.point <- findCutoff(pca$PC[1:computePC],
                               pca$st.dev[1:computePC],
                               method="first",
                               0.01 #derivative cutoff
                              )
    numPC <- round(cutoff.point$x)
    print(paste0(numPC, " significant PCs"))

} else if(numPC > 0){

    print(paste0("Determining significant PCs...User Defined"))
    print(paste0(numPC, " significant PCs"))
}


if(addPC > 0){

    print(paste0("Adding Additional PCs...", addPC, " PCs"))
    pc.use <- numPC + addPC
    print(paste0("PCs used for downstream analysis...", pc.use, " PCs"))

} else if (addPC == 0){

    pc.use <- numPC
    print(paste0("PCs used for downstream analysis...", pc.use, " PCs"))

}


if(normalizationMethod == "scTransform"){

    #add 1.5 times PCs for scTransform to capture more biology
    #Seurat recommendation
    pc.use <- round(numPC * 1.5)
    print(paste0("Multiply significant PCs 1.5x for scTransform..."))
    print(paste0("PCs used for downstream analysis...", pc.use, " PCs"))

}

StopWatchEnd$PCA <- Sys.time()





#########################################
# tSNE and UMAP 
#########################################

print("")
print("********************")
print("Non-linear Dimensionality Reduction")
print(Sys.time())
print("********************")
print("")

StopWatchStart$tSNE <- Sys.time()
print(paste0("Running tSNE..."))
seurat.obj <- RunTSNE(seurat.obj, dims = 1:pc.use, verbose = FALSE, max_iter = 2000, check_duplicates = FALSE)
StopWatchEnd$tSNE <- Sys.time()

StopWatchStart$UMAP <- Sys.time()
print(paste0("Running UMAP..."))
seurat.obj <- RunUMAP(seurat.obj, dims = 1:pc.use, verbose = FALSE)
StopWatchEnd$UMAP <- Sys.time()

print(paste0("Add coordinates to metadata..."))
seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$umap@cell.embeddings)
                          )
seurat.obj <- AddMetaData(seurat.obj,
                          metadata = data.frame(seurat.obj@reductions$tsne@cell.embeddings)
                          )


#########################################
# Clustering
#########################################

print("")
print("********************")
print("Cluster Cells with Seurat")
print(Sys.time())
print("********************")
print("")
StopWatchStart$Clustering <- Sys.time()

# remove clustering results in meta.data from previous analyses
previous.clus.results <- colnames(seurat.obj@meta.data)[grep("_res", colnames(seurat.obj@meta.data))]
seurat.obj@meta.data <- seurat.obj@meta.data[, !colnames(seurat.obj@meta.data) %in% previous.clus.results]


#default in seurat is to use k = 20
print(paste0("Constructing SNN Graph with ", pc.use, " PCs..."))
seurat.obj <- FindNeighbors(seurat.obj,
                            dims = 1:pc.use,
                            verbose = TRUE
                           )

print(paste0("Find Clusters with ", pc.use, " PCs..."))


seurat.obj <- FindClusters(seurat.obj,
                           resolution = selectedResolution, 
                           #method = "igraph", #better for large data
                           algorithm = 1, #default Louvian algorithm,
                           verbose = FALSE,
                           n.start = 50
                          )

seurat.obj@meta.data$Optimal_res <- selectedResolution

##fix naming in meta.data

if (normalizationMethod == "LogNormalize" | normalizationMethod == "scran"){

    colnames(seurat.obj@meta.data) <- gsub("RNA_snn", "Seurat_cluster", colnames(seurat.obj@meta.data))
    seurat.obj@meta.data$seurat_clusters <- NULL

} else if (normalizationMethod == "scTransform"){

    colnames(seurat.obj@meta.data) <- gsub("SCT_snn", "Seurat_cluster", colnames(seurat.obj@meta.data))
    seurat.obj@meta.data$seurat_clusters <- NULL

}


StopWatchEnd$Clustering <- Sys.time()




#########################################
# Save Data
#########################################

print("")
print("********************")
print("Save Data")
print(Sys.time())
print("********************")
print("")
StopWatchStart$SaveData <- Sys.time()

#print("Raw Counts....")
#raw.count.dir <- paste0("./data/", fileName, "_rawCounts")
#raw <- Matrix(seurat.obj@assays$RNA@counts, sparse = TRUE)
#write10xCounts(path = raw.count.dir,
#               x = raw
#               )

#print("Normalized Counts....")
#norm.count.dir <- paste0("./data/", fileName, "_normCounts")
#norm.count <- Matrix(seurat.obj@assays$RNA@data, sparse = TRUE)
#write10xCounts(path = norm.count.dir,
#               x = norm.count
#               )

print("Meta Data....")
meta <- data.frame(seurat.obj@meta.data)
meta.file <- paste0("./data/", sampleNames, "_metaData.csv")
write.csv(meta, file = meta.file)

#print("Scree Plot Information....")
#scree.file <- paste0("./data/", fileName, "_scree.csv")
#write.csv(pca, file = scree.file)

#print("PCA....")
#pca.file <- paste0("./data/", fileName, "_PCA.rds")
#pca.mat <- seurat.obj@reductions$pca
#saveRDS(pca.mat, file = pca.file)

#if (optimalCluster == TRUE){
#print("DE Markers....")
#DE.file <- paste0("./data/", fileName, "_DE_markers.csv")
#write.csv(markers, file = DE.file)
#DE.file.2 <- paste0("./data/", fileName, "_DE_markers_FDR_filered.csv")
#write.csv(markers_fdr, file = DE.file.2)


#print("Silhouette Width....")
#sil.file <- paste0("./data/", fileName, "_silhouette_widths.rds")
#saveRDS(silhouette.width, file = sil.file)
#sil.file.2 <- paste0("./data/", fileName, "_clustering_dashboard.csv")
#write.csv(deGenes, file = sil.file.2)
#}

print("Saving Seurat Object....")
seurat.file <- paste0("./data/", sampleNames, "_seurat.rds")
saveRDS(seurat.obj, file = seurat.file)

StopWatchEnd$SaveData <- Sys.time()



#########################################
# Plotting
#########################################

if(outputPlots == TRUE){

    #############
    print("")
    print("********************")
    print("Output Plots")
    print(Sys.time())
    print("********************")
    print("")
    StopWatchStart$Plotting <- Sys.time()


    if(class(marker.genes) == "data.frame"){

        print("Adding Marker Genes to Metadata...")
        
        seurat.obj <- AddMetaData(seurat.obj, FetchData(seurat.obj, marker.genes$GENE))
        
    }

    meta <- data.frame(seurat.obj@meta.data)

    #######################

    print("1/6....QC plots")

    QC1 <- ggplot(meta, aes(x=nGene_RNA, y=nUMI_RNA)) +
          geom_point(size = 0.6, alpha = 0.7) +
          theme(legend.position="none") + theme_classic()
    QC1 <- ggMarginal(QC1, type="histogram")

    QC2 <- ggplot(meta, aes(x=nGene_RNA, y=percent.mt)) +
          geom_point(size = 0.6, alpha = 0.7) +
          theme(legend.position="none") + theme_classic()
    QC2 <- ggMarginal(QC2, type="histogram")


    qc.plots <- paste0("./figures/",sampleNames, "_QC.pdf")
    print(qc.plots)
    pdf(qc.plots, width = 10, height = 5)
    print(ggarrange(QC1, QC2, ncol = 2, nrow = 1))
    dev.off()

    #######################

    print("2/6....Variable Gene Plot")

    top20 <- head(VariableFeatures(seurat.obj), 20)

    VAR1 <- VariableFeaturePlot(seurat.obj)
    VAR2 <- LabelPoints(plot = VAR1, points = top20, repel = TRUE)

    var.pdf <- paste0("./figures/", sampleNames, "_VarGenes.pdf")
    print(var.pdf)
    pdf(var.pdf, width = 10, height = 8)
    print(VAR2)
    dev.off()

    #######################

    #print("3/6....PCA Scree Plot")

    #pc.pdf <- paste0("./figures/", sampleNames, "_PCA_plots.pdf")
    #print(pc.pdf)
    #pdf(pc.pdf, width = 10, height = 5)
    #par(mfrow=c(1,2))

     #       plot(pca$PC,
     #            pca$st.dev,
     #            xlab = "Principal Components",
     #            ylab = "Standard Deviation",
     #            main = fileName
     #               )
     #       abline(v=round(pc.use), lty = 2)
     #       legend("topright",
     #              legend = c(paste0(round(pc.use), " PCs")),
     #              bty = 'n')

#        plot(pca$PC,
#             pca$prop.var,
#             xlab = "Principal Components",
#             ylab = "Proportion Variance"
#            )
#        abline(v=round(pc.use), lty = 2)
#        legend("topright",
#               legend = c(paste0(round(sum(pca$prop.var[1:round(pc.use)]), 3), " cumulative variance")),
#               bty = 'n')
#
#    dev.off()


    #######################

    print("5/6.....tSNEs and UMAPs")


  #### PLOT tSNEs

    meta <- seurat.obj@meta.data
    meta$cluster <- meta[ ,colnames(meta) == paste0("Seurat_cluster_res.", unique(meta$Optimal_res))]
    #clust.col <- colnames(meta)[grep(unique(meta$cluster), colnames(meta))]
    cent <- meta %>% group_by(cluster) %>% select(tSNE_1, tSNE_2) %>% summarize_all(median)
    
  tsne.c <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "cluster")) +
            geom_point(alpha = 0.5) + theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            labs("Clusters") + ggtitle(paste0(sampleNames, " - Cluster: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) +
            geom_label_repel(aes(label = cluster), data = cent, show.legend = F, size = 2)  +
            theme(legend.position='none') 
  
#  tsne.s <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "Sample_ID")) +
#            geom_point(alpha = 0.5) + theme_classic() +
#            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
#            labs("Samples") + ggtitle(paste0(sampleNames, " - Sample ID: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) 
#            #geom_label_repel(aes(label = Sample_ID), data = cent, show.legend = F, size = 2)  +
#            #theme(legend.position='none') 
  
#  tsne.d <- ggplot(meta, aes_string(x="tSNE_1", y = "tSNE_2", col = "Cohort_ID")) +
#            geom_point(alpha = 0.5) + theme_classic() +
#            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
#            labs("Disease Group") + ggtitle(paste0(sampleNames, " - Disease Group: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) 
#            #geom_label_repel(aes(label = Cohort_ID), data = cent, show.legend = F, size = 2)  +
#            #theme(legend.position='none') 

  tsne.pdf <- paste0("./figures/", sampleNames, "_tSNE.pdf")
  pdf(tsne.pdf, width = 10, height = 10)
  print(tsne.c)
#  print(tsne.s)
#  print(tsne.d)
  dev.off()

  
  
  
  #### PLOT UMAPs
    
    cent <- meta %>% group_by(cluster) %>% select(UMAP_1, UMAP_2) %>% summarize_all(median)
    
  umap.c <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "cluster")) +
            geom_point(alpha = 0.5) + theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            labs("Clusters") + ggtitle(paste0(sampleNames, " - Cluster: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) +
            geom_label_repel(aes(label = cluster), data = cent, show.legend = F, size = 2)  +
            theme(legend.position='none') 
  
 # umap.s <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "Sample_ID")) +
 #           geom_point(alpha = 0.5) + theme_classic() +
 #           theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
 #           labs("Samples") + ggtitle(paste0(sampleNames, " - Sample ID: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) 
 #           #geom_label_repel(aes(label = Sample_ID), data = cent, show.legend = F, size = 2)  +
 #           #theme(legend.position='none') 
  
 # umap.d <- ggplot(meta, aes_string(x="UMAP_1", y = "UMAP_2", col = "Cohort_ID")) +
 #           geom_point(alpha = 0.5) + theme_classic() +
 #           theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
 #           labs("Disease Group") + ggtitle(paste0(sampleNames, " - Disease Group: res ", unique(meta$Optimal_res),  " (", nrow(meta), " cells)")) 
 #           #geom_label_repel(aes(label = Cohort_ID), data = cent, show.legend = F, size = 2)  +
 #           #theme(legend.position='none') 

  umap.pdf <- paste0("./figures/", sampleNames, "_UMAP.pdf")
  pdf(umap.pdf, width = 10, height = 10)
  print(umap.c)
  #print(umap.s)
  #print(umap.d)  
  dev.off()
  
  
  
    tsne.tiff <- paste0("./figures/", sampleNames, "_Cluster_tSNE.tiff")
    tiff(tsne.tiff, width = 700, height = 600)
    print(tsne.c)
    dev.off()
  
#    tsne.tiff <- paste0("./figures/", sampleNames, "_Sample_tSNE.tiff")
#    tiff(tsne.tiff, width = 700, height = 600)
#    print(tsne.s)
#    dev.off()

#    tsne.tiff <- paste0("./figures/", sampleNames, "_Cohort_tSNE.tiff")
#    tiff(tsne.tiff, width = 700, height = 600)
#    print(tsne.d)
#    dev.off()

  
    umap.tiff <- paste0("./figures/", sampleNames, "_Cluster_umap.tiff")
    tiff(umap.tiff, width = 700, height = 600)
    print(umap.c)
    dev.off()
  
 #   umap.tiff <- paste0("./figures/", sampleNames, "_Sample_umap.tiff")
 #   tiff(umap.tiff, width = 700, height = 600)
 #   print(umap.s)
 #   dev.off()

 #   umap.tiff <- paste0("./figures/", sampleNames, "_Cohort_umap.tiff")
 #   tiff(umap.tiff, width = 700, height = 600)
 #   print(umap.d)
 #   dev.off()
  
  
  #######################


    if(class(marker.genes) == "data.frame"){

        #######################

        print("6/6.....Marker Gene Expression")

        sigs <- colnames(meta)[colnames(meta) %in% marker.genes$GENE]
        print(sigs)

        #### UMAP with marker gene expression in hexbin
        markers.pdf <- paste0("./figures/", sampleNames, "_GeneMarkers_UMAP.pdf")
        print(markers.pdf)
        pdf(markers.pdf, height = 5, width = 5)

        for (i in 1:length(sigs)){

                plot.title <-  paste(sigs[i], "\n", sampleNames, "  ", nrow(meta), "cells")
                p <-  ggplot(meta, aes_string(x="UMAP_1", y="UMAP_2", z=sigs[i])) +
                        stat_summary_hex(bins=100, fun = "mean") +
                         theme_classic(base_size=8) +
                         scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                         theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                        ggtitle(plot.title) + labs(z = sigs[i])
                print(p)

            }

        dev.off()

       #### tSNE with marker gene expression in hexbin
      markers.pdf.t <- paste0("./figures/", sampleNames, "_GeneMarkers_tSNE.pdf")
      print(markers.pdf.t)
      pdf(markers.pdf.t, height = 5, width = 5)

        for (i in 1:length(sigs)){

                plot.title <-  paste(sigs[i], "\n", sampleNames, "  ", nrow(meta), "cells")
                p <-  ggplot(meta, aes_string(x="tSNE_1", y="tSNE_2", z=sigs[i])) +
                        stat_summary_hex(bins=100, fun = "mean") +
                         theme_classic(base_size=8) +
                         scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                         theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                        ggtitle(plot.title) + labs(z = sigs[i])
                print(p)

            }

    dev.off()

    StopWatchEnd$Plotting <- Sys.time()

    }

}


#### Tiff files of tSNE with marker gene expression in hexbin

setwd("figures")
dir.create("markers")

for (i in seq_along(sigs)) {
    
    goi <- sigs[i] 
    
    meta_ordered <- meta[order(meta[,goi], decreasing = F), ]

    plot.title <-  paste(goi, "\n", sampleNames, "  ", nrow(meta), "cells")
    
    p <-  ggplot(meta_ordered, aes_string(x="tSNE_1", y="tSNE_2", z=goi)) +
                        stat_summary_hex(bins=100, fun = "mean") +
                         theme_classic(base_size=8) +
                         scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                         theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                        ggtitle(plot.title) + labs(z = goi)
                print(p)
    
    mypath <- paste0("./markers/", sampleNames, "_", goi, "_GeneMarkers_tSNE.tiff")

    tiff(file=mypath)
    print(p)
    dev.off()
    
}


#### Tiff files of umap with marker gene expression in hexbin

for (i in seq_along(sigs)) {
    
    goi <- sigs[i] 
    
    meta_ordered <- meta[order(meta[,goi], decreasing = F), ]

    plot.title <-  paste(goi, "\n", sampleNames, "  ", nrow(meta), "cells")
    
    p <-  ggplot(meta_ordered, aes_string(x="UMAP_1", y="UMAP_2", z=goi)) +
                        stat_summary_hex(bins=100, fun = "mean") +
                         theme_classic(base_size=8) +
                         scale_fill_gradientn("Expression", colours = c(brewer.pal(n = 8, name = "YlOrRd")))+
                         theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                        ggtitle(plot.title) + labs(z = goi)
                print(p)
    
    mypath <- paste0("./markers/", sampleNames, "_", goi, "_GeneMarkers_UMAP.tiff")

    tiff(file=mypath)
    print(p)
    dev.off()
    
  }




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
