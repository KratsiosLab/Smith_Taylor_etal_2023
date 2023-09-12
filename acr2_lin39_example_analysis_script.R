# Sample code for processing young adult C. elegans 10X Genomics single-cell data from motor neurons

# Sample 1 reading in files, Emptydrops to distinguish cells from empty droplets, SoupX background 
# RNA correction, and QC.

library("DropletUtils")
library("scater")
library("dplyr")
library("wbData")
library("Matrix")
library("monocle3")

# loading a gene lookup table for converting between WBGeneIDs and gene symbols
gids <- wb_load_gene_ids("WS273")

setwd("./Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/")

# Loading in custom functions for plotting and other work
source("./Monocle3Functions.txt")

# Loading Sample 1
my_path <- "./3853-ST_new_alignment/3853-ST-1-selected/raw_feature_bc_matrix/"

wt.1 <- read10xCounts(my_path)

wt.1
head(colData(wt.1))
rownames(colData(wt.1)) <- colData(wt.1)$Barcode
colData(wt.1)$Sample <- "wt.1"
colData(wt.1)$stage <- "YA"
wt.1

head(rowData(wt.1))
rowData(wt.1) <- rowData(wt.1)[,1:2]
colnames(rowData(wt.1)) <- c("id", "gene_short_name")
rowData(wt.1)$gene_short_name <- i2s(rowData(wt.1)$id, gids)

# Distinguishing droplets with cells from empty droplets
bcrank <- barcodeRanks(counts(wt.1))

set.seed(100)
e.out <- emptyDrops(counts(wt.1), lower = 50)
summary(e.out$FDR <= 0.01)

# Mode     FALSE    TRUE    NA's 
# logical  3918    6180  682915 

# Looking at the top detected genes in the empty droplets
wt.1.ambient <- e.out@metadata$ambient
dim(wt.1.ambient)
wt.1.ambient[1:10,]
wt.1.ambient.df <- data.frame(row.names = rownames(wt.1.ambient),
                              id = rownames(wt.1.ambient),
                              gene_name = i2s(rownames(wt.1.ambient), gids),
                              expression = wt.1.ambient[,1])

wt.1.ambient.df <- wt.1.ambient.df %>% arrange(desc(expression))
head(wt.1.ambient.df, 15)

# rrn-3.1, mitochondrial genes (nduo-6, atp-6, MTCE.7, ctc-1, ctc-3),
# some neuronal genes (far-1, sbt-1, nlp-21), nspc-9, npsc-14, npsc-20

# Other diagnostics

is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm = TRUE)

table(Limited=e.out$Limited, Significant=is.cell)
#          Significant
# Limited  FALSE  TRUE
#   FALSE 3918  734
# TRUE      0 5446

# Statistics on total UMI counts for cells called as cells

summary(e.out[which(e.out$FDR <= 0.01),]$Total)
# Min.  1st Qu.   Median    Mean    3rd Qu.    Max. 
# 51.0   117.0   781.5  1181.0  1513.0 49659.0 

wt.1$PValue <- e.out$PValue
head(colData(wt.1))

wt.1$FDR <- e.out$FDR
head(colData(wt.1))

is.cell <- e.out$FDR <= 0.01
wt.1 <- wt.1[, which(is.cell), drop = F]

# Making a new filtered gene by barcode matrix for SoupX decontaminationn
wt.1.filt.counts <- counts(wt.1)

# Writing to a new folder for SoupX.
# Generating new filtered matrix files for SoupX
library(Matrix)
writeMM(wt.1.filt.counts, "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/filtered_feature_bc_matrix/matrix.mtx")
system("gzip ./3853-ST_new_alignment/3853-ST-1-selected/SoupX/filtered_feature_bc_matrix/matrix.mtx")
writeMM(wt.1.filt.counts, "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/filtered_feature_bc_matrix/matrix.mtx")
fil.barcodes <- colnames(wt.1)
write.table(fil.barcodes, "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system("gzip ./3853-ST_new_alignment/3853-ST-1-selected/SoupX/filtered_feature_bc_matrix/barcodes.tsv")
write.table(fil.barcodes, "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

sizeFactors(wt.1) <- librarySizeFactors(wt.1)
wt.1 <- logNormCounts(wt.1)

library(Seurat)
wt.1.S <- as.Seurat(wt.1)

# Generating initial clustering with Seurat to aid in SoupX background correction
wt.1.S <- NormalizeData(wt.1.S)
wt.1.S <- FindVariableFeatures(wt.1.S, selection.method = "vst", nfeatures = 4000)
wt.1.S <- ScaleData(wt.1.S)
wt.1.S <- RunPCA(wt.1.S, features = VariableFeatures(wt.1.S), npcs = 75)

ElbowPlot(wt.1.S, ndims = 75)
wt.1.S <- FindNeighbors(wt.1.S, dims = 1:75)
wt.1.S <- FindClusters(wt.1.S, resolution = 1.2)

wt.1.S <- RunUMAP(wt.1.S, dims = 1:75)
DimPlot(wt.1.S, reduction = "umap")

# Seurat identified 34 clusters. 
# This is the sample from the wt.1.

# Checking expression of several characteristic genes
FeaturePlot(wt.1.S, features = s2i("sbt-1", gids))
FeaturePlot(wt.1.S, features = s2i("unc-25", gids))
FeaturePlot(wt.1.S, features = s2i("pat-10", gids))
FeaturePlot(wt.1.S, features = s2i("unc-4", gids))
FeaturePlot(wt.1.S, features = s2i("myo-3", gids))
FeaturePlot(wt.1.S, features = s2i("bnc-1", gids))
FeaturePlot(wt.1.S, features = s2i("flp-13", gids))
FeaturePlot(wt.1.S, features = s2i("flp-1", gids))
FeaturePlot(wt.1.S, features = s2i("oig-1", gids))
FeaturePlot(wt.1.S, features = s2i("nspc-20", gids))

# Extracting reduced dimensions
class(wt.1.S@reductions$umap)
umap <- wt.1.S@reductions$umap

Embeddings(umap)[1:3,1:2]
wt.1_DR <- as.data.frame(Embeddings(umap))
wt.1_DR$Cluster <- Idents(wt.1.S)
head(wt.1_DR)

# Now running SoupX, with soupRange = c(0, 25)
library(SoupX)
my_dir_S <- "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/"
scl.wt.1 <- load10X(my_dir_S, soupRange = c(0, 25))

scl.wt.1$toc[1:5,1:5]
head(wt.1_DR)

head(scl.wt.1$soupProfile[order(scl.wt.1$soupProfile$est, decreasing = TRUE), ], n = 30)
scl.wt.1$soupProfile$gene_short_name <- i2s(rownames(scl.wt.1$soupProfile), gids)
head(scl.wt.1$soupProfile[order(scl.wt.1$soupProfile$est, decreasing = TRUE), ], n = 30)

# top genes are ribosomal RNA (rrn-3.1), mitochondrial genes (nduo-6, atp-6, ctc genes), 
# excretory gland (nspc-9, -20, -14, -19, -10, -13), some neuronal (nlp-21, flp-28)

# Looking the distribution of several possible contaminating genes
plotMarkerMap(scl.wt.1, s2i("pat-10", gids), wt.1_DR)
plotMarkerMap(scl.wt.1, s2i("mlc-3", gids), wt.1_DR)
scl.wt.1 <- setDR(scl.wt.1, wt.1_DR)
plotMarkerMap(scl.wt.1, s2i("flp-28", gids))
plotMarkerMap(scl.wt.1, s2i("nspc-20", gids))
plotMarkerMap(scl.wt.1, s2i("cup-4", gids))
plotMarkerMap(scl.wt.1, s2i("pdf-1", gids))
plotMarkerMap(scl.wt.1, s2i("clik-1", gids))
plotMarkerMap(scl.wt.1, s2i("nlp-12", gids))

# Using the following genes to estimate overall contamination
est.cont <- s2i(c("flp-28", "nspc-20", 'nspc-9', "nspc-14", "nspc-10", "nspc-13", "nspc-19", "pat-10", "clik-3", "cpn-3"), gids)

useToEst = estimateNonExpressingCells(scl.wt.1, nonExpressedGeneList = list(wt.1 = est.cont), 
                                      clusters = setNames(wt.1_DR$Cluster, rownames(wt.1_DR)))

plotMarkerMap(scl.wt.1, geneSet = est.cont, DR = wt.1_DR, useToEst = useToEst)
# Cells used to estimate the contamination

scl.wt.1 <- setClusters(scl.wt.1, wt.1_DR$Cluster)

scl.wt.1 <- calculateContaminationFraction(scl.wt.1, list(wt.1 = est.cont), useToEst = useToEst)
# Estimated global contamination fraction of 2.03%

# correcting the expressing, rounding to integers for downstream monocle3
system.time(out <- adjustCounts(scl.wt.1, roundToInt = TRUE))

cntSoggy = rowSums(scl.wt.1$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

tail(sort(rowSums(scl.wt.1$toc > out)/rowSums(scl.wt.1$toc > 0)), n = 20)

plotChangeMap(scl.wt.1, out, s2i("col-124", gids))
plotChangeMap(scl.wt.1, out, s2i("nspc-20", gids))
plotChangeMap(scl.wt.1, out, s2i("pat-10", gids))
plotChangeMap(scl.wt.1, out, s2i("flp-28", gids))
plotChangeMap(scl.wt.1, out, s2i("flp-1", gids))
plotChangeMap(scl.wt.1, out, s2i("flp-12", gids))

out[1:5,1:5]
scl.wt.1$toc[1:5,1:5]

saveRDS(scl.wt.1, "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/091622_soupX_object_wt.1.rds")
saveRDS(out, "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/091622_soupx_corrected_matrix_wt.1.rds")

dim(out)
head(colnames(out))
head(rownames(out))

head(colnames(wt.1.filt.counts))
head(rownames(wt.1.filt.counts))

wt.1_corrected <- SingleCellExperiment(assays = list(counts = out), colData = colData(wt.1))  

rowData(wt.1_corrected) <- rowData(wt.1)

# Using scater to calculate some quality control metrics

mt.genes <- c("WBGene00010959", "WBGene00010961", "WBGene00010966", "WBGene00010963", "WBGene00010957", "WBGene00010967", "WBGene00010958", "WBGene00010960", "WBGene00010962","WBGene00000829","WBGene00010965","WBGene00010964","WBGene00014454","WBGene00014472")
stress.genes <- c("WBGene00009692", "WBGene00009691", "WBGene00002018", "WBGene00002016", "WBGene00002017",
                  "WBGene00006959", "WBGene00002026", "WBGene00004622", "WBGene00235164", "WBGene00004131",
                  "WBGene00012813", "WBGene00000097", "WBGene00015168", "WBGene00009180", "WBGene00013161",
                  "WBGene00016250", "WBGene00003903", "WBGene00008852", "WBGene00001496", "WBGene00005663",
                  "WBGene00014472", "WBGene00003148", "WBGene00004744", "WBGene00011300", "WBGene00011564",
                  "WBGene00010470", "WBGene00019006", "WBGene00009799", "WBGene00001031", "WBGene00007980",
                  "WBGene00019983", "WBGene00003954", "WBGene00008405", "WBGene00045366", "WBGene00019457",
                  "WBGene00000502", "WBGene00018748", "WBGene00012004", "WBGene00007514", "WBGene00011303",
                  "WBGene00015176", "WBGene00019664", "WBGene00201458", "WBGene00021945", "WBGene00198236",
                  "WBGene00002019", "WBGene00002020", "WBGene00196396", "WBGene00002015")

stress.genes <- setdiff(stress.genes, mt.genes)

mt.genes <- intersect(mt.genes, rowData(wt.1)$id)
stress.genes <- intersect(stress.genes, rowData(wt.1)$id)

# From scater vignette
wt.1 <- addPerCellQC(wt.1, subsets = list(Mito = mt.genes, stress = stress.genes))
wt.1_corrected <- addPerCellQC(wt.1_corrected, subsets = list(Mito = mt.genes, stress = stress.genes))

colnames(colData(wt.1))
colnames(rowData(wt.1))

# Using calculateQCMetrics
summary(wt.1$subsets_Mito_percent)
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.000   1.003   1.932   5.984   7.369  70.283

summary(wt.1_corrected$subsets_Mito_percent)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   0.515   1.484   5.561   6.784  70.490

table(wt.1$subsets_Mito_percent < 20)
# FALSE  TRUE
#  538  5642

table(wt.1_corrected$subsets_Mito_percent < 20)
# FALSE  TRUE 
#  524  5656 

sizeFactors(wt.1) <- librarySizeFactors(wt.1)
wt.1 <- logNormCounts(wt.1)
wt.1
saveRDS(wt.1, "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/091622_wt.1_sce_uncorrected.rds")

sizeFactors(wt.1_corrected) <- librarySizeFactors(wt.1_corrected)
wt.1_corrected <- logNormCounts(wt.1_corrected)
wt.1_corrected
saveRDS(wt.1_corrected, "./3853-ST_new_alignment/3853-ST-1-selected/SoupX/091622_wt.1_sce_SoupX_corrected.rds")

# Sample 2

my_path <- "./3853-ST_new_alignment/3853-ST-2-selected/raw_feature_bc_matrix/"
wt.2 <- read10xCounts(my_path)

wt.2
head(colData(wt.2))
rownames(colData(wt.2)) <- colData(wt.2)$Barcode
colData(wt.2)$Sample <- "wt.2"
colData(wt.2)$Sex <- "Herm"
colData(wt.2)$Genotype <- "wt"
colData(wt.2)$stage <- "YA"
wt.2

head(rowData(wt.2))
rowData(wt.2) <- rowData(wt.2)[,1:2]
colnames(rowData(wt.2)) <- c("id", "gene_short_name")
rowData(wt.2)$gene_short_name <- i2s(rowData(wt.2)$id, gids)

bcrank <- barcodeRanks(counts(wt.2))

set.seed(100)
e.out <- emptyDrops(counts(wt.2), lower = 50)
summary(e.out$FDR <= 0.01)

# Mode     FALSE    TRUE    NA's 
# logical  4662    6215  660791

wt.2.ambient <- e.out@metadata$ambient
dim(wt.2.ambient)
wt.2.ambient[1:10,]
wt.2.ambient.df <- data.frame(row.names = rownames(wt.2.ambient),
                              id = rownames(wt.2.ambient),
                              gene_name = i2s(rownames(wt.2.ambient), gids),
                              expression = wt.2.ambient[,1])

wt.2.ambient.df <- wt.2.ambient.df %>% arrange(desc(expression))
head(wt.2.ambient.df, 15)

# rrn-3.1, mitochondrial genes (nduo-6, atp-6, MTCE.7, ctc-1, ctc-3), 
# nspc- genes (-20, -9, -14, -13, -19)

is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm = TRUE)

table(Limited=e.out$Limited, Significant=is.cell)
#          Significant
# Limited  FALSE  TRUE
#   FALSE  4662  639
# TRUE      0 5576

# Statistics on total UMI counts for cells called as cells

summary(e.out[which(e.out$FDR <= 0.01),]$Total)
# Min.  1st Qu.   Median    Mean    3rd Qu.    Max. 
#   51     109     669    1300    1655   61094 

wt.2$PValue <- e.out$PValue
head(colData(wt.2))

head(colData(wt.2))

is.cell <- e.out$FDR <= 0.01
wt.2 <- wt.2[, which(is.cell), drop = F]

# For SoupX
wt.2.filt.counts <- counts(wt.2)

# Writing to a new folder for SoupX.
# Generating new filtered matrix files for SoupX

writeMM(wt.2.filt.counts, "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/filtered_feature_bc_matrix/matrix.mtx")
system("gzip ./3853-ST_new_alignment/3853-ST-2-selected/SoupX/filtered_feature_bc_matrix/matrix.mtx")
writeMM(wt.2.filt.counts, "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/filtered_feature_bc_matrix/matrix.mtx")
fil.barcodes <- colnames(wt.2)
write.table(fil.barcodes, "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system("gzip ./3853-ST_new_alignment/3853-ST-2-selected/SoupX/filtered_feature_bc_matrix/barcodes.tsv")
write.table(fil.barcodes, "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

sizeFactors(wt.2) <- librarySizeFactors(wt.2)
wt.2 <- logNormCounts(wt.2)

wt.2.S <- as.Seurat(wt.2)

# Generating initial clustering with Seurat to aid in SoupX background correction
wt.2.S <- NormalizeData(wt.2.S)
wt.2.S <- FindVariableFeatures(wt.2.S, selection.method = "vst", nfeatures = 4000)
wt.2.S <- ScaleData(wt.2.S)
wt.2.S <- RunPCA(wt.2.S, features = VariableFeatures(wt.2.S), npcs = 75)

ElbowPlot(wt.2.S, ndims = 75)
wt.2.S <- FindNeighbors(wt.2.S, dims = 1:75)
wt.2.S <- FindClusters(wt.2.S, resolution = 1.2)

wt.2.S <- RunUMAP(wt.2.S, dims = 1:75)
DimPlot(wt.2.S, reduction = "umap")

# Seurat identified 33 clusters. This is from the 2nd wt sample.

FeaturePlot(wt.2.S, features = "WBGene00011392")
FeaturePlot(wt.2.S, features = "WBGene00000668")
FeaturePlot(wt.2.S, features = "WBGene00003934")
FeaturePlot(wt.2.S, features = s2i("col-140", gids))
FeaturePlot(wt.2.S, features = "WBGene00006762")
FeaturePlot(wt.2.S, features = "WBGene00006756")
FeaturePlot(wt.1.S, features = "WBGene00003024")
FeaturePlot(wt.2.S, features = "WBGene00003024")
FeaturePlot(wt.2.S, features = "WBGene00006744")

# Checking expression of npr-2
FeaturePlot(wt.2.S, features = "WBGene00003808")

# Checking expression of unc-129
FeaturePlot(wt.2.S, features = "WBGene00006852")

# Checking expression of bnc-1
FeaturePlot(wt.2.S, features = "WBGene00044791")

# Checking expression of ceh-12
FeaturePlot(wt.2.S, features = "WBGene00000436")

# vab-7
FeaturePlot(wt.2.S, features = "WBGene00006873")

class(wt.2.S@reductions$umap)
umap <- wt.2.S@reductions$umap

Embeddings(umap)[1:3,1:2]
wt.2_DR <- as.data.frame(Embeddings(umap))
wt.2_DR$Cluster <- Idents(wt.2.S)
head(wt.2_DR)

# Now running SoupX, with soupRange = c(0, 25)

my_dir_S <- "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/"
scl.wt.2 <- load10X(my_dir_S, soupRange = c(0, 25))

scl.wt.2$toc[1:5,1:5]
head(wt.2_DR)

head(scl.wt.2$soupProfile[order(scl.wt.2$soupProfile$est, decreasing = TRUE), ], n = 30)

# Now the genes are listed as WBGeneIDs.

scl.wt.2$soupProfile$gene_short_name <- i2s(rownames(scl.wt.2$soupProfile), gids)

head(scl.wt.2$soupProfile[order(scl.wt.2$soupProfile$est, decreasing = TRUE), ], n = 30)

# top genes are ribosomal RNA (rrn-3.1), mitochondrial genes (nduo-6, atp-6, ctc genes), 
# flp-1, Y73F4A.1, pat-10, cup-4, nlp-17

plotMarkerMap(scl.wt.2, s2i("pat-10", gids), wt.2_DR)
plotMarkerMap(scl.wt.2, s2i("R102.2", gids), wt.2_DR)
scl.wt.2 <- setDR(scl.wt.2, wt.2_DR)
plotMarkerMap(scl.wt.2, s2i("flp-28", gids))
plotMarkerMap(scl.wt.2, s2i("nspc-20", gids))
plotMarkerMap(scl.wt.2, s2i("flp-12", gids))

est.cont <- s2i(c("pat-10", "cpn-3", "flp-28", "nspc-20", "nspc-19", "nspc-10", "nspc-13", "nspc-14", "nspc-9", "nspc-18", "nspc-15"), gids)

useToEst = estimateNonExpressingCells(scl.wt.2, nonExpressedGeneList = list(wt.2 = est.cont), 
                                      clusters = setNames(wt.2_DR$Cluster, rownames(wt.2_DR)))

plotMarkerMap(scl.wt.2, geneSet = est.cont, DR = wt.2_DR, useToEst = useToEst)
# Cells used to estimate the contamination

scl.wt.2 <- setClusters(scl.wt.2, wt.2_DR$Cluster)

scl.wt.2 <- calculateContaminationFraction(scl.wt.2, list(wt.2 = est.cont), useToEst = useToEst)
# Estimated global contamination fraction of 1.85%

system.time(out <- adjustCounts(scl.wt.2, roundToInt = TRUE))

cntSoggy = rowSums(scl.wt.2$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed
i2s(names(mostZeroed),gids)
tail(sort(rowSums(scl.wt.2$toc > out)/rowSums(scl.wt.2$toc > 0)), n = 20)

plotChangeMap(scl.wt.2, out, s2i("col-124", gids))
plotChangeMap(scl.wt.2, out, s2i("nlp-17", gids))
plotChangeMap(scl.wt.2, out, s2i("pat-10", gids))
plotChangeMap(scl.wt.2, out, s2i("flp-28", gids))
plotChangeMap(scl.wt.2, out, s2i("nspc-20", gids))
plotChangeMap(scl.wt.2, out, s2i("flp-12", gids))

out[1:5,1:5]
scl.wt.2$toc[1:5,1:5]

saveRDS(scl.wt.2, "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/091622_soupX_object_wt.2.rds")
saveRDS(out, "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/091622_soupX_corrected_matrix_wt.2.rds")

dim(out)
head(colnames(out))
head(rownames(out))

head(colnames(wt.2.filt.counts))
head(rownames(wt.2.filt.counts))

wt.2_corrected <- SingleCellExperiment(assays = list(counts = out), colData = colData(wt.2))  

rowData(wt.2_corrected) <- rowData(wt.2)

# Using scater to calculate some quality control metrics

# From scater vignette
wt.2 <- addPerCellQC(wt.2, subsets = list(Mito = mt.genes, stress = stress.genes))
wt.2_corrected <- addPerCellQC(wt.2_corrected, subsets = list(Mito = mt.genes, stress = stress.genes))

colnames(colData(wt.2))
colnames(rowData(wt.2))

# Using calculateQCMetrics
summary(wt.2$subsets_Mito_percent)
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.000   1.089   2.261   6.973   9.813  92.596 

summary(wt.2_corrected$subsets_Mito_percent)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.6081  1.7857  6.6119  9.4340 93.5455 

table(wt.2$subsets_Mito_percent < 20)
# FALSE  TRUE
#   733  5482  

table(wt.2_corrected$subsets_Mito_percent < 20)
# FALSE  TRUE 
#   718  5497

sizeFactors(wt.2) <- librarySizeFactors(wt.2)
wt.2 <- logNormCounts(wt.2)
wt.2
saveRDS(wt.2, "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/091622_wt.2_sce_uncorrected.rds")

sizeFactors(wt.2_corrected) <- librarySizeFactors(wt.2_corrected)
wt.2_corrected <- logNormCounts(wt.2_corrected)
wt.2_corrected
saveRDS(wt.2_corrected, "./3853-ST_new_alignment/3853-ST-2-selected/SoupX/091622_wt.2_sce_SoupX_corrected.rds")

# Sample 3 
# 6973-ST

my_path <- "./6973-ST/6973-ST-1/raw_feature_bc_matrix/"
wt.3 <- read10xCounts(my_path)

wt.3
head(colData(wt.3))
rownames(colData(wt.3)) <- colData(wt.3)$Barcode
colData(wt.3)$Sample <- "wt"
colData(wt.3)$Genotype <- "wt"
wt.3

head(rowData(wt.3))
rowData(wt.3) <- rowData(wt.3)[,1:2]
colnames(rowData(wt.3)) <- c("id", "gene_short_name")
rowData(wt.3)$gene_short_name <- i2s(rowData(wt.3)$id, gids)

bcrank <- barcodeRanks(counts(wt.3))

# Using droplets with fewer than 50 UMIs as the background. 

set.seed(100)
e.out <- emptyDrops(counts(wt.3), lower = 50)
summary(e.out$FDR <= 0.01)

# Mode     FALSE    TRUE    NA's 
# logical  17984   12214  516424

wt.3.ambient <- e.out@metadata$ambient
dim(wt.3.ambient)
wt.3.ambient[1:10,]
wt.3.ambient.df <- data.frame(row.names = rownames(wt.3.ambient),
                            id = rownames(wt.3.ambient),
                            gene_name = i2s(rownames(wt.3.ambient), gids),
                            expression = wt.3.ambient[,1])

wt.3.ambient.df <- wt.3.ambient.df %>% arrange(desc(expression))
head(wt.3.ambient.df, 15)

# rrn-3.1, mitochondrial genes (nduo-6, atp-6, MTCE.7, ctc-1, ctc-3), pat-10 (muscle), 
# some neuronal genes (flp-12, unc-25)

is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm = TRUE)

table(Limited=e.out$Limited, Significant=is.cell)
#          Significant
# Limited  FALSE  TRUE
#   FALSE  17984  1807
#    TRUE      0 10407

# Statistics on total UMI counts for cells called as cells

summary(e.out[which(e.out$FDR <= 0.01),]$Total)
# Min.  1st Qu.   Median    Mean    3rd Qu.    Max. 
#   51     156       861    1319    1755      34173 

wt.3$PValue <- e.out$PValue
head(colData(wt.3))

wt.3$FDR.lower75 <- e.out$FDR
head(colData(wt.3))

is.cell <- e.out$FDR <= 0.01
wt.3.3 <- wt.3[, which(is.cell), drop = F]

# For SoupX
wt.3.filt.counts <- counts(wt.3)

# Writing to a new folder for SoupX.
# Generating new filtered matrix files for SoupX
library(Matrix)
writeMM(wt.3.filt.counts, "./6973-ST/SoupX/6973-ST-1/filtered_feature_bc_matrix/matrix.mtx")
system("gzip ./6973-ST/SoupX/6973-ST-1/filtered_feature_bc_matrix/matrix.mtx")
writeMM(wt.3.filt.counts, "./6973-ST/SoupX/6973-ST-1/filtered_feature_bc_matrix/matrix.mtx")
fil.barcodes <- colnames(wt.3)
write.table(fil.barcodes, "./6973-ST/SoupX/6973-ST-1/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system("gzip ./6973-ST/SoupX/6973-ST-1/filtered_feature_bc_matrix/barcodes.tsv")
write.table(fil.barcodes, "./6973-ST/SoupX/6973-ST-1/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

sizeFactors(wt.3) <- librarySizeFactors(wt.3)
wt.3 <- logNormCounts(wt.3)

library(Seurat)
wt.3.S <- as.Seurat(wt.3)

# Generating initial clustering with Seurat to aid in SoupX background correction
wt.3.S <- NormalizeData(wt.3.S)
wt.3.S <- FindVariableFeatures(wt.3.S, selection.method = "vst", nfeatures = 4000)
wt.3.S <- ScaleData(wt.3.S)
wt.3.S <- RunPCA(wt.3.S, features = VariableFeatures(wt.3.S), npcs = 100)

ElbowPlot(wt.3.S, ndims = 100)
wt.3.S <- FindNeighbors(wt.3.S, dims = 1:75)
wt.3.S <- FindClusters(wt.3.S, resolution = 1.2)

wt.3.S <- RunUMAP(wt.3.S, dims = 1:75)
DimPlot(wt.3.S, reduction = "umap")

# Seurat identified 40 clusters. There is good separation between many clusters. 
# This is the sample from the wt.3 lin-39::tagRFP.

# Checking expression of sbt-1 (neuronal marker)
FeaturePlot(wt.3.S, features = "WBGene00011392")

# Many clusters are sbt-1 positive, but several clusters are not. B

# Checking expression of col-93 (epidermal marker)
FeaturePlot(wt.3.S, features = "WBGene00000668")
# Clear epidermal clusters

# pat-10 (muscle marker)
FeaturePlot(wt.3.S, features = "WBGene00003934")

# several distinct clusters positive for pat-10

FeaturePlot(wt.3.S, features = "WBGene00000779")
FeaturePlot(wt.3.S, features = "WBGene00006789")
# This is unc-54, it's widespread, more intense in the neuronal clusters, as expected from
# the use of the 3' UTR in the reporter.

# Checking expression of unc-25
FeaturePlot(wt.3.S, features = "WBGene00006762")

# Checking expression of unc-17
FeaturePlot(wt.3.S, features = "WBGene00006756")

# Checking expression of unc-4
FeaturePlot(wt.3.S, features = "WBGene00006744")

# Checking expression of npr-2
FeaturePlot(wt.3.S, features = "WBGene00003808")

# Checking expression of unc-129
FeaturePlot(wt.3.S, features = "WBGene00006852")

# Checking expression of bnc-1
FeaturePlot(wt.3.S, features = "WBGene00044791")

# Checking expression of ceh-12
FeaturePlot(wt.3.S, features = "WBGene00000436")

# vab-7
FeaturePlot(wt.3.S, features = "WBGene00006873")

# Widespread low expression, no clear clustering

class(wt.3.S@reductions$umap)
umap <- wt.3.S@reductions$umap

Embeddings(umap)[1:3,1:2]
wt.3_DR <- as.data.frame(Embeddings(umap))
wt.3_DR$Cluster <- Idents(wt.3.S)
head(wt.3_DR)

# Now I will try running SoupX, with soupRange = c(0, 25)
library(SoupX)
my_dir_S <- "./6973-ST/SoupX/6973-ST-1"
scl.wt.3 <- load10X(my_dir_S, soupRange = c(0, 25))

scl.wt.3$toc[1:5,1:5]
head(wt.3_DR)

head(scl.wt.3$soupProfile[order(scl.wt.3$soupProfile$est, decreasing = TRUE), ], n = 30)

# Now the genes are listed as WBGeneIDs.

scl.wt.3$soupProfile$gene_short_name <- i2s(rownames(scl.wt.3$soupProfile), gids)

head(scl.wt.3$soupProfile[order(scl.wt.3$soupProfile$est, decreasing = TRUE), ], n = 30)

# top genes are ribosomal RNA (rrn-3.1), mitochondrial genes (nduo-6, atp-6, ctc genes), 
# muscle genes (cpn-3, mlc-3, pat-10, lev-11, unc-54) and collagens (col-140, col-124)

plotMarkerMap(scl.wt.3, s2i("pat-10", gids), wt.3_DR)

plotMarkerMap(scl.wt.3, s2i("mlc-3", gids), wt.3_DR)

scl.wt.3.3 <- setDR(scl.wt.3, wt.3_DR)

plotMarkerMap(scl.wt.3, s2i("col-140", gids))
plotMarkerMap(scl.wt.3, s2i("nduo-6", gids))
plotMarkerMap(scl.wt.3, s2i("cpn-3", gids))

plotMarkerMap(scl.wt.3, s2i("flp-12", gids))
plotMarkerMap(scl.wt.3, s2i("clik-1", gids))
plotMarkerMap(scl.wt.3, s2i("unc-25", gids))

est.cont <- s2i(c("pat-10", "lev-11", "mlc-3", "col-140", "col-20", "col-124", 
                  "cpn-3", "clik-1", "mlc-2", "flp-12"), gids)

useToEst = estimateNonExpressingCells(scl.wt.3, nonExpressedGeneList = list(wt.3.3 = est.cont), 
                                      clusters = setNames(wt.3_DR$Cluster, rownames(wt.3_DR)))

plotMarkerMap(scl.wt.3, geneSet = est.cont, DR = wt.3_DR, useToEst = useToEst)
# Cells used to estimate the contamination

scl.wt.3 <- setClusters(scl.wt.3, wt.3_DR$Cluster)

scl.wt.3 <- calculateContaminationFraction(scl.wt.3, list(wt.3 = est.cont), useToEst = useToEst)
# Estimated global contamination fraction of 3.63%


system.time(out <- adjustCounts(scl.wt.3, roundToInt = TRUE))

cntSoggy = rowSums(scl.wt.3$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

tail(sort(rowSums(scl.wt.3$toc > out)/rowSums(scl.wt.3$toc > 0)), n = 20)

plotChangeMap(scl.wt.3, out, s2i("col-124", gids))
plotChangeMap(scl.wt.3, out, s2i("col-140", gids))
plotChangeMap(scl.wt.3, out, s2i("pat-10", gids))
plotChangeMap(scl.wt.3, out, s2i("mlc-3", gids))
plotChangeMap(scl.wt.3, out, s2i("cpn-3", gids))
plotChangeMap(scl.wt.3, out, s2i("flp-12", gids))

out[1:5,1:5]
scl.wt.3$toc[1:5,1:5]

saveRDS(scl.wt.3, "./6973-ST/SoupX/6973-ST-1/111621_soupX_object_wt.3.rds")
saveRDS(out, "./6973-ST/SoupX/6973-ST-1/111621_soupX_corrected_matrix_wt.3.rds")

dim(out)
head(colnames(out))
head(rownames(out))

head(colnames(wt.3.filt.counts))
head(rownames(wt.3.filt.counts))

# I need to make the rownames consistent with the uncorrected counts (WBGeneIDs)

wt.3_corrected <- SingleCellExperiment(assays = list(counts = out), colData = colData(wt.3))  

rowData(wt.3_corrected) <- rowData(wt.3)

# Using scater to calculate some quality control metrics

# Because MTCE.7 and MTCE.33 are both high in the soup, I'm going to add them in.

mt.genes <- c("WBGene00010959", "WBGene00010961", "WBGene00010966", "WBGene00010963", "WBGene00010957", "WBGene00010967", "WBGene00010958", "WBGene00010960", "WBGene00010962","WBGene00000829","WBGene00010965","WBGene00010964","WBGene00014454","WBGene00014472")
stress.genes <- c("WBGene00009692", "WBGene00009691", "WBGene00002018", "WBGene00002016", "WBGene00002017",
                  "WBGene00006959", "WBGene00002026", "WBGene00004622", "WBGene00235164", "WBGene00004131",
                  "WBGene00012813", "WBGene00000097", "WBGene00015168", "WBGene00009180", "WBGene00013161",
                  "WBGene00016250", "WBGene00003903", "WBGene00008852", "WBGene00001496", "WBGene00005663",
                  "WBGene00014472", "WBGene00003148", "WBGene00004744", "WBGene00011300", "WBGene00011564",
                  "WBGene00010470", "WBGene00019006", "WBGene00009799", "WBGene00001031", "WBGene00007980",
                  "WBGene00019983", "WBGene00003954", "WBGene00008405", "WBGene00045366", "WBGene00019457",
                  "WBGene00000502", "WBGene00018748", "WBGene00012004", "WBGene00007514", "WBGene00011303",
                  "WBGene00015176", "WBGene00019664", "WBGene00201458", "WBGene00021945", "WBGene00198236",
                  "WBGene00002019", "WBGene00002020", "WBGene00196396", "WBGene00002015")

stress.genes <- setdiff(stress.genes, mt.genes)

cellQC <- perCellQCMetrics(wt.3, subsets = list(Mito = mt.genes, stress = stress.genes))

str(cellQC)

# From scater vignette
wt.3 <- addPerCellQC(wt.3, subsets = list(Mito = mt.genes, stress = stress.genes))
wt.3_corrected <- addPerCellQC(wt.3_corrected, subsets = list(Mito = mt.genes, stress = stress.genes))

colnames(colData(wt.3))
colnames(rowData(wt.3))

# Using calculateQCMetrics
summary(wt.3$subsets_Mito_percent)
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    0.0000  0.7956  1.7958  4.1995  5.0940 48.4849 

summary(wt.3_corrected$subsets_Mito_percent)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.7565  3.3803  4.0573 48.8189

table(wt.3$subsets_Mito_percent < 20)
# FALSE  TRUE
#   372 11842  

table(wt.3_corrected$subsets_Mito_percent < 20)
# FALSE  TRUE 
#   324 11890 

# Plotting UMI (x) by genes (y)
plotColData(wt.3, x = "sum", y = "detected")

plotColData(wt.3_corrected, x = "sum", y = "detected")

# Also using logNormCounts instead of normalize()
sizeFactors(wt.3) <- librarySizeFactors(wt.3)
wt.3 <- logNormCounts(wt.3)
wt.3
saveRDS(wt.3, "./6973-ST/SoupX/6973-ST-1/111621_wt.3_sce_uncorrected.rds")

sizeFactors(wt.3_corrected) <- librarySizeFactors(wt.3_corrected)
wt.3_corrected <- logNormCounts(wt.3_corrected)
wt.3_corrected
saveRDS(wt.3_corrected, "./6973-ST/SoupX/6973-ST-1/111621_wt.3_sce_SoupX_corrected.rds")

# Sample 4 (lin-39::RFP rep 2), 7309-ST-1

my_path <- "./7309-ST/7309-ST-1/raw_feature_bc_matrix/"
wt.4 <- read10xCounts(my_path)

wt.4
head(colData(wt.4))
rownames(colData(wt.4)) <- colData(wt.4)$Barcode
colData(wt.4)$Sample <- "wt"
colData(wt.4)$Genotype <- "wt"
wt.4

head(rowData(wt.4))
rowData(wt.4) <- rowData(wt.4)[,1:2]
colnames(rowData(wt.4)) <- c("id", "gene_short_name")
rowData(wt.4)$gene_short_name <- i2s(rowData(wt.4)$id, gids)

bcrank <- barcodeRanks(counts(wt.4))

# Using droplets with fewer than 50 UMIs as the background. 

set.seed(100)
e.out <- emptyDrops(counts(wt.4), lower = 50)
summary(e.out$FDR <= 0.01)

# Mode     FALSE    TRUE    NA's 
# logical  8056    8183  513939 

wt.4.ambient <- e.out@metadata$ambient
dim(wt.4.ambient)
wt.4.ambient[1:10,]
wt.4.ambient.df <- data.frame(row.names = rownames(wt.4.ambient),
                            id = rownames(wt.4.ambient),
                            gene_name = i2s(rownames(wt.4.ambient), gids),
                            expression = wt.4.ambient[,1])

wt.4.ambient.df <- wt.4.ambient.df %>% arrange(desc(expression))
head(wt.4.ambient.df, 15)

# rrn-3.1, mitochondrial genes (nduo-6, atp-6, MTCE.7, ctc-1, ctc-3), pat-10 (muscle), 
# some neuronal genes (flp-12, flp-11, unc-25)

# Other diagnostics

is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm = TRUE)

table(Limited=e.out$Limited, Significant=is.cell)
#          Significant
# Limited  FALSE  TRUE
#   FALSE  8056  762
#   TRUE      0 7421

# Because Limited == TRUE and Significant == FALSE, it means the number of permutations 
# was not limiting.

# Statistics on total UMI counts for cells called as cells

summary(e.out[which(e.out$FDR <= 0.01),]$Total)
# Min.  1st Qu.   Median    Mean    3rd Qu.    Max. 
#   51     296    1177    1563    2069   39289  

wt.4$PValue <- e.out$PValue
head(colData(wt.4))

wt.4$FDR.lower75 <- e.out$FDR
head(colData(wt.4))

is.cell <- e.out$FDR <= 0.01
wt.4.4 <- wt.4[, which(is.cell), drop = F]

# For SoupX
wt.4.filt.counts <- counts(wt.4)

# Writing to a new folder for SoupX.
# Generating new filtered matrix files for SoupX
library(Matrix)
writeMM(wt.4.filt.counts, "./7309-ST/SoupX/7309-ST-1/filtered_feature_bc_matrix/matrix.mtx")
system("gzip ./7309-ST/SoupX/7309-ST-1/filtered_feature_bc_matrix/matrix.mtx")
writeMM(wt.4.filt.counts, "./7309-ST/SoupX/7309-ST-1/filtered_feature_bc_matrix/matrix.mtx")
fil.barcodes <- colnames(wt.4)
write.table(fil.barcodes, "./7309-ST/SoupX/7309-ST-1/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system("gzip ./7309-ST/SoupX/7309-ST-1/filtered_feature_bc_matrix/barcodes.tsv")
write.table(fil.barcodes, "./7309-ST/SoupX/7309-ST-1/filtered_feature_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

sizeFactors(wt.4) <- librarySizeFactors(wt.4)
wt.4 <- logNormCounts(wt.4)

library(Seurat)
wt.4.S <- as.Seurat(wt.4)

# Generating initial clustering with Seurat to aid in SoupX background correction
wt.4.S <- NormalizeData(wt.4.S)
wt.4.S <- FindVariableFeatures(wt.4.S, selection.method = "vst", nfeatures = 4000)
wt.4.S <- ScaleData(wt.4.S)
wt.4.S <- RunPCA(wt.4.S, features = VariableFeatures(wt.4.S), npcs = 100)

ElbowPlot(wt.4.S, ndims = 100)
wt.4.S <- FindNeighbors(wt.4.S, dims = 1:75)
wt.4.S <- FindClusters(wt.4.S, resolution = 1.2)

wt.4.S <- RunUMAP(wt.4.S, dims = 1:75)
DimPlot(wt.4.S, reduction = "umap")

# Seurat identified 37 clusters. There is good separation between many clusters. 
# This is the sample from the wt.4 lin-39::tagRFP.

# Checking expression of sbt-1 (neuronal marker)
FeaturePlot(wt.4.S, features = "WBGene00011392")

# Many clusters are sbt-1 positive, but several clusters are not. But I did have my
# fluorescent gates set pretty low.

# Checking expression of col-93 (epidermal marker)
FeaturePlot(wt.4.S, features = "WBGene00000668")
# Clear epidermal clusters

# pat-10 (muscle marker)
FeaturePlot(wt.4.S, features = "WBGene00003934")

# several distinct clusters positive for pat-10

FeaturePlot(wt.4.S, features = "WBGene00000779")
FeaturePlot(wt.4.S, features = "WBGene00006789")
# This is unc-54, it's widespread, more intense in the neuronal clusters, as expected from
# the use of the 3' UTR in the reporter.

# Checking expression of unc-25
FeaturePlot(wt.4.S, features = "WBGene00006762")

# Checking expression of unc-17
FeaturePlot(wt.4.S, features = "WBGene00006756")

# Checking expression of unc-4
FeaturePlot(wt.4.S, features = "WBGene00006744")

# Checking expression of npr-2
FeaturePlot(wt.4.S, features = "WBGene00003808")

# Checking expression of unc-129
FeaturePlot(wt.4.S, features = "WBGene00006852")

# Checking expression of bnc-1
FeaturePlot(wt.4.S, features = "WBGene00044791")

# Checking expression of ceh-12
FeaturePlot(wt.4.S, features = "WBGene00000436")

# vab-7
FeaturePlot(wt.4.S, features = "WBGene00006873")

# Widespread low expression, no clear clustering

class(wt.4.S@reductions$umap)
umap <- wt.4.S@reductions$umap

Embeddings(umap)[1:3,1:2]
wt.4_DR <- as.data.frame(Embeddings(umap))
wt.4_DR$Cluster <- Idents(wt.4.S)
head(wt.4_DR)

# Now I will try running SoupX, with soupRange = c(0, 25)
library(SoupX)
my_dir_S <- "./7309-ST/SoupX/7309-ST-1"
scl.wt.4 <- load10X(my_dir_S, soupRange = c(0, 25))

scl.wt.4$toc[1:5,1:5]
head(wt.4_DR)

head(scl.wt.4$soupProfile[order(scl.wt.4$soupProfile$est, decreasing = TRUE), ], n = 30)

# Now the genes are listed as WBGeneIDs.

scl.wt.4$soupProfile$gene_short_name <- i2s(rownames(scl.wt.4$soupProfile), gids)

head(scl.wt.4$soupProfile[order(scl.wt.4$soupProfile$est, decreasing = TRUE), ], n = 30)

# top genes are ribosomal RNA (rrn-3.1), mitochondrial genes (nduo-6, atp-6, ctc genes), 
# muscle genes (cpn-3, mlc-3, pat-10, lev-11, unc-54) and collagens (col-140, col-124)

plotMarkerMap(scl.wt.4, s2i("pat-10", gids), wt.4_DR)

plotMarkerMap(scl.wt.4, s2i("mlc-3", gids), wt.4_DR)

scl.wt.4 <- setDR(scl.wt.4, wt.4_DR)

plotMarkerMap(scl.wt.4, s2i("col-140", gids))
plotMarkerMap(scl.wt.4, s2i("nduo-6", gids))
plotMarkerMap(scl.wt.4, s2i("cpn-3", gids))

plotMarkerMap(scl.wt.4, s2i("flp-12", gids))
plotMarkerMap(scl.wt.4, s2i("clik-1", gids))
plotMarkerMap(scl.wt.4, s2i("unc-25", gids))

est.cont <- s2i(c("pat-10", "lev-11", "mlc-3", "col-140", "col-20", "col-124", 
                  "cpn-3", "clik-1", "mlc-2", "flp-12"), gids)

useToEst = estimateNonExpressingCells(scl.wt.4, nonExpressedGeneList = list(wt.4 = est.cont), 
                                      clusters = setNames(wt.4_DR$Cluster, rownames(wt.4_DR)))

plotMarkerMap(scl.wt.4, geneSet = est.cont, DR = wt.4_DR, useToEst = useToEst)
# Cells used to estimate the contamination

scl.wt.4 <- setClusters(scl.wt.4, wt.4_DR$Cluster)

scl.wt.4 <- calculateContaminationFraction(scl.wt.4, list(wt.4 = est.cont), useToEst = useToEst)
# Estimated global contamination fraction of 2.92%

system.time(out <- adjustCounts(scl.wt.4, roundToInt = TRUE))

cntSoggy = rowSums(scl.wt.4$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

tail(sort(rowSums(scl.wt.4$toc > out)/rowSums(scl.wt.4$toc > 0)), n = 20)

plotChangeMap(scl.wt.4, out, s2i("col-124", gids))
plotChangeMap(scl.wt.4, out, s2i("col-140", gids))
plotChangeMap(scl.wt.4, out, s2i("pat-10", gids))
plotChangeMap(scl.wt.4, out, s2i("mlc-3", gids))
plotChangeMap(scl.wt.4, out, s2i("cpn-3", gids))
plotChangeMap(scl.wt.4, out, s2i("flp-12", gids))

out[1:5,1:5]
scl.wt.4$toc[1:5,1:5]

saveRDS(scl.wt.4, "./7309-ST/SoupX/7309-ST-1/112921_soupX_object_wt.4.rds")
saveRDS(out, "./7309-ST/SoupX/7309-ST-1/112921_soupX_corrected_matrix_wt.4.rds")

dim(out)
head(colnames(out))
head(rownames(out))

head(colnames(wt.4.filt.counts))
head(rownames(wt.4.filt.counts))

# I need to make the rownames consistent with the uncorrected counts (WBGeneIDs)

wt.4_corrected <- SingleCellExperiment(assays = list(counts = out), colData = colData(wt.4))  

rowData(wt.4_corrected) <- rowData(wt.4)

# Using scater to calculate some quality control metrics

# Because MTCE.7 and MTCE.33 are both high in the soup, I'm going to add them in.

mt.genes <- c("WBGene00010959", "WBGene00010961", "WBGene00010966", "WBGene00010963", "WBGene00010957", "WBGene00010967", "WBGene00010958", "WBGene00010960", "WBGene00010962","WBGene00000829","WBGene00010965","WBGene00010964","WBGene00014454","WBGene00014472")
stress.genes <- c("WBGene00009692", "WBGene00009691", "WBGene00002018", "WBGene00002016", "WBGene00002017",
                  "WBGene00006959", "WBGene00002026", "WBGene00004622", "WBGene00235164", "WBGene00004131",
                  "WBGene00012813", "WBGene00000097", "WBGene00015168", "WBGene00009180", "WBGene00013161",
                  "WBGene00016250", "WBGene00003903", "WBGene00008852", "WBGene00001496", "WBGene00005663",
                  "WBGene00014472", "WBGene00003148", "WBGene00004744", "WBGene00011300", "WBGene00011564",
                  "WBGene00010470", "WBGene00019006", "WBGene00009799", "WBGene00001031", "WBGene00007980",
                  "WBGene00019983", "WBGene00003954", "WBGene00008405", "WBGene00045366", "WBGene00019457",
                  "WBGene00000502", "WBGene00018748", "WBGene00012004", "WBGene00007514", "WBGene00011303",
                  "WBGene00015176", "WBGene00019664", "WBGene00201458", "WBGene00021945", "WBGene00198236",
                  "WBGene00002019", "WBGene00002020", "WBGene00196396", "WBGene00002015")

stress.genes <- setdiff(stress.genes, mt.genes)

cellQC <- perCellQCMetrics(wt.4, subsets = list(Mito = mt.genes, stress = stress.genes))

str(cellQC)

# From scater vignette
wt.4 <- addPerCellQC(wt.4, subsets = list(Mito = mt.genes, stress = stress.genes))
wt.4_corrected <- addPerCellQC(wt.4_corrected, subsets = list(Mito = mt.genes, stress = stress.genes))

colnames(colData(wt.4))
colnames(rowData(wt.4))

# Using calculateQCMetrics
summary(wt.4$subsets_Mito_percent)
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    0.0000  0.4236  0.9630  2.3988  2.2806 47.7612 

summary(wt.4_corrected$subsets_Mito_percent)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0402  0.4027  1.9408  1.7019 48.0122 

table(wt.4$subsets_Mito_percent < 20)
# FALSE  TRUE
#   106  8077  

table(wt.4_corrected$subsets_Mito_percent < 20)
# FALSE  TRUE 
#    95  8088 

# Plotting UMI (x) by genes (y)
plotColData(wt.4, x = "sum", y = "detected")

plotColData(wt.4_corrected, x = "sum", y = "detected")

# Also using logNormCounts instead of normalize()
sizeFactors(wt.4) <- librarySizeFactors(wt.4)
wt.4 <- logNormCounts(wt.4)
wt.4
saveRDS(wt.4, "./7309-ST/SoupX/7309-ST-1/112921_wt.4_sce_uncorrected.rds")

sizeFactors(wt.4_corrected) <- librarySizeFactors(wt.4_corrected)
wt.4_corrected <- logNormCounts(wt.4_corrected)
wt.4_corrected
saveRDS(wt.4_corrected, "./7309-ST/SoupX/7309-ST-1/112921_wt.4_sce_SoupX_corrected.rds")

# Merging the datasets together
# Coding an experiment ID
colData(wt.1_corrected)$Experiment <- "acr-2_YA_wt_1"
colData(wt.2_corrected)$Experiment <- "acr-2_YA_wt_2"
colData(wt.3_corrected)$Experiment <- "lin-39_YA_wt_1"
colData(wt.4_corrected)$Experiment <- "lin-39_YA_wt_2"

colData(wt.1_corrected)$strain <- "acr-2::GFP"
colData(wt.2_corrected)$strain <- "acr-2::GFP"
colData(wt.3_corrected)$strain <- "lin-39::RFP"
colData(wt.4_corrected)$strain <- "lin-39::RFP"

# generating unique barcodes including experiment ID
colData(wt.1_corrected)$unique_barcode <- paste(colData(wt.1_corrected)$Barcode, 
                                                colData(wt.1_corrected)$Experiment, sep = "_")
colData(wt.2_corrected)$unique_barcode <- paste(colData(wt.2_corrected)$Barcode, 
                                                colData(wt.2_corrected)$Experiment, sep = "_")
colData(wt.3_corrected)$unique_barcode <- paste(colData(wt.3_corrected)$Barcode, 
                                                colData(wt.3_corrected)$Experiment, sep = "_")
colData(wt.4_corrected)$unique_barcode <- paste(colData(wt.4_corrected)$Barcode, 
                                                colData(wt.4_corrected)$Experiment, sep = "_")

colnames(wt.1_corrected) <- colData(wt.1_corrected)$unique_barcode
colnames(wt.2_corrected) <- colData(wt.2_corrected)$unique_barcode
colnames(wt.3_corrected) <- colData(wt.3_corrected)$unique_barcode
colnames(wt.4_corrected) <- colData(wt.4_corrected)$unique_barcode

acr2.lin39 <- combine_cds(list(wt.1_corrected, wt.2_corrected,
                               wt.3_corrected, wt.4_corrected),
                          cell_names_unique = T,
                          keep_all_genes = T)

# Keeping cells with fewer than 20% of UMIs from mitochondrial genes

acr2.lin39 <- acr2.lin39[,colData(acr2.lin39)$subsets_Mito_percent < 20]
acr2.lin39 <- detect_genes(acr2.lin39)

# Keeping only genes found in more than 5 single cells
acr2.lin39 <- acr2.lin39[rowData(acr2.lin39)$num_cells_expressed > 5,]
acr2.lin39

# Performing dimensionality reduction

acr2.lin39 <- preprocess_cds(acr2.lin39, num_dim = 50)
plot_pc_variance_explained(acr2.lin39)

acr2.lin39 <- align_cds(acr2.lin39, alignment_group = "strain", alignment_k = 5)
acr2.lin39 <- reduce_dimension(acr2.lin39, 
                               reduction_method = "UMAP",
                               preprocess_method = "Aligned",
                               umap.min_dist = 0.3,
                               umap.n_neighbors = 75)

acr2.lin39 <- cluster_cells(acr2.lin39,
                            res = 3e-4)

# Storing the UMAP coordinates in the cell metadata table for custom functions
colData(acr2.lin39)$UMAP_1 <- reducedDims(acr2.lin39)[["UMAP"]][,1]
colData(acr2.lin39)$UMAP_2 <- reducedDims(acr2.lin39)[["UMAP"]][,2]

plot_cells(acr2.lin39, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 0.5)

plot_cells(acr2.lin39, 
           color_cells_by = "partition",
           label_groups_by_cluster = F,
           group_label_size = 3,
           cell_size = 0.5)

# storing the cluster and partition identities in the cell metadata for subsetting and other functions
colData(acr2.lin39)$cluster <- monocle3::clusters(acr2.lin39)
colData(acr2.lin39)$partition <- monocle3::partitions(acr2.lin39)

plot_cells(acr2.lin39,
           color_cells_by = "Experiment",
           label_cell_groups = F,
           cell_size = 0.5)

# Individual clusters were annotated based on expression of known marker genes.

# Example
colData(acr2.lin39)$Cell.type <- "Unannotated"

# VA neurons are characterized by the co-expression of bnc-1 and unc-4
plot.expr.UMAP(acr2.lin39, "bnc-1", coexpr_gene = "unc-4", size = 0.5)

colData(acr2.lin39)$Cell.type <- ifelse(
  colData(acr2.lin39)$cluster == #cluster number,
    "VA",
  as.character(colData(acr2.lin39)$Cell.type)
)

# This was done iteratively to label as many clusters as possible.

# Subsetting VNC MNs (AS, DA, DB, DD, VA, VB, VC, VD) based on cell IDS

al.MN.wt <- acr2.lin39[,colData(acr2.lin39)$cluster %in% c()] #input cluster numbers of VNC MNs
al.MN.wt <- detect_genes(al.MN.wt)
al.MN.wt <- al.MN.wt[rowData(al.MN.wt)$num_cells_expressed > 5,]

al.MN.wt <- preprocess_cds(al.MN.wt, num_dim = 65)
plot_pc_variance_explained(al.MN.wt)

al.MN.wt <- align_cds(al.MN.wt, alignment_group = "Strain", alignment_k = 5)

al.MN.wt <- reduce_dimension(al.MN.wt, 
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.2,
                             umap.n_neighbors = 50)

al.MN.wt <- cluster_cells(al.MN.wt,
                          res = 3e-3)

colData(al.MN.wt)$UMAP_1 <- reducedDims(al.MN.wt)[["UMAP"]][,1]
colData(al.MN.wt)$UMAP_2 <- reducedDims(al.MN.wt)[["UMAP"]][,2]

plot_cells(al.MN.wt, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 0.5)

# Cluster IDs added to the colData metadata table for later subsetting.
colData(al.MN.wt)$cluster <- monocle3::clusters(al.MN.wt)

plot_cells(al.MN.wt,
           color_cells_by = "Cell.type",
           group_label_size = 3,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(al.MN.wt,
           color_cells_by = "strain",
           label_cell_groups = F,
           cell_size = 0.5)

# Some clusters were of low quality, they can be removed as illustrated below, followed by 
# subsequent re-running of the dimensionality reductions.

al.MN.wt <- al.MN.wt[,!(colData(al.MN.wt)$cluster %in% c())] # Input cluster numbers of low quality cells

al.MN.wt <- detect_genes(al.MN.wt)
al.MN.wt <- al.MN.wt[rowData(al.MN.wt)$num_cells_expressed > 5,]
al.MN.wt
# 10695 features in 13926 cells

al.MN.wt <- preprocess_cds(al.MN.wt, num_dim = 40)
plot_pc_variance_explained(al.MN.wt)

al.MN.wt <- align_cds(al.MN.wt, alignment_group = "Strain", alignment_k = 5)

al.MN.wt <- reduce_dimension(al.MN.wt, 
                             reduction_method = "UMAP",
                             preprocess_method = "Aligned",
                             umap.min_dist = 0.3,
                             umap.n_neighbors = 75)

al.MN.wt <- cluster_cells(al.MN.wt,
                          res = 3e-3)

colData(al.MN.wt)$UMAP_1 <- reducedDims(al.MN.wt)[["UMAP"]][,1]
colData(al.MN.wt)$UMAP_2 <- reducedDims(al.MN.wt)[["UMAP"]][,2]

plot_cells(al.MN.wt, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 0.5)

colData(al.MN.wt)$cluster <- monocle3::clusters(al.MN.wt)

plot_cells(al.MN.wt,
           color_cells_by = "Cell.type",
           group_label_size = 3,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(al.MN.wt,
           color_cells_by = "Strain",
           label_cell_groups = F,
           cell_size = 0.5)

plot.expr.UMAP(al.MN.wt, "mab-5", size = 0.5)

# Creating a subset of just the GABAergic MNs

GABA.MNs <- c() # clusters corresponding to GABA neurons (marked by expression of unc-25, unc-30, unc-47)

GABA.wt <- al.MN.wt[,colData(al.MN.wt)$cluster %in% GABA.MNs]
GABA.wt <- detect_genes(GABA.wt)
GABA.wt <- GABA.wt[rowData(GABA.wt)$num_cells_expressed > 5,]

GABA.wt <- preprocess_cds(GABA.wt, num_dim = 40)
plot_pc_variance_explained(GABA.wt)

GABA.wt <- reduce_dimension(GABA.wt, 
                            reduction_method = "UMAP",
                            preprocess_method = "PCA",
                            umap.min_dist = 0.2,
                            umap.n_neighbors = 50)

GABA.wt <- cluster_cells(GABA.wt,
                         res = 8e-3)

colData(GABA.wt)$UMAP_1 <- reducedDims(GABA.wt)[["UMAP"]][,1]
colData(GABA.wt)$UMAP_2 <- reducedDims(GABA.wt)[["UMAP"]][,2]

plot_cells(GABA.wt, 
           color_cells_by = "cluster",
           group_label_size = 3,
           cell_size = 0.5)

colData(GABA.wt)$cluster <- monocle3::clusters(GABA.wt)

plot_cells(GABA.wt,
           color_cells_by = "Cell.type",
           group_label_size = 3,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot.expr.UMAP(GABA.wt, "oig-1", size = 0.5)

# Identifying subcluster specific markers within GABA neurons
GABA.markers <- top_markers(GABA.wt, group_cells_by = "cluster", genes_to_test_per_group = 30)

# Looking at the top markers for clusters 20 and 19
GABA.markers %>% filter(cell_group == 20) %>% arrange(desc(marker_score)) %>% head(20)
GABA.markers %>% filter(cell_group == 19) %>% arrange(desc(marker_score)) %>% head(20)

plot.expr.UMAP(GABA.wt, "ceh-75", size = 0.5)

# Assigning subclass identity (iin this case to VD12)
colData(GABA.wt)$Cell.type <- ifelse(
  colData(GABA.wt)$cluster %in% c(19),
  "VD12",
  as.character(colData(GABA.wt)$Cell.type)
)

GABA.markers %>% filter(cell_group == 22) %>% arrange(desc(marker_score)) %>% head(20)
GABA.markers %>% filter(cell_group == 23) %>% arrange(desc(marker_score)) %>% head(20)
GABA.markers %>% filter(cell_group == 21) %>% arrange(desc(marker_score)) %>% head(20)

# Passing annotations from neuron subsets up to the larger datasets
colData(al.MN.wt)[colnames(GABA.wt),]$Cell.type <- colData(GABA.wt)$Cell.type
colData(al.MN)[colnames(al.MN.wt),]$Cell.type <- colData(al.MN.wt)$Cell.type

# Saving files
saveRDS(al.MN.wt, "./7595-ST/SoupX/091922_acr2_lin39_combined_wt_VNC_neurons_cds.rds")
saveRDS(al.MN, "./7595-ST/SoupX/091922_acr2_lin39_combined_all_genotypes_VNC_neurons_cds.rds")
saveRDS(GABA.wt, "./7595-ST/SoupX/091922_acr2_lin39_combined_GABA_wt_VNC_neurons_cds.rds")

# Example of subsetting the cholinergic MNs
ACh.wt <- al.MN.wt[,!colData(al.MN.wt)$cluster %in% GABA.MNs]
ACh.wt <- detect_genes(ACh.wt)
ACh.wt <- ACh.wt[rowData(ACh.wt)$num_cells_expressed > 5,]
ACh.wt
# 10115 features in 9666 cells

ACh.wt <- preprocess_cds(ACh.wt, num_dim = 30)
plot_pc_variance_explained(ACh.wt)

ACh.wt <- align_cds(ACh.wt, alignment_group = "Strain", alignment_k = 5)

ACh.wt <- reduce_dimension(ACh.wt, 
                           reduction_method = "UMAP",
                           preprocess_method = "Aligned",
                           umap.min_dist = 0.2,
                           umap.n_neighbors = 50)

ACh.wt <- cluster_cells(ACh.wt,
                        res = 8e-3)

colData(ACh.wt)$UMAP_1 <- reducedDims(ACh.wt)[["UMAP"]][,1]
colData(ACh.wt)$UMAP_2 <- reducedDims(ACh.wt)[["UMAP"]][,2]

plot_cells(ACh.wt, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

colData(ACh.wt)$cluster <- monocle3::clusters(ACh.wt)

plot_cells(ACh.wt,
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot.expr.UMAP(ACh.wt, "oig-1", size = 0.5)

ACh.markers <- top_markers(ACh.wt, group_cells_by = "cluster", genes_to_test_per_group = 30)
ACh.markers %>% filter(cell_group == 20) %>% arrange(desc(specificity)) %>% head(20)

ggplot(as.data.frame(colData(ACh.wt)), aes(x = cluster, y = num_genes_expressed)) + 
  geom_jitter(width = 0.05, height = 0.05) +
  geom_boxplot(outlier.shape = NA) + 
  theme_classic()

# Examples of annotating subclasses
colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(58, 39),
  "VC_4_5",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(35, 38, 18, 19, 25, 21, 29, 22, 20),
  "VC",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(1, 52, 42),
  "VB02",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(54, 14, 13),
  "VB01",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(10),
  "VB01.UPR",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(24),
  "DB02",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(33,11),
  "DB",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(17,59),
  "DB01",
  as.character(colData(ACh.wt)$Cell.type)
)

plot_cells(ACh.wt,
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(8, 53),
  "AS11",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(44),
  "VA01",
  as.character(colData(ACh.wt)$Cell.type)
)

colData(ACh.wt)$Cell.type <- ifelse(
  colData(ACh.wt)$cluster %in% c(56),
  "SAB",
  as.character(colData(ACh.wt)$Cell.type)
)

# Saving files
saveRDS(acr2.lin39, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/7595-ST/SoupX/092022_acr2_lin39_all_cells_cds.rds")
saveRDS(ACh.wt, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/7595-ST/SoupX/092022_acr2_lin39_combined_ACh_wt_VNC_neurons_cds.rds")
saveRDS(all.MN.wt, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/7595-ST/SoupX/091922_acr2_lin39_combined_wt_VNC_neurons_cds.rds")
saveRDS(GABA.wt, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/7595-ST/SoupX/091922_acr2_lin39_combined_GABA_wt_VNC_neurons_cds.rds")

