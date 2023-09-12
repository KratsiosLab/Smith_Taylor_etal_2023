# Examining mouse spinal cord scRNA-seq

library(Seurat)
library(dplyr)

library(ggplot2)

msp <- readRDS("~/Downloads/lepichon_gitler_integrated.RDS")
msp
colnames(msp@meta.data)
table(msp$cholinergictypes, exclude = NULL)
# LePichon labels
table(msp$group, exclude = NULL)
# lab
table(msp$cell_class, exclude = NULL)
# Gitler labels

# Unifying labels
msp@meta.data$cell_class <- ifelse(
  is.na(msp@meta.data$cell_class) & msp@meta.data$cholinergictypes == "Cholinergic Interneurons",
  "Cholinergic INs",
  ifelse(
    is.na(msp@meta.data$cell_class) & msp@meta.data$cholinergictypes == "Skeletal Motor Neurons",
    "Skeletal MNs",
    ifelse(
      is.na(msp@meta.data$cell_class) & msp@meta.data$cholinergictypes == "Visceral Motor Neurons",
      "Visceral MNs",
      as.character(msp@meta.data$cell_class)
    )
  )
)

table(msp@meta.data$cell_class, msp@meta.data$group,
        exclude = NULL)

# The subtypes may be delineated by the Seurat clusters
table(msp@meta.data$seurat_clusters, msp@meta.data$cell_class, exclude = NULL)

# This is true, but there are only 4 main clusters that contain Skeletal MNs.
DimPlot(msp, group.by = "cell_class")
FeaturePlot(msp, features = c("Tns1", "Bcl6", "Pax2"))

mMN <- subset(msp, subset = cell_class != "Cholinergic INs")
remove(msp)
mMN <- NormalizeData(mMN)
mMN <- ScaleData(mMN)
mMN <- FindVariableFeatures(mMN)
mMN <- RunPCA(mMN, npcs = 15, verbose = F)
mMN <- RunUMAP(mMN, reduction = "pca", dims = 1:15)
mMN <- FindNeighbors(mMN, reduction = "pca", dims = 1:15)
mMN <- FindClusters(mMN, resolution = 0.5)

DimPlot(mMN, group.by = "group")
FeaturePlot(mMN, features = c("Bcl6", "Tns1", "Nos1", "Fbn2"))
DimPlot(mMN, group.by = "cell_class")

v.MN <- subset(mMN, subset = cell_class == "Visceral MNs")

# Integrating at this level

# split the dataset into a list of two seurat objects (stim and CTRL)
v.list <- SplitObject(v.MN, split.by = "group")

# normalize and identify variable features for each dataset independently
v.list <- lapply(X = v.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = v.list)

v.anchors <- FindIntegrationAnchors(object.list = v.list, anchor.features = features)
v.MN <- IntegrateData(anchorset = v.anchors)
DefaultAssay(v.MN) <- "integrated"
v.git <- subset(v.MN, subset = group == "Gitler")
DefaultAssay(v.git) <- "RNA"

v.git <- ScaleData(v.git)
v.git <- FindVariableFeatures(v.git)
v.git <- RunPCA(v.git, npscs = 30)
ElbowPlot(v.git, ndims = 30)
v.git <- RunUMAP(v.git, reduction = "pca", dims=c(1:30), seed.use=5, n.neighbors = 40L)
v.git <- FindNeighbors(v.git, reduction = "pca", dims = 1:30)
v.git <- FindClusters(v.git, resolution = 0.5)

DimPlot(v.git, label = T)
# 18 clusters identified. How well do they align with expression of markers from Blum, et al?

DotPlot(v.git, features = rev(c("Trhr", "Nts", "Gm10248", "Chodl", "Mdga1", "Frem1", 
                           "Slc26a7", "Sst", "Bnc2", "Fras1", "Gpc3", "Slc16a12","Gas7", "Lypd6",
                           "Tll2", "Col24a1",
                           "Etv1", "Lgr5", "Mamdc2", "Gm20754", "Gm26911", "Gm43948", "Col12a1",
                           "Fn1", "Mpped2", "Dscaml1", "Fstl4", "Sox5", "Cdh20", "Pde11a",
                           "Tnc"))) + coord_flip()

# Can I manually merge clusters?

v.git$git_clusters <- paste("g", test$seurat_clusters, sep = "")

v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g17" = "g3")
v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g16" = "g11")
v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g12" = "g4")

v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g11" = "go5")
v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g6" = "go12")
v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g14" = "go13")
v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g8" = "go11")
v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g5" = "go6")
v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g13" = "go8")
v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "g15" = "go8")

v.git$UMAP_1 <- v.git@reductions$umap@cell.embeddings[,1]
v.git$UMAP_2 <- v.git@reductions$umap@cell.embeddings[,2]

DimPlot(v.git, group.by = "git_clusters", label = T) + geom_hline(yintercept = 0.5)

v.git$git_clusters <- ifelse(
  v.git$git_clusters == "g4" & v.git$UMAP_2 > 0,
  "g14",
  as.character(v.git$git_clusters)
)

table(v.git$git_clusters)

v.git$git_clusters <- dplyr::recode(v.git$git_clusters, "go5" = "g5",
                                    "go6" = "g6",
                                    "go8" = "g8",
                                    "go11" = "g11",
                                    "go12" = "g12", 
                                    "go13" = "g13")
table(v.git$git_clusters)

v.git$git_clusters <- factor(v.git$git_clusters, levels = c("g0", "g1", "g2", "g3", "g4",
                                                             "g5", "g6", "g7", "g8", "g9", 
                                                             "g10", "g11", "g12", "g13", "g14"))

FeaturePlot(v.git, c("Sox5", "Cdh20", "Fras1", "Gas7", "Lypd6"))

DimPlot(v.git, group.by = "git_clusters", label = T)
DimPlot(v.git, label = T)

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(c("Trhr", "Nts", "Gm10248", "Chodl", "Mdga1", "Frem1", 
                                "Slc26a7", "Sst", "Bnc2", "Fras1", "Gpc3", "Slc16a12","Gas7", "Lypd6",
                                "Tll2", "Col24a1",
                                "Etv1", "Lgr5", "Mamdc2", "Gm20754", "Gm26911", "Gm43948", "Col12a1",
                                "Fn1", "Mpped2", "Dscaml1", "Fstl4", "Pde11a", "Sox5", "Cdh20",
                                "Tnc"))) + coord_flip()

# the git_cluster annotations are quite close to what is shown in the paper. I'm not sure about 
# cluster 15

# Checking HD TFs

hd <- read.table("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/HD TFs Mouse AnimalTFDB3.csv",
                 header = T, sep = ",", stringsAsFactors = F)
head(hd)
hd.mouse <- hd$Column2
hd.mouse.use <- intersect(hd.mouse, rownames(v.git))
# 157 of the 240 are detected

library(pheatmap)
v.git.mat <- v.git@assays$RNA@counts
dim(v.git.mat)
hd.v.git.mat <- v.git.mat[hd.mouse.use, ]
hd.v.git.mat <- as.matrix(hd.v.git.mat)

# Some of these may not be expressed.
hd.expr <- data.frame(gene = hd.mouse.use,
                      expr = Matrix::rowSums(hd.v.git.mat),
                      n = Matrix::rowSums(hd.v.git.mat!=0))
table(hd.expr$expr == 0)
table(hd.expr$n == 0)

head(hd.expr)

DotPlot(v.git, group.by = "git_clusters", 
        features = hd.mouse.use, cluster.idents = F) + coord_flip()

test <- DotPlot(v.git, group.by = "git_clusters", 
                features = hd.mouse.use, cluster.idents = F) + coord_flip()

pct.hd <- test$data$pct.exp

head(pct.hd)

pct.hd <- matrix(pct.hd, nrow = 157)
pct.hd[1:5,1:5]
colnames(pct.hd) <- c("g0", "g1", 'g2', "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10", "g11", "g12",
                      "g13", "g14")

rownames(pct.hd) <- hd.mouse.use
avg.data <- data.frame(gene = hd.mouse.use,
                       max.pct = apply(pct.hd, 1, max),
                       min.pct = apply(pct.hd, 1, min))

head(avg.data)
table(avg.data$max.pct > 5)

# 49 of the genes are not found in more than 5% of the cells in any cluster.

hd.use <- avg.data %>% filter(max.pct > 5) %>% pull(gene)

ggplot(hd.expr, aes(x = n, y = expr)) + geom_point() 

ggplot(hd.expr[which(hd.expr$gene %in% hd.use),], aes(x = n, y = expr)) + geom_point() 

hd.v.git.mat <- hd.v.git.mat[hd.use,]
hd.v.git.s <- subset(v.git, features = hd.use)
hd.v.git.s
agg.hd.v.git.mat <- AverageExpression(hd.v.git.s, assays = "RNA", group.by = "git_clusters")
agg.hd.v.git.mat <- agg.hd.v.git.mat[[1]]
agg.hd.v.git.mat[1:5,1:5]
agg.scaled <- t(scale(t(agg.hd.v.git.mat)))

hd.avg.heat <- pheatmap(agg.hd.v.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                 treeheight_row = 100, treeheight_col = 50, 
                 color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

#library(viridis)
#hd.avg.heat <- pheatmap(agg.hd.v.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
#                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
#                    treeheight_row = 30, treeheight_col = 30, 
#                    color = inferno(10), border_color = NA)

hd.avg.scaled.heat <- pheatmap(agg.scaled, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                    treeheight_row = 30, treeheight_col = 30, 
                    color = colorRampPalette(colors = c("white", "red", "darkred"))(100), border_color = NA)

gene.order <- hd.avg.heat[["tree_row"]]$order
hd.order <- hd.use[gene.order]

DotPlot(v.git, group.by = "git_clusters", 
        features = hd.order, cluster.idents = F) + coord_flip()

# ordering by expression level
hd.expr <- hd.expr %>% arrange(desc(expr))
head(hd.expr)

new.hd.order <- hd.expr %>% filter(gene %in% hd.use) %>% pull(gene)

v.git$git_clusters <- factor(v.git$git_clusters, levels = c("g0", "g1", "g2", "g3", "g4",
                                                            "g5", "g6", "g7", "g8", "g9", 
                                                            "g10", "g11", "g12", "g13", "g14"))

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(hd.order), cluster.idents = T, dot.min = 0.05) + coord_flip()

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(hd.order), cluster.idents = F, dot.min = 0.05) + coord_flip()

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(new.hd.order), cluster.idents = F, dot.min = 0.05) + coord_flip()

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(new.hd.order), cluster.idents = T, dot.min = 0.05) + coord_flip()

# NHR 
nhr <- read.table("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/NHRs mouse JAX.csv",
                  sep = ",", header = T, stringsAsFactors = F)

nhr.gene <- unique(nhr$Symbol)
# 53 genes
nhr.mouse.use <- intersect(nhr.gene, rownames(v.git.mat))

nhr.v.git.mat <- v.git.mat[nhr.mouse.use, ]
nhr.v.git.mat <- as.matrix(nhr.v.git.mat)

# Some of these may not be expressed.
nhr.expr <- data.frame(gene = nhr.mouse.use,
                      expr = Matrix::rowSums(nhr.v.git.mat),
                      n = Matrix::rowSums(nhr.v.git.mat!=0))
table(nhr.expr$expr == 0)
table(nhr.expr$n == 0)

head(nhr.expr)
DotPlot(v.git, group.by = "git_clusters", 
        features = nhr.mouse.use, cluster.idents = F) + coord_flip()

nhr.test <- DotPlot(v.git, group.by = "git_clusters", 
                features = nhr.mouse.use, cluster.idents = F) + coord_flip()

pct.nhr <- nhr.test$data$pct.exp

head(pct.nhr)

pct.nhr <- matrix(pct.nhr, nrow = 53)
pct.nhr[1:5,1:5]
colnames(pct.nhr) <- c("g0", "g1", 'g2', "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10", "g11", "g12",
                      "g13", "g14")

rownames(pct.nhr) <- nhr.mouse.use
nhr.avg.data <- data.frame(gene = nhr.mouse.use,
                       max.pct = apply(pct.nhr, 1, max),
                       min.pct = apply(pct.nhr, 1, min))

head(nhr.avg.data)
table(nhr.avg.data$max.pct > 5)

# 37 of the genes are not found in more than 5% of the cells in any cluster.

nhr.use <- nhr.avg.data %>% filter(max.pct > 5) %>% pull(gene)

ggplot(nhr.expr, aes(x = n, y = expr)) + geom_point() 

ggplot(nhr.expr[which(nhr.expr$gene %in% nhr.use),], aes(x = n, y = expr)) + geom_point() 

nhr.v.git.mat <- nhr.v.git.mat[nhr.use,]
nhr.v.git.s <- subset(v.git, features = nhr.use)
nhr.v.git.s
agg.nhr.v.git.mat <- AverageExpression(nhr.v.git.s, assays = "RNA", group.by = "git_clusters")
agg.nhr.v.git.mat <- agg.nhr.v.git.mat[[1]]
agg.nhr.v.git.mat[1:5,1:5]
agg.nhr.scaled <- t(scale(t(agg.nhr.v.git.mat)))

nhr.avg.heat <- pheatmap(agg.nhr.v.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                        clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                        treeheight_row = 100, treeheight_col = 50, 
                        color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

#library(viridis)
#nhr.avg.heat <- pheatmap(agg.nhr.v.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
#                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
#                    treeheight_row = 30, treeheight_col = 30, 
#                    color = inferno(10), border_color = NA)

nhr.avg.scaled.heat <- pheatmap(agg.scaled, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                               clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                               treeheight_row = 30, treeheight_col = 30, 
                               color = colorRampPalette(colors = c("white", "red", "darkred"))(100), border_color = NA)

nhr.gene.order <- nhr.avg.heat[["tree_row"]]$order
nhr.order <- nhr.use[nhr.gene.order]

# ordering by expression level
nhr.expr <- nhr.expr %>% arrange(desc(expr))
head(nhr.expr)

new.nhr.order <- nhr.expr %>% filter(gene %in% nhr.use) %>% pull(gene)

v.git$git_clusters <- factor(v.git$git_clusters, levels = c("g0", "g1", "g2", "g3", "g4",
                                                            "g5", "g6", "g7", "g8", "g9", 
                                                            "g10", "g11", "g12", "g13", "g14"))

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(nhr.order), cluster.idents = T, dot.min = 0.05) + coord_flip()

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(nhr.order), cluster.idents = F, dot.min = 0.05) + coord_flip()

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(new.nhr.order), cluster.idents = F, dot.min = 0.05) + coord_flip()

DotPlot(v.git, group.by = "git_clusters", 
        features = rev(new.nhr.order), cluster.idents = T, dot.min = 0.05) + coord_flip()

saveRDS(v.git, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/072723_gitler_visceral_MN_subset.rds")

# passing back annotations to supersets
v.git$git_clusters <- as.character(v.git$git_clusters)
table(v.git$git_clusters)

v.MN$git_clusters <- "Unknown"
v.MN@meta.data[colnames(v.git),]$git_clusters <- v.git$git_clusters
table(v.MN$git_clusters)

mMN$v_git_clusters <- "Unknown"
mMN@meta.data[colnames(v.git),]$v_git_clusters <- v.git$git_clusters
table(mMN$v_git_clusters)

saveRDS(v.MN, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/072823_mouse_spinal_cord_visceral_integrated_scRNA.rds")
saveRDS(mMN, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/072823_mouse_spinal_cord_integrated_scRNA.rds")

# Now for skeletal MNs

sk.MN <- subset(mMN, subset = cell_class == "Skeletal MNs")
sk.git.MN <- subset(sk.MN, subset = group == "Gitler")

sk.git.MN <- ScaleData(sk.git.MN)
sk.git.MN <- FindVariableFeatures(sk.git.MN)
sk.git.MN <- RunPCA(sk.git.MN, npcs = 30)
sk.git.MN <- RunUMAP(sk.git.MN, reduction = "pca", dims = 1:20)
sk.git.MN <- FindNeighbors(sk.git.MN, reduction = "pca", dims = 1:20)
sk.git.MN <- FindClusters(sk.git.MN, resolution = 0.5)

DimPlot(sk.git.MN, label = T)

# 15 clusters identified

FeaturePlot(sk.git.MN, features = c("Htr1d", "Npas1", "Rbfox3", "Vipr2"))
FeaturePlot(sk.git.MN, features = c("Creb5", "Pard3b", "Kitl", "Stxbp6"))

DotPlot(sk.git.MN, features = c("Htr1d", "Npas1", "Creb5", "Pard3b", "Klhl1", "Sema3c", "Nrp2",
                            "Stxbp6", "Plch1", "Kitl", "Crtac1", "Rbfox3", "Vipr2", "Col5a1",
                            "Stk32a", "Sv2b"))

sk.git.MN$git_clusters <- "Unknown"
sk.git.MN$git_clusters <- ifelse(
  sk.git.MN@meta.data$seurat_cluster %in% c(1,7,10),
  "gamma",
  ifelse(
    sk.git.MN@meta.data$seurat_cluster %in% c(2),
    "gamma*",
    "alpha"
  )
)

DimPlot(sk.git.MN, group.by = "git_clusters", label = T)

mMN@meta.data[colnames(sk.git.MN),]$git_clusters <- sk.git.MN$git_clusters
table(mMN$git_clusters)
table(mMN$group)

a.clusters <- c(0,3,4,5,6,8,9,11,12,13,14)
a.git.MN <- subset(sk.git.MN, subset = seurat_clusters %in% a.clusters)
a.git.MN <- RunPCA(a.git.MN, npcs = 30)
ElbowPlot(a.git.MN, ndims = 30)
a.git.MN <- RunUMAP(a.git.MN, reduction = "pca", dims = 1:20)
a.git.MN <- FindNeighbors(a.git.MN, reduction = "pca", dims = 1:20)
a.git.MN <- FindClusters(a.git.MN, resolution = 0.5)
# 14 clusters. The publication has 12

DimPlot(a.git.MN, label = T)

DotPlot(a.git.MN, features = rev(c("Rreb1", "Sdk2", "Htr1f", "Pdgfd", "Ptchd4", "Etl4", "Grm5", 
                                   "Cacna1e", "Kcnh7", "Cpne4", "Grik1", "Hcrtr2", "Megf11",
                                   "Kcnq5", "Chodl", "Erbb4", "Gm2164", "Cdh9", "Gab1", "Prom1",
                                   "Clca3a1", "Cdh8", "A330008L17Rik", "Sema3e", "Ppp1r1c",
                                   "Mpped2", "Shisa9", "Hmcn1", "Slc44a5"))) + coord_flip()

a.git.MN$git_clusters <- ifelse(
  a.git.MN$seurat_clusters == 0,
  "sk0",
  ifelse(a.git.MN$seurat_clusters == 1,
         "sk1",
         ifelse(
           a.git.MN$seurat_clusters == 2,
           "sk2",
           ifelse(
             a.git.MN$seurat_clusters == 3,
             "sk3",
             ifelse(
               a.git.MN$seurat_clusters == 7,
               "sk4",
               ifelse(
                 a.git.MN$seurat_clusters == 4,
                 "sk5",
                 ifelse(
                   a.git.MN$seurat_clusters == 5,
                   "sk6",
                   ifelse(
                     a.git.MN$seurat_clusters == 6,
                     "sk7",
                     ifelse(
                       a.git.MN$seurat_clusters == 8,
                       "sk0",
                       ifelse(
                         a.git.MN$seurat_clusters == 10,
                         "sk8",
                         ifelse(
                           a.git.MN$seurat_clusters == 9,
                           "sk9",
                           ifelse(
                             a.git.MN$seurat_clusters == 11,
                             "sk10",
                             ifelse(
                               a.git.MN$seurat_clusters == 12,
                               "sk11",
                               ifelse(
                                 a.git.MN$seurat_clusters == 13,
                                 "sk8",
                                 as.character(a.git.MN$git_clusters)
                               )
                             )
                           )
                         )
                       )
                     )
                   )
                 )
               )
             )
           )
         )
  )
)

table(a.git.MN$git_clusters)

a.git.MN$git_clusters <- factor(a.git.MN$git_clusters, 
                                levels = c("sk0", "sk1", "sk2", "sk3", "sk4", "sk5", "sk6",
                                           "sk7", "sk8", "sk9", "sk10", "sk11"))

DotPlot(a.git.MN, group.by = "git_clusters",
        features = rev(c("Rreb1", "Sdk2", "Htr1f", "Pdgfd", "Ptchd4", "Etl4", "Grm5", 
                                   "Cacna1e", "Kcnh7", "Cpne4", "Grik1", "Hcrtr2", "Megf11",
                                   "Kcnq5", "Chodl", "Erbb4", "Gm2164", "Cdh9", "Gab1", "Prom1",
                                   "Clca3a1", "Cdh8", "A330008L17Rik", "Sema3e", "Ppp1r1c",
                                   "Mpped2", "Shisa9", "Hmcn1", "Slc44a5"))) + coord_flip()

DimPlot(a.git.MN, group.by = "git_clusters", label = T)

a.git.MN$git_clusters <- as.character(a.git.MN$git_clusters)

sk.git.MN@meta.data[colnames(a.git.MN),]$git_clusters <- a.git.MN$git_clusters
mMN@meta.data[colnames(a.git.MN),]$git_clusters <- a.git.MN$git_clusters
table(mMN$git_clusters, exclude = NULL)

a.git.MN$git_clusters <- factor(a.git.MN$git_clusters, 
                                levels = c("sk0", "sk1", "sk2", "sk3", "sk4", "sk5", "sk6",
                                           "sk7", "sk8", "sk9", "sk10", "sk11"))

# Now to look at HD and NHRs in skeletal muscle groups, both alpha and gamma

sk.hd.mouse.use <- intersect(hd.mouse, rownames(sk.git.MN))
# 157 of the 240 are detected

sk.git.mat <- sk.git.MN@assays$RNA@counts
dim(sk.git.mat)
hd.sk.git.mat <- sk.git.mat[sk.hd.mouse.use, ]
hd.sk.git.mat <- as.matrix(hd.sk.git.mat)

# Some of these may not be expressed.
hd.sk.expr <- data.frame(gene = sk.hd.mouse.use,
                      expr = Matrix::rowSums(hd.sk.git.mat),
                      n = Matrix::rowSums(hd.sk.git.mat!=0))
table(hd.sk.expr$expr == 0)
table(hd.sk.expr$n == 0)
# 141 are detected

head(hd.sk.expr)

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = sk.hd.mouse.use, cluster.idents = F) + coord_flip()

sk.test <- DotPlot(sk.git.MN, group.by = "git_clusters", 
                features = sk.hd.mouse.use, cluster.idents = F) + coord_flip()

sk.pct.hd <- sk.test$data$pct.exp

head(sk.pct.hd)

sk.pct.hd <- matrix(sk.pct.hd, nrow = 157)
sk.pct.hd[1:5,1:5]
colnames(sk.pct.hd) <- c("gamma", "gamma*", "sk0", "sk1", 'sk2', "sk3", "sk4", "sk5", 
                         "sk6", "sk7", "sk8", "sk9", "sk10", "sk11")

rownames(sk.pct.hd) <- sk.hd.mouse.use
sk.avg.data <- data.frame(gene = sk.hd.mouse.use,
                       max.pct = apply(sk.pct.hd, 1, max),
                       min.pct = apply(sk.pct.hd, 1, min))

head(sk.avg.data)
table(sk.avg.data$max.pct > 5)

# 72 of the genes are not found in more than 5% of the cells in any cluster.

hd.sk.use <- sk.avg.data %>% filter(max.pct > 5) %>% pull(gene)

ggplot(hd.sk.expr, aes(x = n, y = expr)) + geom_point() 

ggplot(hd.sk.expr[which(hd.sk.expr$gene %in% hd.sk.use),], aes(x = n, y = expr)) + geom_point() 

hd.sk.git.mat <- hd.sk.git.mat[hd.sk.use,]
hd.sk.git.s <- subset(sk.git.MN, features = hd.sk.use)
hd.sk.git.s
agg.hd.sk.git.mat <- AverageExpression(hd.sk.git.s, assays = "RNA", group.by = "git_clusters")
agg.hd.sk.git.mat <- agg.hd.sk.git.mat[[1]]
agg.hd.sk.git.mat[1:5,1:5]
agg.hd.sk.scaled <- t(scale(t(agg.hd.sk.git.mat)))

hd.sk.avg.heat <- pheatmap(agg.hd.sk.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                        clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                        treeheight_row = 100, treeheight_col = 50, 
                        color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

#library(viridis)
#hd.sk.avg.heat <- pheatmap(agg.hd.sk.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
#                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
#                    treeheight_row = 30, treeheight_col = 30, 
#                    color = inferno(10), border_color = NA)

# hd.sk.avg.scaled.heat <- pheatmap(agg.hd.sk.scaled, cluster_rows = T, cluster_cols = T, show_colnames = T, 
#                                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
#                                treeheight_row = 30, treeheight_col = 30, 
#                                color = colorRampPalette(colors = c("white", "red", "darkred"))(100), border_color = NA)

gene.order <- hd.sk.avg.heat[["tree_row"]]$order
hd.sk.order <- hd.sk.use[gene.order]

# ordering by expression level
hd.sk.expr <- hd.sk.expr %>% arrange(desc(expr))
head(hd.sk.expr)

new.hd.sk.order <- hd.sk.expr %>% filter(gene %in% hd.sk.use) %>% pull(gene)

sk.git.MN$git_clusters <- factor(sk.git.MN$git_clusters, 
                                 levels = c("gamma", "gamma*", "sk0", "sk1", 'sk2', "sk3", "sk4", "sk5", 
                                           "sk6", "sk7", "sk8", "sk9", "sk10", "sk11"))

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = rev(hd.sk.order), cluster.idents = T, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = rev(hd.sk.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = rev(new.hd.sk.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = rev(new.hd.sk.order), cluster.idents = T, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

# NHR 
nhr.mouse.use <- intersect(nhr.gene, rownames(sk.git.mat))

nhr.sk.git.mat <- sk.git.mat[nhr.mouse.use, ]
nhr.sk.git.mat <- as.matrix(nhr.sk.git.mat)

# Some of these may not be expressed.
sk.nhr.expr <- data.frame(gene = nhr.mouse.use,
                       expr = Matrix::rowSums(nhr.sk.git.mat),
                       n = Matrix::rowSums(nhr.sk.git.mat!=0))
table(sk.nhr.expr$expr == 0)
table(sk.nhr.expr$n == 0)

head(sk.nhr.expr)
DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = nhr.mouse.use, cluster.idents = F) + coord_flip()

sk.nhr.test <- DotPlot(sk.git.MN, group.by = "git_clusters", 
                    features = nhr.mouse.use, cluster.idents = F) + coord_flip()

sk.pct.nhr <- sk.nhr.test$data$pct.exp

head(sk.pct.nhr)

sk.pct.nhr <- matrix(sk.pct.nhr, nrow = 53)
sk.pct.nhr[1:5,1:5]
colnames(sk.pct.nhr) <- c("gamma", "gamma*", "sk0", "sk1", 'sk2', "sk3", "sk4", "sk5", 
                          "sk6", "sk7", "sk8", "sk9", "sk10", "sk11")

rownames(sk.pct.nhr) <- nhr.mouse.use
sk.nhr.avg.data <- data.frame(gene = nhr.mouse.use,
                           max.pct = apply(sk.pct.nhr, 1, max),
                           min.pct = apply(sk.pct.nhr, 1, min))

head(sk.nhr.avg.data)
table(sk.nhr.avg.data$max.pct > 5)

# 41 of the genes are good to keep.

sk.nhr.use <- sk.nhr.avg.data %>% filter(max.pct > 5) %>% pull(gene)

ggplot(sk.nhr.expr, aes(x = n, y = expr)) + geom_point() 

ggplot(sk.nhr.expr[which(sk.nhr.expr$gene %in% sk.nhr.use),], aes(x = n, y = expr)) + geom_point() 

nhr.sk.git.mat <- nhr.sk.git.mat[sk.nhr.use,]
nhr.sk.git.s <- subset(sk.git.MN, features = sk.nhr.use)
nhr.sk.git.s
agg.nhr.sk.git.mat <- AverageExpression(nhr.sk.git.s, assays = "RNA", group.by = "git_clusters")
agg.nhr.sk.git.mat <- agg.nhr.sk.git.mat[[1]]
agg.nhr.sk.git.mat[1:5,1:5]
agg.nhr.sk.scaled <- t(scale(t(agg.nhr.sk.git.mat)))

sk.nhr.avg.heat <- pheatmap(agg.nhr.sk.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                         treeheight_row = 100, treeheight_col = 50, 
                         color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

sk.nhr.avg.heat <- pheatmap(agg.nhr.sk.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                            clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                            treeheight_row = 100, treeheight_col = 50, scale = "row",
                            color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

#library(viridis)
#sk.nhr.avg.heat <- pheatmap(agg.nhr.sk.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
#                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
#                    treeheight_row = 30, treeheight_col = 30, 
#                    color = inferno(10), border_color = NA)

sk.nhr.avg.scaled.heat <- pheatmap(agg.nhr.sk.scaled, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                                treeheight_row = 30, treeheight_col = 30, 
                                color = colorRampPalette(colors = c("white", "red", "darkred"))(100), border_color = NA)

sk.nhr.gene.order <- sk.nhr.avg.heat[["tree_row"]]$order
sk.nhr.order <- sk.nhr.use[sk.nhr.gene.order]

# ordering by expression level
sk.nhr.expr <- sk.nhr.expr %>% arrange(desc(expr))
head(sk.nhr.expr)

new.sk.nhr.order <- sk.nhr.expr %>% filter(gene %in% sk.nhr.use) %>% pull(gene)

sk.git.MN$git_clusters <- factor(sk.git.MN$git_clusters, 
                                 levels = c("gamma", "gamma*", "sk0", "sk1", 'sk2', "sk3", "sk4", "sk5", 
                                            "sk6", "sk7", "sk8", "sk9", "sk10", "sk11"))

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = rev(sk.nhr.order), cluster.idents = T, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = rev(sk.nhr.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = rev(new.sk.nhr.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(sk.git.MN, group.by = "git_clusters", 
        features = rev(new.sk.nhr.order), cluster.idents = T, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

saveRDS(sk.git.MN, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/072723_gitler_visceral_MN_subset.rds")

all.git.MN <- subset(mMN, subset = group == "Gitler")

table(all.git.MN$git_clusters, exclude = NULL)

all.hd.mouse.use <- intersect(hd.mouse, rownames(all.git.MN))
# 157 of the 240 are detected

all.git.mat <- all.git.MN@assays$RNA@counts
dim(all.git.mat)
hd.all.git.mat <- all.git.mat[all.hd.mouse.use, ]
hd.all.git.mat <- as.matrix(hd.all.git.mat)

# Some of these may not be expressed.
hd.all.expr <- data.frame(gene = all.hd.mouse.use,
                         expr = Matrix::rowSums(hd.all.git.mat),
                         n = Matrix::rowSums(hd.all.git.mat!=0))
table(hd.all.expr$expr == 0)
table(hd.all.expr$n == 0)
# 151 are detected

head(hd.all.expr)

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = all.hd.mouse.use, cluster.idents = F) + coord_flip() + 
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

all.test <- DotPlot(all.git.MN, group.by = "git_clusters", 
                   features = all.hd.mouse.use, cluster.idents = F) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))
  
all.pct.hd <- all.test$data$pct.exp

head(all.pct.hd)

all.pct.hd <- matrix(all.pct.hd, nrow = 157)
dim(all.pct.hd)
all.pct.hd[1:5,1:5]
colnames(all.pct.hd) <- c("g0", "g1", "g10", "g11", "g12", "g13", "g14", "g2", "g3", "g4", "g5", "g6",
                          "g7", "g8", "g9", "gamma", "gamma*", "sk0", "sk1", "sk10", "sk11", "sk2",
                          "sk3", "sk4", "sk5", "sk6", "sk7", "sk8", "sk9") 

rownames(all.pct.hd) <- all.hd.mouse.use
all.avg.data <- data.frame(gene = all.hd.mouse.use,
                          max.pct = apply(all.pct.hd, 1, max),
                          min.pct = apply(all.pct.hd, 1, min))

head(all.avg.data)
table(all.avg.data$max.pct > 5)

# 76 of the genes are found in more than 5% of the cells in any cluster.

hd.all.use <- all.avg.data %>% filter(max.pct > 5) %>% pull(gene)

ggplot(hd.all.expr, aes(x = n, y = expr)) + geom_point() 

ggplot(hd.all.expr[which(hd.all.expr$gene %in% hd.all.use),], aes(x = n, y = expr)) + geom_point() 

hd.all.git.mat <- hd.all.git.mat[hd.all.use,]
hd.all.git.s <- subset(all.git.MN, features = hd.all.use)
hd.all.git.s
agg.hd.all.git.mat <- AverageExpression(hd.all.git.s, assays = "RNA", group.by = "git_clusters")
agg.hd.all.git.mat <- agg.hd.all.git.mat[[1]]
agg.hd.all.git.mat[1:5,1:5]

colnames(agg.hd.all.git.mat) <- c("g0", "g1", "g10", "g11", "g12", "g13", "g14", "g2", "g3", "g4", 
                                  "g5", "g6", "g7", "g8", "g9", "gamma", "gamma*", "sk0",
                                  "sk1", "sk10", "sk11", "sk2", "sk3", "sk4", "sk5", "sk6", "sk7", "sk8", "sk9")

agg.hd.all.scaled <- t(scale(t(agg.hd.all.git.mat)))

type.ann.df <- data.frame(row.names = c("g0", "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9",
                                       "g10", "g11", "g12", "g13", "g14", "gamma", "gamma*", "sk0",
                                       "sk1", "sk2", "sk3", "sk4", "sk5", "sk6", "sk7", "sk8", "sk9",
                                       "sk10", "sk11"), 
                          class = c(rep("Visceral", 15), rep("Skeletal", 14)))

type.ann.df

colnames(agg.hd.all.git.mat)

ann_colors = list(class = c(Visceral = "lightgreen", Skeletal = "steelblue"))

hd.all.avg.heat <- pheatmap(agg.hd.all.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                           treeheight_row = 100, treeheight_col = 50, annotation_col = type.ann.df,
                           annotation_colors = ann_colors, clustering_method = "ward.D2",
                           color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

#library(viridis)
#hd.all.avg.heat <- pheatmap(agg.hd.all.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
#                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
#                    treeheight_row = 30, treeheight_col = 30, 
#                    color = inferno(10), border_color = NA)

# hd.all.avg.scaled.heat <- pheatmap(agg.hd.all.scaled, cluster_rows = T, cluster_cols = T, show_colnames = T, 
#                                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
#                                treeheight_row = 30, treeheight_col = 30, 
#                                color = colorRampPalette(colors = c("white", "red", "darkred"))(100), border_color = NA)

gene.order <- hd.all.avg.heat[["tree_row"]]$order
hd.all.order <- hd.all.use[gene.order]

# ordering by expression level
hd.all.expr <- hd.all.expr %>% arrange(desc(expr))
head(hd.all.expr)

new.hd.all.order <- hd.all.expr %>% filter(gene %in% hd.all.use) %>% pull(gene)

hclust.hd.cell.order <- hd.all.avg.heat[["tree_col"]]$order
hd.cell.order <- colnames(agg.hd.all.git.mat)[hclust.hd.cell.order]

all.git.MN$git_clusters <- factor(all.git.MN$git_clusters, levels = hd.cell.order)

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(hd.all.order), cluster.idents = T, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(hd.all.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(new.hd.all.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(new.hd.all.order), cluster.idents = T, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

# NHR 
nhr.mouse.use <- intersect(nhr.gene, rownames(all.git.mat))

nhr.all.git.mat <- all.git.mat[nhr.mouse.use, ]
nhr.all.git.mat <- as.matrix(nhr.all.git.mat)

# Some of these may not be expressed.
all.nhr.expr <- data.frame(gene = nhr.mouse.use,
                          expr = Matrix::rowSums(nhr.all.git.mat),
                          n = Matrix::rowSums(nhr.all.git.mat!=0))
table(all.nhr.expr$expr == 0)
table(all.nhr.expr$n == 0)

head(all.nhr.expr)
DotPlot(all.git.MN, group.by = "git_clusters", 
        features = nhr.mouse.use, cluster.idents = F) + coord_flip()

all.nhr.test <- DotPlot(all.git.MN, group.by = "git_clusters", 
                       features = nhr.mouse.use, cluster.idents = F) + coord_flip()

all.pct.nhr <- all.nhr.test$data$pct.exp

head(all.pct.nhr)

all.pct.nhr <- matrix(all.pct.nhr, nrow = 53)
all.pct.nhr[1:5,1:5]

colnames(all.pct.nhr) <- c("g0", "g1", "g10", "g11", "g12", "g13", "g14", "g2", "g3", "g4", "g5", "g6",
                           "g7", "g8", "g9", "gamma", "gamma*", "sk0", "sk1", "sk10", "sk11", "sk2",
                           "sk3", "sk4", "sk5", "sk6", "sk7", "sk8", "sk9")

rownames(all.pct.nhr) <- nhr.mouse.use
all.nhr.avg.data <- data.frame(gene = nhr.mouse.use,
                              max.pct = apply(all.pct.nhr, 1, max),
                              min.pct = apply(all.pct.nhr, 1, min))

head(all.nhr.avg.data)
table(all.nhr.avg.data$max.pct > 5)

# 43 of the genes are good to keep.

all.nhr.use <- all.nhr.avg.data %>% filter(max.pct > 5) %>% pull(gene)

ggplot(all.nhr.expr, aes(x = n, y = expr)) + geom_point() 

ggplot(all.nhr.expr[which(all.nhr.expr$gene %in% all.nhr.use),], aes(x = n, y = expr)) + geom_point() 

nhr.all.git.mat <- nhr.all.git.mat[all.nhr.use,]
nhr.all.git.s <- subset(all.git.MN, features = all.nhr.use)
nhr.all.git.s
agg.nhr.all.git.mat <- AverageExpression(nhr.all.git.s, assays = "RNA", group.by = "git_clusters")
agg.nhr.all.git.mat <- agg.nhr.all.git.mat[[1]]
agg.nhr.all.git.mat[1:5,1:5]
agg.nhr.all.scaled <- t(scale(t(agg.nhr.all.git.mat)))

all.nhr.avg.heat <- pheatmap(agg.nhr.all.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                            clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                            treeheight_row = 150, treeheight_col = 100, clustering_method = "ward.D2",
                            annotation_col = type.ann.df, annotation_colors = ann_colors,
                            color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

#library(viridis)
#all.nhr.avg.heat <- pheatmap(agg.nhr.all.git.mat, cluster_rows = T, cluster_cols = T, show_colnames = T, 
#                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
#                    treeheight_row = 30, treeheight_col = 30, 
#                    color = inferno(10), border_color = NA)

all.nhr.avg.scaled.heat <- pheatmap(agg.nhr.all.scaled, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                                   clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                                   treeheight_row = 30, treeheight_col = 30, clustering_method = "ward.D2",
                                   color = colorRampPalette(colors = c("white", "red", "darkred"))(100), border_color = NA)

all.nhr.gene.order <- all.nhr.avg.heat[["tree_row"]]$order
all.nhr.order <- all.nhr.use[all.nhr.gene.order]

hclust.nhr.cell.order <- all.nhr.avg.heat[["tree_col"]]$order
nhr.cell.order <- colnames(agg.nhr.all.git.mat)[hclust.nhr.cell.order]

all.git.MN$git_clusters <- factor(all.git.MN$git_clusters, levels = nhr.cell.order)

# ordering by expression level
all.nhr.expr <- all.nhr.expr %>% arrange(desc(expr))
head(all.nhr.expr)

new.all.nhr.order <- all.nhr.expr %>% filter(gene %in% all.nhr.use) %>% pull(gene)

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(all.nhr.order), cluster.idents = T, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

all.git.MN$git_clusters <- factor(all.git.MN$git_clusters, levels = c("g0", "g1", "g2", "g3", "g4", 
                                                                      "g5", "g6", "g7", "g8", "g9","g10", "g11", "g12", "g13", "g14",
                                                                      "gamma", "gamma*", "sk0",
                                                                      "sk1", "sk2", "sk3", "sk4", "sk5", "sk6", "sk7", "sk8", "sk9",
                                                                      "sk10", "sk11"))

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(all.nhr.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(new.all.nhr.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(new.all.nhr.order), cluster.idents = T, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

saveRDS(all.git.MN, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/072823_mouse_spinal_cord_gitler_scRNA.rds")

saveRDS(mMN, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/072823_mouse_spinal_cord_integrated_scRNA.rds")

# I can try just clustering without a heatmap. Generating these to save


# clustering on all.git.MN with hclust for both HD and NHR families

agg.hd.all.git.mat

hd.clust <- hclust(dist(t(agg.hd.all.git.mat)), method = "ward.D2")
plot(hd.clust)
rect.hclust(hd.clust, k = 2)

agg.nhr.all.git.mat

nhr.clust <- hclust(dist(t(agg.nhr.all.git.mat)), method = "ward.D2")
plot(nhr.clust)

all.avg.expr.mat <- AverageExpression(all.git.MN, assays = "RNA", group.by = "git_clusters")
all.avg.expr.mat <- all.avg.expr.mat[[1]]
colnames(all.avg.expr.mat)
colnames(all.avg.expr.mat) <- as.character(unname(nhr.cell.order))

dim(all.avg.expr.mat)
all.gene.data <- data.frame(gene = rownames(all.avg.expr.mat),
                            expr = Matrix::rowSums(all.avg.expr.mat),
                            n = Matrix::rowSums(all.avg.expr.mat!=0))
table(all.gene.data$expr == 0)
table(all.gene.data$n == 0)
# 310 are not detected

all.genes.use <- all.gene.data %>% filter(expr > 0) %>% pull(gene)

all.test <- DotPlot(all.git.MN, group.by = "git_clusters", 
                        features = all.genes.use, cluster.idents = F) + coord_flip()

all.pct <- all.test$data$pct.exp

head(all.pct)

all.pct <- matrix(all.pct, nrow = 24788)
dim(all.pct)
all.pct[1:5,1:5]
table(all.git.MN$git_clusters)
colnames(all.pct) <- colnames(all.avg.expr.mat)

rownames(all.pct) <- all.genes.use
all.avg.data <- data.frame(gene = rownames(all.pct),
                               max.pct = apply(all.pct, 1, max),
                               min.pct = apply(all.pct, 1, min))

head(all.avg.data)
table(all.avg.data$max.pct > 5)

# 16468 of the genes are good to keep.

all.genes.use <- all.avg.data %>% filter(max.pct > 5) %>% pull(gene)
all.avg.expr.mat <- all.avg.expr.mat[all.genes.use,]
all.avg.heat <- pheatmap(all.avg.expr.mat, cluster_rows = F, cluster_cols = T, show_colnames = T, 
                         show_rownames = F, 
                             clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                             treeheight_row = 0, treeheight_col = 100, clustering_method = "ward.D2",
                             annotation_col = type.ann.df, annotation_colors = ann_colors,
                             color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

# What happens if I try to run dimensionality reduction solely on HD TFs?

hd.all.git.s <- RunPCA(hd.all.git.s, features = rownames(hd.all.git.s), npcs = 15)
ElbowPlot(hd.all.git.s)
hd.all.git.s <- RunUMAP(hd.all.git.s, dims = 1:10)
hd.all.git.s <- FindNeighbors(hd.all.git.s, dims = 1:10)
hd.all.git.s <- FindClusters(hd.all.git.s, resolution = 0.5)
# 11 clusters

DimPlot(hd.all.git.s, label = T)
DimPlot(hd.all.git.s, group.by = "git_clusters", label = T)
DimPlot(hd.all.git.s, group.by = "cell_class", label = T)


# NHR
nhr.all.git.s <- RunPCA(nhr.all.git.s, features = rownames(nhr.all.git.s), npcs = 15)
ElbowPlot(nhr.all.git.s)
nhr.all.git.s <- RunUMAP(nhr.all.git.s, dims = 1:10)
nhr.all.git.s <- FindNeighbors(nhr.all.git.s, dims = 1:10)
nhr.all.git.s <- FindClusters(nhr.all.git.s, resolution = 0.5)
# 10 clusters

DimPlot(nhr.all.git.s, label = T)
DimPlot(nhr.all.git.s, group.by = "git_clusters", label = T)
DimPlot(nhr.all.git.s, group.by = "cell_class", label = T)

# I could try using both HD and NHRs together

hd.nhr.all.git.s <- subset(all.git.MN, features = c(all.nhr.use, hd.all.use))
hd.nhr.all.git.s <- RunPCA(hd.nhr.all.git.s, features = rownames(hd.nhr.all.git.s), npcs = 15)
ElbowPlot(hd.nhr.all.git.s)
hd.nhr.all.git.s <- RunUMAP(hd.nhr.all.git.s, dims = 1:10)
hd.nhr.all.git.s <- FindNeighbors(hd.nhr.all.git.s, dims = 1:10)
hd.nhr.all.git.s <- FindClusters(hd.nhr.all.git.s, resolution = 0.5)
# 8 clusters

DimPlot(hd.nhr.all.git.s, label = T)
DimPlot(hd.nhr.all.git.s, group.by = "git_clusters", label = T)
DimPlot(hd.nhr.all.git.s, group.by = "cell_class", label = T)

# For the worm data, what happens if you try to run dimensionality reduction just off HD? NHR?, both?

MN.cds <- readRDS("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/060223_wt_MN_VNC_updated_subtype_names_cleaned_cds.rds")
library(monocle3)
library(wbData)
library(dplyr)
gids <- wb_load_gene_ids("WS273")

tf.jay <- read.table("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/Analysis/MN_TFs_FromGoAnalysis.csv",
                     sep = ",", header = T)

head(tf.jay)
tfs.jay <- tf.jay$Wormbase.ID
tfs.class <- tf.jay$Category.2 %>% strsplit(": ") %>% sapply("[", 2)
tf.df <- data.frame(id = tfs.jay,
                    gene_short_name = i2s(tfs.jay, gids),
                    class = tfs.class)

table(tf.df$class)
hd.use <- intersect(tf.df[which(tf.df$class == "homeodomain"),]$id, rownames(MN.cds))
i2s(hd.use, gids)

MN.t2 <- read.table("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/Files_for_website/072723_MN_subclass_thresholded_expression.csv",
                    header = T, row.names = 1, stringsAsFactors = F, sep = ",")
colnames(MN.t2)
MN.t2["ceh-9",]
hd.use <- intersect(tf.df[which(tf.df$class == "homeodomain"),]$gene_short_name, rownames(MN.t2))
MN.t2.df <- data.frame(gene = rownames(MN.t2),
                       n = Matrix::rowSums(MN.t2!=0))
table(MN.t2.df$n == 0)
hd.use <- MN.t2.df %>% filter(gene %in% hd.use & n > 0)
hd.use <- hd.use$gene
hd.use

library(Seurat)
rownames(MN.cds) <- rowData(MN.cds)$gene_short_name 

MN.s <- CreateSeuratObject(counts = exprs(MN.cds),
                           meta.data = as.data.frame(colData(MN.cds)))

MN.s <- NormalizeData(MN.s)
MN.s <- ScaleData(MN.s)

mn.avg <- AverageExpression(MN.s, group.by = "Cell.type")
mn.avg <- mn.avg[[1]]
table(MN.s$Cell.type)
colnames(mn.avg) <- c("AS11", "AS2_3", "AS4_8", "AS9_10", "DA1", "DA2_5", "DA6_8", "DA9", "DB1",
                      "DB2", "DB3_7", "DD2_3", "DD4_5", "VA1", "VA11", "VA12", "VA2", "VA3_8",
                      "VA9_10", "VB1", "VB10_11", "VB2", "VB3", "VB4_9", "VC1_3_6", "VC4_5", 
                      "VD12", "VD3_7", "VD8_11")
mn.avg[c("unc-4", "unc-129", "npr-2", "ceh-12", "vab-7"),]
mn.avg[c("mab-5", "lin-39", "egl-5", "ceh-13", "vab-7"),]

avg.hd <- intersect(hd.use, rownames(mn.avg))
length(intersect(hd.use, rownames(mn.avg)))
length(intersect(hd.use, rownames(MN.t2)))
# 39. Which are missing? 

setdiff(hd.use, avg.hd)

mn.hd.avg <- mn.avg[avg.hd,]

library(pheatmap)
ce.avg.hd.heat <- pheatmap(mn.hd.avg, cluster_rows = T, cluster_cols = T, show_colnames = T, 
                         show_rownames = T, 
                         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                         treeheight_row = 100, treeheight_col = 100, clustering_method = "ward.D2",
                         color = colorRampPalette(colors = c("white","red", "red", "darkred", "darkred"))(100), border_color = NA)

# what about just clustering?
MN.hd.s <- subset(MN.s, features = avg.hd)
MN.hd.s <- RunPCA(MN.hd.s, features = rownames(MN.s), npcs = 25)
ElbowPlot(MN.hd.s, ndims = 25)
MN.hd.s <- RunUMAP(MN.hd.s, dims = 1:15)
MN.hd.s <- FindNeighbors(MN.hd.s, dims = 1:15)
MN.hd.s <- FindClusters(MN.hd.s, resolution = 0.5)
# 32 clusters

DimPlot(MN.hd.s, label = T)
DimPlot(MN.hd.s, group.by = "Cell.type", label = T)
DimPlot(MN.hd.s, group.by = "anatomical_class", label = T)
DimPlot(MN.hd.s, group.by = "strain", label = T)

dim(mn.avg)
# 10695

MN.cds
dim(MN.t2)
# 10710. There is a difference in the number of genes

length(setdiff(rownames(mn.avg), rownames(MN.t2)))
# 1306

length(setdiff(rownames(MN.t2), rownames(mn.avg)))
# 1321

length(intersect(rownames(MN.t2), rowData(MN.cds)$gene_short_name))
# only 9390. Why is there a difference? There are genes in the thresholded data that are
# not present in the dataset. And there are genes in the dataset that are not present in the 
# thresholded data. I don't understand this difference. 

wt.MN <- readRDS("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/7595-ST/SoupX/091922_acr2_lin39_combined_wt_VNC_neurons_cds.rds")
wt.MN
# Also 10695. Maybe this is because the thresholding was performed on the larger dataset
# with additional neurons. It makes sense more genes would be present. But why would there be 
# genes in the data that were not included in the thresholding? Perhaps they are very lowly 
# expressed and get removed?

data.only <- setdiff(rowData(MN.cds)$gene_short_name, rownames(MN.t2))
raw.prop <- readRDS("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/Analysis/060823_YA_proportion_matrix.rds")
dim(raw.prop)
MN.prop <- raw.prop[, c("AS1_3", "AS4_8", "AS9_10", "AS11", "DA1", "DA2_5", "DA6_8", "DA9", 
                        "DB1", "DB2", "DB3_7", "DD2_3", "DD4_5", "VA1", "VA2", "VA3_8", "VA9_10",
                        "VA11", "VA12", "VB1", "VB2", "VB3", "VB4_9", "VB10_11", "VC1_3_6", 
                        "VC4_5", "VD3_7", "VD8_11", "VD12")]

MN.prop[s2i(data.only[5:20], gids),]

# These are probably genes that did not pass the threshold for being included (i.e., below 2% of 
# the cells in every cluster.) So they were removed from the thresholded data.
dim(MN.prop)
# 12595
dim(MN.t2)
# 10710

length(intersect(rownames(MN.prop), s2i(rownames(MN.t2), gids)))
# 10710. So all of the MN.t2 are present in the MN.prop

length(intersect(rownames(MN.prop), s2i(data.only, gids)))
# 1285 of the 1305. So there are some (20) of the data.only set that were likely too low even 
# among all cells.

test.max <- apply(MN.prop[intersect(rownames(MN.prop), s2i(data.only, gids)),], 1, max)
max(test.max)
# 0.0198915. Just as I suspected.

# now about the other dataset that are in the thresholded data but not in the dataset. 
thres.only <- setdiff(rownames(MN.t2), rowData(MN.cds)$gene_short_name)
test2.max <- apply(MN.prop[intersect(rownames(MN.prop), s2i(thres.only, gids)),], 1, max)
max(test2.max)
# 0.0526                           

summary(test2.max)

all.test2.max <- apply(raw.prop[intersect(rownames(raw.prop), s2i(thres.only, gids)), ], 1, max)
summary(all.test2.max)
# these genes are much more abundant in the whole dataset. They are likely not detected in the 
# MN.cds because of removing genes that are found in <5 cells when the MNs are subset. That 
# is my best guess. 

# 7/31/23 renaming clusters for the newer dotplots with different gene order.

all.git.MN@meta.data$git_clusters <- as.character(all.git.MN@meta.data$git_clusters)
all.git.MN@meta.data$git_clusters <- dplyr::recode(all.git.MN@meta.data$git_clusters,
                                                   "g0" = "v0",
                                                   "g1" = "v1",
                                                   "g2" = "v2",
                                                   "g3" = "v3",
                                                   "g4" = "v4",
                                                   "g5" = "v5", 
                                                   "g6" = "v6",
                                                   "g7" = "v7",
                                                   "g8" = "v8",
                                                   "g9" = "v9",
                                                   "g10" = "v10",
                                                   "g11" = "v11",
                                                   "g12" = "v12",
                                                   "g13" = "v13",
                                                   "g14" = "v14")

hd <- read.table("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/HD TFs Mouse AnimalTFDB3.csv",
                                header = T, sep = ",", stringsAsFactors = F)
hd.mouse <- hd$Column2

all.hd.mouse.use <- intersect(hd.mouse, rownames(all.git.MN))
# 157 of the 240 are detected

all.git.mat <- all.git.MN@assays$RNA@counts
dim(all.git.mat)
hd.all.git.mat <- all.git.mat[all.hd.mouse.use, ]
hd.all.git.mat <- as.matrix(hd.all.git.mat)

# Some of these may not be expressed.
hd.all.expr <- data.frame(gene = all.hd.mouse.use,
                          expr = Matrix::rowSums(hd.all.git.mat),
                          n = Matrix::rowSums(hd.all.git.mat!=0))
table(hd.all.expr$expr == 0)
table(hd.all.expr$n == 0)
# 151 are detected

head(hd.all.expr)
library(ggplot2)
DotPlot(all.git.MN, group.by = "git_clusters", 
        features = all.hd.mouse.use, cluster.idents = F) + coord_flip() + 
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

all.test <- DotPlot(all.git.MN, group.by = "git_clusters", 
                    features = all.hd.mouse.use, cluster.idents = F) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

all.pct.hd <- all.test$data$pct.exp

head(all.pct.hd)

all.pct.hd <- matrix(all.pct.hd, nrow = 157)
dim(all.pct.hd)
all.pct.hd[1:5,1:5]

colnames(all.pct.hd) <- c("gamma", "gamma*", "sk0", "sk1", "sk10", "sk11", "sk2",
                          "sk3", "sk4", "sk5", "sk6", "sk7", "sk8", "sk9","v0", "v1", "v10", 
                          "v11", "v12", "v13", "v14", "v2", "v3", "v4", "v5", "v6",
                          "v7", "v8", "v9") 

rownames(all.pct.hd) <- all.hd.mouse.use
all.avg.data <- data.frame(gene = all.hd.mouse.use,
                           max.pct = apply(all.pct.hd, 1, max),
                           min.pct = apply(all.pct.hd, 1, min))

head(all.avg.data)
table(all.avg.data$max.pct > 5)

# 76 of the genes are found in more than 5% of the cells in any cluster.

hd.all.use <- all.avg.data %>% filter(max.pct > 5) %>% pull(gene)

ggplot(hd.all.expr, aes(x = n, y = expr)) + geom_point() 

ggplot(hd.all.expr[which(hd.all.expr$gene %in% hd.all.use),], aes(x = n, y = expr)) + geom_point() 

hd.all.git.mat <- hd.all.git.mat[hd.all.use,]
hd.all.git.s <- subset(all.git.MN, features = hd.all.use)
hd.all.git.s
agg.hd.all.git.mat <- AverageExpression(hd.all.git.s, assays = "RNA", group.by = "git_clusters")
agg.hd.all.git.mat <- agg.hd.all.git.mat[[1]]
agg.hd.all.git.mat[1:5,1:5]

colnames(agg.hd.all.git.mat) <- c("gamma", "gamma*", "sk0", "sk1", "sk10", "sk11", "sk2",
                                  "sk3", "sk4", "sk5", "sk6", "sk7", "sk8", "sk9","v0", "v1", "v10", 
                                  "v11", "v12", "v13", "v14", "v2", "v3", "v4", "v5", "v6",
                                  "v7", "v8", "v9")

hd.all.avg.clust <- hclust(dist(t(agg.hd.all.git.mat)), method = "ward.D2")
plot(hd.all.avg.clust)
rect.hclust(hd.all.avg.clust, k = 2)

hd.all.use <- sort(hd.all.use)
gene.order.paschalis <- c(13:34,56:59,47:49,35:38,41:45,39:40,52,1:12, 46,50:51,53:55,60:76)
hd.all.order <- hd.all.use[gene.order.paschalis]
hd.all.order

new.order <- c(3,4,1,2,5:12,16,17,13:15,20:22,18:19,23:76)
hd.all.order <- hd.all.order[new.order]
hd.all.order

ant.post.order <- c("Hoxb2", "Hoxb3", "Hoxd3", "Hoxb4", "Hoxc4", "Hoxb5", "Hoxc5", "Hoxb6", "Hoxa7",
                    "Hoxb7", "Hoxb8", "Hoxd8", "Hoxa9", "Hoxb9", "Hoxd9", "Hoxa10", "Hoxd10", 
                    "Hoxa11", "Hoxc11", "Hoxd11", "Hoxc12", "Hoxc13", hd.all.order[23:76])

ant.post.order

hclust.hd.cell.order <- hd.all.avg.clust$order
hd.cell.order <- colnames(agg.hd.all.git.mat)[hclust.hd.cell.order]

all.git.MN$git_clusters <- factor(all.git.MN$git_clusters, levels = hd.cell.order)

DotPlot(all.git.MN, group.by = "git_clusters", 
        features = rev(hd.all.order), cluster.idents = F, dot.min = 0.05) + coord_flip() +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))



