# Generating figures for adult MN diversity manuscript

setwd("./Dropbox (VU Basic Sciences)/Miller lab/10X Genomics")
library(monocle3)
library(wbData)
library(dplyr)
library(ggplot2)

gids <- wb_load_gene_ids("WS273")

wt.MN <- readRDS("./7595-ST/SoupX/091922_acr2_lin39_combined_wt_VNC_neurons_cds.rds")
wt.GABA <- readRDS("./7595-ST/SoupX/091922_acr2_lin39_combined_GABA_wt_VNC_neurons_cds.rds")
wt.ACh <- readRDS("./7595-ST/SoupX/092022_acr2_lin39_combined_ACh_wt_VNC_neurons_cds.rds")
wt.all <- readRDS("./7595-ST/SoupX/092022_acr2_lin39_all_cells_cds.rds")
MN.all <- readRDS("./7595-ST/SoupX/091922_acr2_lin39_combined_all_genotypes_VNC_neurons_cds.rds")

plot_cells(wt.MN, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

source("./Monocle3Functions.txt")

plot_cells(wt.MN, 
           color_cells_by = "anatomical_class",
           label_cell_groups = F,
           cell_size = 0.25) + guides(color = guide_legend(title="Neuron Class",
                                                           override.aes = list(size = 4)))

colData(wt.MN)$Strain <- dplyr::recode(colData(wt.MN)$Strain, "acr-2" = "acr-2::GFP",
                                       "lin-39" = "lin-39::RFP")

plot_cells(wt.MN, 
           color_cells_by = "Strain",
           label_cell_groups = F,
           cell_size = 0.25) 

colnames(colData(wt.MN))

plot_genes_by_group(wt.MN, markers = c("unc-3", "unc-17", "cha-1", "cho-1", "unc-25", "unc-47", "unc-46"),
                    group_cells_by = "anatomical_class", ordering_type = "none") + labs(x = "Neuron Class", y= "Gene") +
  guides(size = guide_legend(title = "Percent of Cells\nExpressing"))

wt.ACh <- wt.ACh[,colData(wt.ACh)$Cell.type != 'SAB']

plot_cells(wt.ACh, 
           color_cells_by = "anatomical_class",
           label_cell_groups = F,
           cell_size = 0.25) + guides(color = guide_legend(title="Neuron Class",
                                                           override.aes = list(size = 4)))

plot_cells(wt.ACh, 
           color_cells_by = "Cell.type",
           label_cell_groups = F,
           cell_size = 0.5) + guides(color = guide_legend(title="Neuron Class",
                                                          override.aes = list(size = 4)))

plot_cells(wt.ACh,
           color_cells_by = "Cell.type",
           label_groups_by_cluster = F,
           group_label_size = 3,
           cell_size = 0.5)

source("./Monocle3Functions.txt")
plot.expr.UMAP(wt.ACh, "ceh-13", size = 0.5)
plot.expr.UMAP(wt.ACh, "lin-39", size = 0.5)
plot.expr.UMAP(wt.ACh, "mab-5", size = 0.5)
plot.expr.UMAP(wt.ACh, "egl-5", size = 0.5)
plot.expr.UMAP(wt.ACh, "cho-1", size = 0.5)

plot.expr.UMAP(wt.ACh, "mab-5", coexpr_gene = "lin-39", size = 0.5)
plot.expr.UMAP(wt.ACh, "ceh-13", coexpr_gene = "lin-39", size = 0.5)
plot.expr.UMAP(wt.ACh, "ceh-13", coexpr_gene = "mab-5", size = 0.5)

hox <- c("ceh-13", "mab-5", "lin-39", "egl-5")
setwd("..")
for(i in 1:length(hox)){
  gene <- hox[[i]]
  p <- plot.expr.UMAP(wt.ACh, gene, size = 0.25)
  ggsave(p, filename = paste("112222_wt_ACh_VNC_", gene, "_UMAP.pdf", sep = ""), width = 4.5, height = 3, 
         unit ="in")
}

# plotting subsets of each anatomical class
wt.AS <- wt.ACh[,colData(wt.ACh)$anatomical_class == "AS"]
wt.DA <- wt.ACh[,colData(wt.ACh)$anatomical_class == "DA"]
wt.DB <- wt.ACh[,colData(wt.ACh)$anatomical_class == "DB"]
wt.VA <- wt.ACh[,colData(wt.ACh)$anatomical_class == "VA"]
wt.VB <- wt.ACh[,colData(wt.ACh)$anatomical_class == "VB"]
wt.VC <- wt.ACh[,colData(wt.ACh)$anatomical_class == "VC"]

# Re-running dimensionality reduction on each class

dim.reduce <- function(cds) {
  cds <- detect_genes(cds)
  cds <- cds[rowData(cds)$num_cells_expressed > 5,]
  cds <- preprocess_cds(cds, num_dim = 25)
  plot_pc_variance_explained(cds)
  cds <- align_cds(cds, alignment_group = "Strain", alignment_k = 3)
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "Aligned", umap.min_dist = 0.3, umap.n_neighbors = 75)
  cds <- cluster_cells(cds, res = 3e-3)
  colData(cds)$UMAP_1 <- reducedDims(cds)[["UMAP"]][,1]
  colData(cds)$UMAP_2 <- reducedDims(cds)[["UMAP"]][,2]
  cds
}

wt.AS <- dim.reduce(wt.AS)
wt.DA <- dim.reduce(wt.DA)
wt.DB <- dim.reduce(wt.DB)
wt.VA <- dim.reduce(wt.VA)
wt.VB <- dim.reduce(wt.VB)
wt.VC <- dim.reduce(wt.VC)

plot_cells(wt.AS, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(wt.DA, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(wt.DB, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(wt.VA, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(wt.VB, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(wt.VC, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

plot_cells(wt.VC, 
           color_cells_by = "cluster",
           group_label_size = 4,
           label_groups_by_cluster = T,
           cell_size = 0.5)
colData(wt.VC)$cluster <- monocle3::clusters(wt.VC)

VC.clust <- top_markers(wt.VC)

ggplot(as.data.frame(colData(wt.VC)), aes(x = cluster, y = total)) + geom_boxplot()

VC.clust %>% filter(cell_group == 6) %>% arrange(desc(marker_score)) %>% head(10)
# likely AS/VC doublets.

colData(wt.VC)$Cell.type <- ifelse(
  colData(wt.VC)$cluster == 6,
  "Doublets",
  as.character(colData(wt.VC)$Cell.type)
)

colData(wt.ACh)[colnames(wt.VC),]$Cell.type <- colData(wt.VC)$Cell.type
colData(wt.MN)[colnames(wt.VC),]$Cell.type <- colData(wt.VC)$Cell.type
colData(wt.all)[colnames(wt.VC),]$Cell.type <- colData(wt.VC)$Cell.type
colData(MN.all)[colnames(wt.VC),]$Cell.type <- colData(wt.VC)$Cell.type

# Removing doublets

plot_cells(wt.ACh, 
           color_cells_by = "Cell.type",
           label_cell_groups = F,
           cell_size = 0.5) + guides(color = guide_legend(title="Neuron Class",
                                                          override.aes = list(size = 4)))

plot_cells(wt.MN, 
           color_cells_by = "Cell.type",
           label_cell_groups = F,
           cell_size = 0.5) + guides(color = guide_legend(title="Neuron Class",
                                                          override.aes = list(size = 4)))

wt.VC <- wt.VC[,colData(wt.VC)$Cell.type != "Doublets"]
wt.ACh <- wt.ACh[,colData(wt.VC)$Cell.type != "Doublets"]
wt.MN <- wt.MN[,colData(wt.VC)$Cell.type != "Doublets"]

table(colData(MN.all)$Genotype, colData(MN.all)$Cell.type)

# Without doublets
plot_cells(wt.MN, 
           color_cells_by = "anatomical_class",
           label_cell_groups = F,
           cell_size = 0.5) + guides(color = guide_legend(title="Neuron Class",
                                                          override.aes = list(size = 4)))

colData(wt.MN)$Strain <- dplyr::recode(colData(wt.MN)$Strain, "acr-2" = "acr-2::GFP",
                                       "lin-39" = "lin-39::RFP")

plot_cells(wt.MN, 
           color_cells_by = "Strain",
           label_cell_groups = F,
           cell_size = 0.5) 

colnames(colData(wt.MN))

plot_genes_by_group(wt.MN, markers = c("unc-3", "unc-17", "cha-1", "cho-1", "unc-25", "unc-47", "unc-46"),
                    group_cells_by = "anatomical_class", ordering_type = "none") + labs(x = "Neuron Class", y= "Gene") +
  guides(size = guide_legend(title = "Percent of Cells\nExpressing"))

plot_genes_by_group(wt.MN, markers = c("ceh-12", "bnc-1", "unc-4", "unc-129", "vab-7",  "unc-55", "ttr-39", "unc-25"),
                    group_cells_by = "anatomical_class", ordering_type = "none") + labs(x = "Neuron Class", y= "Gene") +
  guides(size = guide_legend(title = "Percent of Cells\nExpressing"))

wt.ACh <- wt.ACh[,colData(wt.ACh)$Cell.type != 'SAB']

plot_cells(wt.ACh, 
           color_cells_by = "anatomical_class",
           label_cell_groups = F,
           cell_size = 0.5) + guides(color = guide_legend(title="Neuron Class",
                                                          override.aes = list(size = 4)))

plot_cells(wt.ACh, 
           color_cells_by = "Cell.type",
           label_cell_groups = F,
           cell_size = 0.5) + guides(color = guide_legend(title="Neuron Class",
                                                          override.aes = list(size = 4)))

source("./Monocle3Functions.txt")
plot.expr.UMAP(wt.ACh, "ceh-13", size = 0.5)
plot.expr.UMAP(wt.ACh, "lin-39", size = 0.5)
plot.expr.UMAP(wt.ACh, "mab-5", size = 0.5)
plot.expr.UMAP(wt.ACh, "egl-5", size = 0.5)

plot.expr.UMAP(wt.ACh, "mab-5", coexpr_gene = "lin-39", size = 0.5)
plot.expr.UMAP(wt.ACh, "ceh-13", coexpr_gene = "lin-39", size = 0.5)
plot.expr.UMAP(wt.ACh, "ceh-13", coexpr_gene = "mab-5", size = 0.5)

# plotting subsets of each anatomical class
wt.AS <- wt.ACh[,colData(wt.ACh)$anatomical_class == "AS"]
wt.DA <- wt.ACh[,colData(wt.ACh)$anatomical_class == "DA"]
wt.DB <- wt.ACh[,colData(wt.ACh)$anatomical_class == "DB"]
wt.VA <- wt.ACh[,colData(wt.ACh)$anatomical_class == "VA"]
wt.VB <- wt.ACh[,colData(wt.ACh)$anatomical_class == "VB"]
wt.VC <- wt.ACh[,colData(wt.ACh)$anatomical_class == "VC"]

# Re-running dimensionality reduction on each class

dim.reduce <- function(cds) {
  cds <- detect_genes(cds)
  cds <- cds[rowData(cds)$num_cells_expressed > 5,]
  cds <- preprocess_cds(cds, num_dim = 25)
  plot_pc_variance_explained(cds)
  cds <- align_cds(cds, alignment_group = "Strain", alignment_k = 3)
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "Aligned", umap.min_dist = 0.3, umap.n_neighbors = 75)
  cds <- cluster_cells(cds, res = 3e-3)
  colData(cds)$UMAP_1 <- reducedDims(cds)[["UMAP"]][,1]
  colData(cds)$UMAP_2 <- reducedDims(cds)[["UMAP"]][,2]
  cds
}

wt.VC <- dim.reduce(wt.VC)

plot_cells(wt.VC, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

# GABA

wt.GABA <- wt.GABA[,colData(wt.GABA)$Cell.type != 'VD_unc-3?']

plot_cells(wt.GABA, 
           color_cells_by = "anatomical_class",
           label_cell_groups = F,
           cell_size = 0.5) + guides(color = guide_legend(title="Neuron Class",
                                                          override.aes = list(size = 4)))

# wt.VA plots
genes <- c("nlp-7", "flp-11", "sfrp-1", "nlp-11", "vab-3", "avr-15", "flp-18", "sra-36", "mig-13",
           "egl-44", "nlp-13")

for(i in 1:length(genes)) {
  gene <- genes[[i]]
  p <- plot.expr.UMAP(wt.VA, gene, size = 0.25)
  ggsave(plot = p, filename = paste("~/Dropbox (VU Basic Sciences)/sc-RNA-Seq MN paper/FIGURES/Seth_Plots/UMAP_and_expression_plots/Jan_23_pdfs/042423_VA_only_wt_", gene, "_UMAP.pdf", sep = ""),
         width = 4.5, height = 3, unit = "in")
}

# Plotting among all cholinergic MNs
genes <- c("nlp-7", "flp-11", "sfrp-1", "nlp-11", "vab-3", "avr-15", "flp-18", "sra-36", "mig-13",
           "egl-44")

for(i in 1:length(genes)) {
  gene <- genes[[i]]
  p <- plot.expr.UMAP(wt.ACh, gene, size = 0.25)
  ggsave(plot = p, filename = paste("~/Dropbox (VU Basic Sciences)/sc-RNA-Seq MN paper/FIGURES/Seth_Plots/UMAP_and_expression_plots/Jan_23_pdfs/030823_ACh_only_wt_", gene, "_UMAP.pdf", sep = ""),
         width = 4.5, height = 3, unit = "in")
}

# Running differential expression between subtypes within each anatomical class

table(colData(wt.MN)$Cell.type, colData(wt.MN)$anatomical_class)

colData(wt.MN)$Cell.type <- dplyr::recode(colData(wt.MN)$Cell.type, "VD_10_11?" = "VD_8_11")

colData(wt.MN)$anatomical_class <- ifelse(
  colData(wt.MN)$Cell.type %in% c("VD_8_11", "VD"),
  "VD",
  ifelse(
    colData(wt.MN)$Cell.type %in% c("DD", "DD_4_5"),
    "DD",
    as.character(colData(wt.MN)$anatomical_class)
  )
)
MN.wt.S <- CreateSeuratObject(counts = exprs(wt.MN),
                              meta.data = as.data.frame(colData(wt.MN)))

MN.wt.S <- NormalizeData(MN.wt.S)

# Now to split and run FindMarkers
Idents(MN.wt.S) <- MN.wt.S@meta.data$anatomical_class
classes <- c("AS", "DA", "DB", "DD", "VA", "VB", "VC", "VD")
subclass.DE <- list()
for(i in 1:length(classes)){
  class <- classes[[i]]
  sub.S <- subset(MN.wt.S, subset = anatomical_class == class)
  Idents(sub.S) <- sub.S@meta.data$Cell.type
  sub.markers <- FindAllMarkers(sub.S, only.pos = T)
  sub.markers$gene_name <- i2s(sub.markers$gene, gids)
  sub.markers$class <- class
  subclass.DE[[i]] <- sub.markers
}

sub.markers.df <- do.call("rbind", subclass.DE)
head(sub.markers.df)

sub.markers.df %>% arrange(desc(avg_log2FC)) %>% head(10)

sub.markers.df <- sub.markers.df[,c(9,6,7,8,2,3,4,1,5)]
sub.markers.df %>% arrange(desc(avg_log2FC)) %>% head(10)
colnames(sub.markers.df) <- c("anatomical_class", "subclass", "WBGeneID", "gene_symbol", "avg_log2FC",
                              "pct_of_subclass", "pct_of_other_subclasses", "p-value", "adj_p-value")

sub.markers.df <- sub.markers.df %>% arrange(anatomical_class, subclass, `adj_p-value`, desc(avg_log2FC))
head(sub.markers.df)

sub.markers.df %>% filter(subclass == "AS_2_3") %>% head(10)
sub.markers.df %>% filter(subclass == "DD_4_5") %>% head(10)
sub.markers.df %>% filter(subclass == "VD12") %>% head(10)

table(sub.markers.df$`adj_p-value` < 0.05, sub.markers.df$subclass)

sub.markers.df %>% filter(subclass == "VA10") %>% head(10)

saveRDS(sub.markers.df, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/7595-ST/SoupX/012323_acr2_lin39_wt_VNC_MN_within_class_subclass_markers.rds")
write.table(sub.markers.df, "./UMAP_and_expression_plots/Jan_23_pdfs/012323_all_wt_VNC_MN_within_class_subclass_DE.csv",
            sep = ",", quote = F, row.names = F)

