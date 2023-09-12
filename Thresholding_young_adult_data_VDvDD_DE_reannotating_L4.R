# Thresholding the wt data from the acr-2/lin-39 combined data

library(monocle3)
library(ggplot2)
library(wbData)
gids <- wb_load_gene_ids("WS273")

setwd("./Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/")

YA.all <- readRDS("./Kratsios_scRNAseq/7595-ST/SoupX/092022_acr2_lin39_all_cells_cds.rds")
table(colData(YA.all)$Genotype, colData(YA.all)$Cell.type)

YA.wt <- YA.all[,colData(YA.all)$Genotype == "wt"]
table(colData(YA.wt)$Experiment, colData(YA.wt)$Cell.type)
YA.wt <- YA.wt[,colData(YA.wt)$Experiment != "lin-39_YA_wt_3"]

table(colData(YA.wt)$Tissue)

YA.wt <- detect_genes(YA.wt)
YA.wt <- YA.wt[rowData(YA.wt)$num_cells_expressed > 5,]
YA.wt

# 15049 features in 29378 cells

head(colData(YA.wt))

YA.wt <- preprocess_cds(YA.wt, method = "PCA", num_dim = 50)
# Removing doublets and low-quality stressed cells.
YA.wt <- YA.wt[,!colData(YA.wt)$Cell.type %in% c("Doublets", "Doublets_lowq", "Doublets_2")]
YA.wt <- YA.wt[,!colData(YA.wt)$Cell.type %in% c("VA.VB.UPR", "VB01.UPR", "VB02.UPR")]
YA.wt <- detect_genes(YA.wt)
YA.wt <- YA.wt[rowData(YA.wt)$num_cells_expressed > 5,]
YA.wt
# 14993 features in 26524 cells

YA.wt <- preprocess_cds(YA.wt, method = "PCA", num_dim = 50)
plot_pc_variance_explained(YA.wt)

YA.wt <- reduce_dimension(YA.wt, 
                          reduction_method = "UMAP",
                          umap.min_dist = 0.3,
                          umap.n_neighbors = 75)

YA.wt <- cluster_cells(YA.wt, res = 3e-3)

plot_cells(YA.wt, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

colData(YA.wt)$UMAP_1 <- reducedDims(YA.wt)[["UMAP"]][,1]
colData(YA.wt)$UMAP_2 <- reducedDims(YA.wt)[["UMAP"]][,2]

colData(YA.wt)$cluster <- monocle3::clusters(YA.wt)

table(colData(YA.wt)$Cell.type)

plot_cells(YA.wt, 
           color_cells_by = "Tissue",
           label_cell_groups = F,
           cell_size = 0.5)

source("./Monocle3Functions.txt")

plot.expr.UMAP(YA.wt, "sbt-1", size = 0.5)

table(colData(YA.wt)$Tissue, exclude = NULL)

plot_cells(YA.wt, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

saveRDS(YA.wt, "./Kratsios_scRNAseq/030823_YA_wt_neuron_cds.rds")

plot.expr.UMAP(YA.wt, "egl-21", size = 0.5)
plot.expr.UMAP(YA.wt, "sbt-1", size = 0.5)
plot.cell.type.m3(YA.wt, "CAN")
library(dplyr)

# Subsetting out just the neuronal clusters
YA.neuron <- YA.wt[,colData(YA.wt)$cluster %in% c(25,78,77,70,45,48,12,44,1,40,17,29,14,26,20,
                                                  53,52,60,75,7,71,8,22,74,47,72,2,6,67,27,13,
                                                  61,73,68,10,16,49,50,30,46,39,9,24,4,11,69,
                                                  76,84,3,63,28,41) | colData(YA.wt)$Cell.type == "CAN"]

colData(YA.neuron)$Tissue <- ifelse(
  colData(YA.neuron)$Cell.type %in% c("Spermathecal-uterine_junction_or_uterine_toroid",
                                      "Germline", "Intestine", "Pharyngeal_gland_cell", "Epidermis", 
                                      "Uterine_cell"),
  "Nonneuronal",
  "Neuron"
)

colData(YA.wt)[colnames(YA.neuron),]$Tissue <- colData(YA.neuron)$Tissue
colData(YA.wt)$Tissue <- ifelse(
  colData(YA.wt)$Tissue == "Neuron",
  "Neuron",
  "Nonneuronal"
)

YA.neuron <- YA.neuron[,!colData(YA.neuron)$Cell.type %in% 
                         c("Spermathecal-uterine_junction_or_uterine_toroid", "Unknown_unc3_MN_3",
                           "Germline", "Intestine", "Pharyngeal_gland_cell", "Epidermis", 
                           "Uterine_cell", "Unknown_unc3_MN_1", "VD_unc-3?")]

YA.neuron <- detect_genes(YA.neuron)
YA.neuron <- YA.neuron[rowData(YA.neuron)$num_cells_expressed > 5,]
YA.neuron
# 12528 features in 17420 cells

YA.neuron <- preprocess_cds(YA.neuron, method = "PCA", num_dim = 50)
plot_pc_variance_explained(YA.neuron)

YA.neuron <- reduce_dimension(YA.neuron, 
                              reduction_method = "UMAP",
                              umap.min_dist = 0.3,
                              umap.n_neighbors = 75)

YA.neuron <- cluster_cells(YA.neuron, res = 3e-3)

plot_cells(YA.neuron, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

colData(YA.neuron)$UMAP_1 <- reducedDims(YA.neuron)[["UMAP"]][,1]
colData(YA.neuron)$UMAP_2 <- reducedDims(YA.neuron)[["UMAP"]][,2]

colData(YA.neuron)$cluster <- monocle3::clusters(YA.neuron)

# Standardizing names
colData(YA.neuron)$Cell.type <- dplyr::recode(colData(YA.neuron)$Cell.type,
                                              "VD10_11?" = "VD_8_11",
                                              "VD_10_11?" = "VD_8_11",
                                              "DD4_5" = "DD_4_5",
                                              "AS02" = "AS_2_3",
                                              "VD_DD" = "VD")

table(colData(YA.neuron)$Cell.type)

plot_cells(YA.neuron, 
           color_cells_by = "Tissue",
           label_cell_groups = F,
           cell_size = 0.5)

source("./Monocle3Functions.txt")

plot.expr.UMAP(YA.neuron, "sbt-1", size = 0.5)

table(colData(YA.neuron)$Tissue, exclude = NULL)

plot_cells(YA.neuron, 
           color_cells_by = "cluster",
           group_label_size = 4,
           cell_size = 0.5)

cell.table <- as.data.frame(table(colData(YA.neuron)$Cell.type))

colnames(cell.table) <- c("Cell.type", "n")
head(cell.table)
keep <- cell.table %>% filter(n > 20) %>% select(Cell.type)
keep <- as.character(keep$Cell.type)
keep

saveRDS(YA.neuron, "./Kratsios_scRNAseq/7595-ST/SoupX/060823_all_neurons_for_thresholding_cds.rds")

YA.truth <- read.table("./Kratsios_scRNAseq/101822_ground_truth_expression_for_YA_wt.csv", sep =",",
                       header = T, row.names = 1, stringsAsFactors = F)

rownames(YA.truth) <- s2i(rownames(YA.truth), gids)

YA_hbox_truth <- YA.truth[1:111,]
YA_metabo <- YA.truth[112:119,]
YA_inx <- YA.truth[120:133,]
YA_pan.non <- YA.truth[134:168,]

YA_datasets <- list(YA_hbox_truth, YA_metabo, YA_inx, YA_pan.non)

#~ Read datasets ----

keep <- keep[c(1:31,33:49)]
keep
YA.expr <- get.expr.info.by.facet.3(YA.neuron[,colData(YA.neuron)$Cell.type %in% keep], "Cell.type")
YA.new.prop <- YA.expr[[2]]
YA.new.TPM <- YA.expr[[1]]

table(colData(YA.neuron[,colData(YA.neuron)$Cell.type %in% keep])$Cell.type)

# Proportion of cells with at least 1 UMI (generated separately)
YA_prop_by_type <- YA.new.prop
dim(YA_prop_by_type)
colnames(YA_prop_by_type)
YA_prop_by_type <- YA_prop_by_type[, c(1:48)]
dim(YA_prop_by_type)
YA.new.TPM <- YA.new.TPM[,c(1:48)]
colnames(YA.new.TPM)

YA_ground_truth <- bind_rows(YA_datasets)

head(rownames(YA_ground_truth))
head(rownames(YA_prop_by_type))
library(tidyverse)
YA_ground_truth$type <- rep(c("hbox", "metabo", "inx", "pan"), times = map_int(YA_datasets, nrow))
YA_ground_truth <- YA_ground_truth[rownames(YA_ground_truth) %in% rownames(YA_prop_by_type), ] # no expression data

# save type stratifying column with same reshaping
YA_source_types <- YA_ground_truth$type
YA_ground_truth$type <- NULL

all.equal(dim(YA_ground_truth), c(136, 48))

#~ Definitions ----

# hard thresholds fixed, different from YA (low_hard = 0.02, high_hard = 0.01)
low_hard <- 0.02
high_hard <- 0.01

get_thres <- function(vals, level = 0.04, hard_low = 0.02, hard_high = 0.01){
  if(sum(vals) == 0) return(1)
  if(max(vals) < hard_low) return(1)
  if(min(vals) > hard_high) return(0)
  return(level*max(vals))
}


#~~~~~~~~~~~~~~~~~~~~~~~ ----


# Bootstraps ----

# Idea: for each threshold level, compute the bin matrix and run bootstraps

binarize <- function(dyn_lev, props, truth){
  cbind(1L*(props >= apply(props, 1, get_thres, dyn_lev, low_hard, high_hard)),
        truth)
}

#~~ define threshold percentiles ----


levels <- c(0, seq(1e-3,1e-2, 1e-4),
            seq(1.1e-2, 1e-1, 1e-3),
            seq(1.1e-1,1, 1e-2))

length(intersect(rownames(YA_ground_truth), rownames(YA_prop_by_type)))
#138
length(intersect(colnames(YA_ground_truth), colnames(YA_prop_by_type)))
#48

setdiff(colnames(YA_ground_truth), colnames(YA_prop_by_type))
setdiff(colnames(YA_prop_by_type), colnames(YA_ground_truth))

colnames(YA_ground_truth)

YA_all_bins <- map(levels,
                   binarize,
                   YA_prop_by_type[rownames(YA_ground_truth), colnames(YA_ground_truth)],
                   YA_ground_truth)
names(YA_all_bins) <- as.character(levels)

#~~ callback functions ----
boot_tpr <- function(df_boot, i){
  # boot provides the original df and a vector of indexes selected in that bootstrap
  boot_bin <- df_boot[i,1:48]
  boot_truth <- df_boot[i,49:96]
  TPR <- sum(boot_bin * boot_truth)/sum(boot_truth)
}
boot_fdr <- function(df_boot, i){
  # boot provides the original df and a vector of indexes selected in that bootstrap
  boot_bin <- df_boot[i,1:48]
  boot_truth <- df_boot[i,49:96]
  FDR <- sum(boot_bin * (!boot_truth))/(sum(boot_bin))
}
boot_fpr <- function(df_boot, i){
  # boot provides the original df and a vector of indexes selected in that bootstrap
  boot_bin <- df_boot[i,1:48]
  boot_truth <- df_boot[i,49:96]
  FPR <- sum(boot_bin * (!boot_truth))/sum(!(boot_truth))
}


#~~ Loop bootstraps ----
## Running on Windows with snow -Running on mac with "snow", 4 cpus
YA_res_tpr <- map_dfr(YA_all_bins, ~{
  res_boot_tpr <- boot(., boot_tpr, R=5000, strata = as.integer(as.factor(YA_source_types)),
                       parallel="snow", ncpus=4)
  set_names(c(res_boot_tpr$t0,boot.ci(res_boot_tpr, type = "bca")$bca[4:5]),
            c("t0", "lower", "upper"))
})
YA_res_fdr <- map_dfr(YA_all_bins, ~{
  res_boot_tpr <- boot(., boot_fdr, R=5000, strata = as.integer(as.factor(YA_source_types)),
                       parallel="snow", ncpus=4)
  set_names(c(res_boot_tpr$t0,boot.ci(res_boot_tpr, type = "bca")$bca[4:5]),
            c("t0", "lower", "upper"))
})
YA_res_fpr <- map_dfr(YA_all_bins, ~{
  res_boot_tpr <- boot(., boot_fpr, R=5000, strata = as.integer(as.factor(YA_source_types)),
                       parallel="snow", ncpus=4)
  set_names(c(res_boot_tpr$t0,boot.ci(res_boot_tpr, type = "perc")$perc[4:5]),
            c("t0", "lower", "upper"))
})

names(YA_res_tpr) <- c("TPR", "TPR_lower", "TPR_upper")
names(YA_res_fpr) <- c("FPR", "FPR_lower", "FPR_upper")
names(YA_res_fdr) <- c("FDR", "FDR_lower", "FDR_upper")

YA_res_boot <- bind_cols(percentile = names(YA_all_bins), YA_res_tpr, YA_res_fpr, YA_res_fdr)

#~~~~~~~~~~~~~~~~~----

# Plots results ----

# -> Figures S6C, D

YA_chosen_levels <- c(0.02, 0.06, 0.13, 0.18)

ggplot(YA_res_boot, aes(x=1-FDR, y = TPR)) +
  theme_classic() +
  geom_ribbon(aes(ymin=TPR_lower, ymax=TPR_upper), fill = "gray") +
  geom_ribbon(aes(xmin=1-FDR_lower, xmax=1-FDR_upper), fill = "gray") +
  # geom_errorbar(aes(ymin=TPR_lower, ymax=TPR_upper), color = "gray") +
  # geom_errorbarh(aes(xmin=1-FDR_lower, xmax=1-FDR_upper), color = "gray") +
  geom_point() +
  geom_line() +
  # scale_x_continuous(limits=c(0.5,1)) +
  # scale_y_continuous(limits=c(0,1)) +
  xlab("Precision") +
  ylab("Recall") +
  geom_point(data=filter(YA_res_boot, percentile %in% YA_chosen_levels),
             color = "red", size=2)

# ggsave("output/PR_rib.png",
#        width = 5, height = 4, unit = "in")

ggplot(YA_res_boot, aes(x=FPR, y = TPR)) +
  theme_classic() +
  geom_ribbon(aes(ymin=TPR_lower, ymax=TPR_upper), fill = "gray") +
  geom_ribbon(aes(xmin=FPR_lower, xmax=FPR_upper), fill = "gray") +
  # geom_errorbar(aes(ymin=TPR_lower, ymax=TPR_upper), color = "gray") +
  # geom_errorbarh(aes(xmin=FPR_lower, xmax=FPR_upper), color = "gray") +
  geom_point() +
  geom_line() +
  # scale_x_continuous(limits=c(0,.5)) +
  # scale_y_continuous(limits=c(0,1)) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  geom_point(data=filter(YA_res_boot, percentile %in% YA_chosen_levels),
             color = "red", size=2)

# ggsave("output/ROC_rib.png",
#        width = 5, height = 4, unit = "in")

head(YA_res_boot)

# Statistics at certain levels

YA_res_boot %>% filter(percentile %in% c(0.02, 0.06, 0.13, 0.18))
# A tibble: 4 × 10
# percentile   TPR TPR_lower TPR_upper    FPR FPR_lower FPR_upper   FDR FDR_lower FDR_upper
# <chr>      <dbl>     <dbl>     <dbl>  <dbl>     <dbl>     <dbl> <dbl>     <dbl>     <dbl>
# 1 0.02       0.899     0.842     0.929 0.153     0.111     0.198  0.183    0.135      0.249
# 2 0.06       0.864     0.808     0.897 0.107     0.0714    0.147  0.141    0.0959     0.207
# 3 0.13       0.782     0.726     0.824 0.0735    0.0474    0.104  0.110    0.0719     0.173
# 4 0.18       0.726     0.665     0.771 0.0621    0.0384    0.0897 0.102    0.0654     0.167

YA_res_boot %>% filter(percentile %in% c(0.02, 0.04, 0.09, 0.15))

# A tibble: 4 × 10
# percentile   TPR TPR_lower TPR_upper    FPR FPR_lower FPR_upper   FDR FDR_lower FDR_upper
# <chr>      <dbl>     <dbl>     <dbl>  <dbl>     <dbl>     <dbl> <dbl>     <dbl>     <dbl>
# 1 0.02       0.899     0.842     0.929 0.153     0.111     0.198  0.183    0.135      0.249
# 2 0.04       0.883     0.830     0.915 0.123     0.0853    0.166  0.155    0.110      0.222
# 3 0.09       0.832     0.779     0.868 0.0912    0.0589    0.130  0.127    0.0832     0.193
# 4 0.15       0.761     0.706     0.805 0.0671    0.0429    0.0962 0.104    0.0687     0.165

YA_res_boot %>% filter(percentile %in% c(0.04, 0.08, 0.14, 0.2))

chosen_levels <- c(0.015, 0.035, 0.09, 0.15)

#~~ Save csv ----

saveRDS(YA_res_boot, "./Kratsios_scRNAseq/Analysis/060823_thresholding_results.rds")

# Using these settings to generate binary expression matrices
YA.TPM.use <- YA.new.TPM[,colnames(YA_ground_truth)]
YA.prop.use <- YA.new.prop[,colnames(YA_ground_truth)]

dim(YA.prop.use)
dim(YA.TPM.use)

make_expr_matrix <- function(values, dynamic_level, hard_low, hard_high){
  exprs <- matrix(nrow = nrow(YA.prop.use),
                  ncol = ncol(YA.prop.use),
                  dimnames = dimnames(YA.prop.use))
  for(g in rownames(exprs)){
    exprs[g,] <- values[g,] >= get_thres(values[g,],
                                         dynamic_level,
                                         hard_low,
                                         hard_high)
  }
  return(1L*exprs)
}

# make the 3 matrices (cf Word doc for choice of thresholds)
# can take a few dozens of seconds
YA_exprs_liberal <- make_expr_matrix(YA.prop.use,
                                     dynamic_level = 0.02,
                                     hard_low = 0.02,
                                     hard_high = 0.01)

YA_exprs_medium <- make_expr_matrix(YA.prop.use,
                                    dynamic_level = 0.06,
                                    hard_low = 0.02,
                                    hard_high = 0.01)

YA_exprs_conservative <- make_expr_matrix(YA.prop.use,
                                          dynamic_level = 0.13,
                                          hard_low = 0.02,
                                          hard_high = 0.01)

YA_exprs_stringent <- make_expr_matrix(YA.prop.use, 
                                       dynamic_level = 0.18, 
                                       hard_low = 0.02,
                                       hard_high = 0.01)

saveRDS(YA.prop.use, "./Kratsios_scRNAseq/Analysis/060823_YA_proportion_matrix.rds")
saveRDS(YA_exprs_liberal, "./Kratsios_scRNAseq/Analysis/060823_YA_t1_binary_expression_matrix.rds")
saveRDS(YA_exprs_medium, "./Kratsios_scRNAseq/Analysis/060823_YA_t2_binary_expression_matrix.rds")
saveRDS(YA_exprs_conservative, "./Kratsios_scRNAseq/Analysis/060823_YA_t3_binary_expression_matrix.rds")
saveRDS(YA_exprs_stringent, "./Kratsios_scRNAseq/Analysis/060823_YA_t4_binary_expression_matrix.rds")
YA_exprs_liberal[1:5,1:5]

# Applying a thresholded logical matrix to the TPM matrices to set genes below the threshold to 0 and retain
# the continuous values for genes above the threshold.

YA_exprs_liberal<- as.data.frame(YA_exprs_liberal)
YA_exprs_liberal$gene <- rownames(YA_exprs_liberal)
colnames(YA_exprs_liberal)
cell.types <- colnames(YA_exprs_liberal)[1:48]

YA_exprs_medium<- as.data.frame(YA_exprs_medium)
YA_exprs_medium$gene <- rownames(YA_exprs_medium)
YA_exprs_conservative<- as.data.frame(YA_exprs_conservative)
YA_exprs_conservative$gene <- rownames(YA_exprs_conservative)

YA_exprs_stringent <- as.data.frame(YA_exprs_stringent)
YA_exprs_stringent$gene <- rownames(YA_exprs_stringent)

length(rownames(YA_exprs_liberal))
length(rownames(YA.TPM.use))

identical(rownames(YA_exprs_liberal), rownames(YA.TPM.use))

create.prp.04.cells <- function(Neuron.type, data) {
  tmp <- data[which(data[, Neuron.type] > 0), c(Neuron.type, "gene")]
  keep <- tmp[,2]
  tmp.expr <- data.frame(gene = keep, value = YA.TPM.use[keep, Neuron.type])
  tmp.expr
}
colnames(YA.TPM.use)
colnames(YA_exprs_liberal)

setdiff(cell.types, colnames(YA.TPM.use))

YA.TPM.liberal <- lapply(cell.types, create.prp.04.cells, YA_exprs_liberal)
YA.TPM.medium <- lapply(cell.types, create.prp.04.cells, YA_exprs_medium)
YA.TPM.conservative <- lapply(cell.types, create.prp.04.cells, YA_exprs_conservative)
YA.TPM.stringent <- lapply(cell.types, create.prp.04.cells, YA_exprs_stringent)

dim(YA.TPM.liberal[[1]])
dim(YA.TPM.medium[[1]])
dim(YA.TPM.conservative[[1]])

dim(YA.TPM.liberal[[2]])
dim(YA.TPM.medium[[2]])
dim(YA.TPM.conservative[[2]])

# Combining the lists
YA.TPM.liberal.mat <- YA.TPM.liberal %>% Reduce(function(df1,df2) full_join(df1,df2,by = "gene"), .)
YA.TPM.liberal.mat[1:5,1:5]
dim(YA.TPM.liberal.mat)
YA.TPM.liberal.mat[is.na(YA.TPM.liberal.mat)] <- 0
rownames(YA.TPM.liberal.mat) <- YA.TPM.liberal.mat$gene
dim(YA.TPM.liberal.mat)
YA.TPM.liberal.mat <- as.matrix(YA.TPM.liberal.mat[,2:49])
dim(YA.TPM.liberal.mat)
YA.TPM.liberal.mat[1:5,1:5]
colnames(YA.TPM.liberal.mat) <- cell.types

YA.TPM.liberal.mat[1:5,1:5]
YA.TPM.liberal.mat[c("WBGene00022277", "WBGene00022276", "WBGene00022278"), 1:5]
dim(YA.TPM.liberal.mat)

# Medium
YA.TPM.medium.mat <- YA.TPM.medium %>% Reduce(function(df1,df2) full_join(df1,df2,by = "gene"), .)
YA.TPM.medium.mat[1:5,1:5]
YA.TPM.medium.mat[is.na(YA.TPM.medium.mat)] <- 0
rownames(YA.TPM.medium.mat) <- YA.TPM.medium.mat$gene
YA.TPM.medium.mat <- as.matrix(YA.TPM.medium.mat[,2:49])

YA.TPM.medium.mat[1:5,1:5]
colnames(YA.TPM.medium.mat) <- cell.types

YA.TPM.medium.mat[1:5,1:5]

# Conservative
YA.TPM.conservative.mat <- YA.TPM.conservative %>% Reduce(function(df1,df2) full_join(df1,df2,by = "gene"), .)
YA.TPM.conservative.mat[1:5,1:5]
YA.TPM.conservative.mat[is.na(YA.TPM.conservative.mat)] <- 0
rownames(YA.TPM.conservative.mat) <- YA.TPM.conservative.mat$gene
YA.TPM.conservative.mat <- as.matrix(YA.TPM.conservative.mat[,2:49])

YA.TPM.conservative.mat[1:5,1:5]
colnames(YA.TPM.conservative.mat) <- cell.types

YA.TPM.conservative.mat[1:5,1:5]

# Stringent
YA.TPM.stringent.mat <- YA.TPM.stringent %>% Reduce(function(df1,df2) full_join(df1,df2,by = "gene"), .)
YA.TPM.stringent.mat[1:5,1:5]
YA.TPM.stringent.mat[is.na(YA.TPM.stringent.mat)] <- 0
rownames(YA.TPM.stringent.mat) <- YA.TPM.stringent.mat$gene
YA.TPM.stringent.mat <- as.matrix(YA.TPM.stringent.mat[,2:49])

YA.TPM.stringent.mat[1:5,1:5]
colnames(YA.TPM.stringent.mat) <- cell.types

YA.TPM.stringent.mat[1:5,1:5]

YA.TPM.raw.mat <- as(as.matrix(YA.TPM.use), "sparseMatrix")
YA.TPM.liberal.mat <- as(as.matrix(YA.TPM.liberal.mat), "sparseMatrix")
YA.TPM.medium.mat <- as(as.matrix(YA.TPM.medium.mat), "sparseMatrix")
YA.TPM.conservative.mat <- as(as.matrix(YA.TPM.conservative.mat), "sparseMatrix")
YA.TPM.stringent.mat <- as(as.matrix(YA.TPM.stringent.mat), "sparseMatrix")

YA.unthresholded.genes <- data.frame(cell.type = colnames(YA.TPM.raw.mat), ngenes = diff(YA.TPM.raw.mat@p))
YA.liberal.genes <- data.frame(cell.type = colnames(YA.TPM.liberal.mat), ngenes = diff(YA.TPM.liberal.mat@p))
YA.medium.genes <- data.frame(cell.type = colnames(YA.TPM.medium.mat), ngenes = diff(YA.TPM.medium.mat@p))
YA.conservative.genes <- data.frame(cell.type = colnames(YA.TPM.conservative.mat), ngenes = diff(YA.TPM.conservative.mat@p))
YA.stringent.genes <- data.frame(cell.type = colnames(YA.TPM.stringent.mat), ngenes = diff(YA.TPM.stringent.mat@p))

YA.TPM.raw.mat[c("WBGene00022277", "WBGene00022276", "WBGene00022278"), 1:5]
YA.TPM.liberal.mat[c("WBGene00022277", "WBGene00022276", "WBGene00022278"), 1:5]
YA.TPM.medium.mat[c("WBGene00022277", "WBGene00022276", "WBGene00022278"), 1:5]
YA.TPM.conservative.mat[c("WBGene00022277", "WBGene00022276", "WBGene00022278"), 1:5]
YA.TPM.stringent.mat[c("WBGene00022277", "WBGene00022276", "WBGene00022278"), 1:5]

saveRDS(YA.TPM.liberal.mat, "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t1_expression_matrix.rds")
saveRDS(YA.TPM.medium.mat, "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t2_expression_matrix.rds")
saveRDS(YA.TPM.conservative.mat, "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t3_expression_matrix.rds")
saveRDS(YA.TPM.stringent.mat, "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t4_expression_matrix.rds")
saveRDS(YA.TPM.use, "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_raw_expression_matrix.rds")
YA.TPM.liberal.mat[1:4,1:5]
colnames(YA.prop.use)

dim(YA.TPM.use)

dim(YA.TPM.liberal.mat)
VNC.sub <- sort(as.character(c(unique(colData(ACh.wt)$Cell.type), unique(colData(VD.DD.wt)$Cell.type))))
VNC.sub
length(intersect(VNC.sub, keep))
setdiff(keep, VNC.sub)
VNC.sub[3] <- "AS_2_3"
VNC.sub <- sort(VNC.sub)
setdiff(keep, VNC.sub)

write.table(as.matrix(YA.TPM.liberal.mat), "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t1_expression.csv",
            sep = ",", quote = F, row.names = T, col.names = NA)
write.table(as.matrix(YA.TPM.medium.mat), "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t2_expression.csv",
            sep = ",", quote = F, row.names = T, col.names = NA)
write.table(as.matrix(YA.TPM.conservative.mat), "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t3_expression.csv",
            sep = ",", quote = F, row.names = T, col.names = NA)
write.table(as.matrix(YA.TPM.stringent.mat), "./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t4_expression.csv",
            sep = ",", quote = F, row.names = T, col.names = NA)

# Chosen levels = 0.02, 0.06, 0.13, 0.18, hard_low = 0.02, hard_high = 0.01

# Generating long format TPM and prop data for making dotplots

setwd("./Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/")
library(dplyr)
library(ggplot2)
library(monocle3)
library(Seurat)

YA.2.tpm <- readRDS("./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t2_expression_matrix.rds")
YA_prop_by_type <- readRDS("./Kratsios_scRNAseq/Analysis/060823_YA_proportion_matrix.rds")
YA.neuron <- readRDS("./Kratsios_scRNAseq/7595-ST/SoupX/060823_all_neurons_for_thresholding_cds.rds")

# now I need to make these into the long format (one row for each gene in each neuron) and 
# add in the proportion data. Then I can generate the code for this

# For neuropeptides, I should use a scaled version of the data.
neuron.order <- c("VA1", "VA2", "VA3_8", "VA9_10", "VA11", "VA12", 
                  "VB1", "VB2", "VB3", "VB4_9", "VB10_11",
                  "DA1", "DA2_5", "DA6_8", "DA9",
                  "DB1", "DB2", "DB3_7",
                  "AS1_3", "AS4_8", "AS9_10", "AS11",
                  "VC1_3_6", "VC4_5", 
                  "VD3_7", "VD8_11", "VD12",
                  "DD2_3", "DD4_5")

YA.2.MN <- YA.TPM.medium.mat[,colnames(YA.TPM.medium.mat) %in% neuron.order]
med.scaled.center <-t(scale(t(YA.2.MN)))
med.scaled.center[1:5,1:5]
# Z-scores, centered around 0.

# At first, I'll use the scaled, centered matrix
med.scaled.center <- as.data.frame(as.matrix(med.scaled.center))
med.scaled.center <- cbind(id = rownames(med.scaled.center), 
                           gene_name = i2s(rownames(med.scaled.center), gids), med.scaled.center)
med.scaled.center[1:5,1:5]
# Now I need to melt it into long form

library(reshape2)
?melt
med.scaled.long <- melt(med.scaled.center, id.vars = c("id", "gene_name"), 
                        variable.name = "cell.type", value.name = "scaled.expr")
head(med.scaled.long)

# Non-scaled
TPM.df <- data.frame(id = rownames(YA.2.MN), gene_name = i2s(rownames(YA.2.MN), gids), YA.2.MN)
TPM.df[1:5,1:5]
TPM.long <- melt(TPM.df, id.vars = c("id", "gene_name"),
                 variable.name = "cell.type", value.name = "TPM")

YA.prop <- YA_prop_by_type[,colnames(YA_prop_by_type) %in% neuron.order]
dim(YA.prop)
dim(YA.2.MN)
YA.prop <- YA.prop[,colnames(YA.2.MN)]

identical(colnames(YA.prop), colnames(YA.2.MN))
YA.prop[1:5,1:5]
med.prop <- data.frame(id = rownames(YA.prop), gene_name = i2s(rownames(YA.prop), gids), YA.prop)
med.prop <- med.prop[rownames(med.scaled.center),]
med.prop[1:5,1:5]
med.prop.long <- melt(med.prop, id.vars = c("id", "gene_name"), 
                      variable.name = "cell.type", value.name = "prop")
head(med.prop.long)

identical(TPM.long$id, med.prop.long$id)
identical(TPM.long$cell.type, med.prop.long$cell.type)

TPM.long$prop <- med.prop.long$prop*100

head(TPM.long)
TPM.long %>% filter(TPM == 0) %>% arrange(desc(prop)) %>% head(10)

class(TPM.long$prop)

TPM.long[which(TPM.long$TPM ==0), "prop"] <- 0

TPM.long %>% head(10)
TPM.long %>% filter(TPM == 0) %>% arrange(desc(prop)) %>% head(10)

table(TPM.long$TPM == 0 , TPM.long$prop == 0)

identical(med.scaled.long$id, med.prop.long$id)
identical(med.scaled.long$cell.type, med.prop.long$cell.type)

med.scaled.long$prop <- TPM.long$prop
head(med.scaled.long)

# for neuropeptides 
flp <- sort(grep("flp-", unique(med.scaled.long$gene_name), value = T))

flp.TPM <- YA.2.MN[intersect(rownames(YA.2.MN), s2i(flp, gids)),]

flp.TPM[1:5,1:5]
dim(flp.TPM)

flp.TPM <-flp.TPM[which(rowSums(flp.TPM != 0) > 0), ]
dim(flp.TPM)

flp.use <- sort(i2s(rownames(flp.TPM), gids))

flp.use

flp.order <- c("flp-1", "flp-2", "flp-3", "flp-4", "flp-6", "flp-7", "flp-9", "flp-11", "flp-12", "flp-13", 
               "flp-14", "flp-15","flp-16", "flp-17", "flp-18", "flp-19", "flp-20", "flp-22", "flp-23", "flp-24", "flp-26", 
               "flp-27", "flp-28", "flp-32")

# I'm not going to cluster the rows
med.scaled.long$cell.type <- factor(med.scaled.long$cell.type, levels = neuron.order)

flp.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% 
                                           flp.order),]

flp.scaled.long$gene_name <- factor(
  flp.scaled.long$gene_name, levels = rev(flp.order))

head(flp.scaled.long)

ggplot(flp.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# TPM 
TPM.long$cell.type <- factor(TPM.long$cell.type, levels = neuron.order)

flp.TPM.long <- TPM.long[which(TPM.long$gene_name %in% flp.order),]

flp.TPM.long$gene_name <- factor(
  flp.TPM.long$gene_name, levels = rev(flp.order))

head(flp.TPM.long)

ggplot(flp.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = TPM, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\n(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# log
ggplot(flp.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = log(TPM,2), size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\nlog2(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 


# nlp
nlp <- c(sort(grep("nlp-", unique(med.scaled.long$gene_name), value = T)), 
         c("pdf-1", "pdf-2", "ntc-1", "snet-1"))
nlp
nlp.TPM <- YA.2.MN[intersect(rownames(YA.2.MN), s2i(nlp, gids)),]
nlp.TPM[1:5,1:5]
dim(nlp.TPM)
nlp.TPM <- nlp.TPM[which(rowSums(nlp.TPM != 0) > 0), ]
dim(nlp.TPM)
nlp.use <- sort(i2s(rownames(nlp.TPM), gids))
nlp.use
# I'm not going to cluster the rows

nlp.order <- c("nlp-6", "nlp-12", "nlp-14", "nlp-15", "nlp-38", "nlp-58", "nlp-59", "nlp-62",
               "nlp-64", "nlp-71", "nlp-73", "snet-1", "ntc-1",  "pdf-1", "pdf-2", 
               "nlp-1", "nlp-5", "nlp-7", "nlp-8", "nlp-9",
               "nlp-11", "nlp-13", "nlp-16", "nlp-17", "nlp-20", "nlp-21", "nlp-24", "nlp-25",
               "nlp-27", "nlp-28", "nlp-29", "nlp-30", "nlp-33", "nlp-35", "nlp-36", "nlp-45",
               "nlp-48", "nlp-50", "nlp-51", "nlp-52", "nlp-53", "nlp-56", "nlp-57", "nlp-61",
               "nlp-63", "nlp-69", "nlp-70", "nlp-76", "nlp-77", "nlp-78", "nlp-79", "nlp-81")

nlp.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% nlp.order),]
nlp.scaled.long$gene_name <- factor(nlp.scaled.long$gene_name, levels = rev(nlp.order))

head(nlp.scaled.long)

ggplot(nlp.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# TPM 
nlp.TPM.long <- TPM.long[which(TPM.long$gene_name %in% nlp.order),]
nlp.TPM.long$gene_name <- factor(nlp.TPM.long$gene_name, levels = rev(nlp.order))

head(nlp.TPM.long)

ggplot(nlp.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = TPM, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\n(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# log
ggplot(nlp.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = log(TPM,2), size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\nlog2(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

receptors <- c(str_sort(unique(grep("frpr-", med.scaled.long$gene_name, value = T)), numeric = T),
               str_sort(unique(grep("npr-", med.scaled.long$gene_name, value =T)), numeric = T),
               str_sort(unique(grep("dmsr-", med.scaled.long$gene_name, value = T)), numeric = T),
               "tkr-1")
receptors.TPM <- YA.2.MN[intersect(rownames(YA.2.MN), s2i(receptors, gids)),]
receptors.TPM[1:5,1:5]
dim(receptors.TPM)
receptors.TPM <- receptors.TPM[which(rowSums(receptors.TPM != 0) > 0), ]
dim(receptors.TPM)
receptors.use <- i2s(rownames(receptors.TPM), gids)
receptors.use <- intersect(receptors, receptors.use)

# I'm not going to cluster the rows

receptors.order <- receptors.use

receptors.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% receptors.order),]
receptors.scaled.long$gene_name <- factor(receptors.scaled.long$gene_name, levels = rev(receptors.order))

head(receptors.scaled.long)

ggplot(receptors.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# TPM 
receptors.TPM.long <- TPM.long[which(TPM.long$gene_name %in% receptors.order),]
receptors.TPM.long$gene_name <- factor(receptors.TPM.long$gene_name, levels = rev(receptors.order))

head(receptors.TPM.long)

ggplot(receptors.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = TPM, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\n(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# log
ggplot(receptors.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = log(TPM,2), size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\nlog2(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# All receptors from Beets, et al., 
all.recep <- read.table("./Kratsios_scRNAseq/pep_receptors.csv", sep = ",", header = F, 
                        stringsAsFactors = F)
all.recep <- as.character(all.recep$V1)
all.recep.TPM <- YA.2.MN[intersect(rownames(YA.2.MN), s2i(all.recep, gids)),]
all.recep.TPM[1:5,1:5]
dim(all.recep.TPM)
all.recep.TPM <- all.recep.TPM[which(rowSums(all.recep.TPM != 0) > 0), ]
dim(all.recep.TPM)
all.recep.use <- i2s(rownames(all.recep.TPM), gids)
all.recep.use <- intersect(all.recep, all.recep.use)

# I'm not going to cluster the rows

all.recep.order <- str_sort(all.recep.use, numeric = T)

all.recep.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% all.recep.order),]
all.recep.scaled.long$gene_name <- factor(all.recep.scaled.long$gene_name, levels = rev(all.recep.order))

head(all.recep.scaled.long)

ggplot(all.recep.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# TPM 
all.recep.TPM.long <- TPM.long[which(TPM.long$gene_name %in% all.recep.order),]
all.recep.TPM.long$gene_name <- factor(all.recep.TPM.long$gene_name, levels = rev(all.recep.order))

head(all.recep.TPM.long)

ggplot(all.recep.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = TPM, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\n(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# log
ggplot(all.recep.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = log(TPM,2), size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\nlog2(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 


# All receptors from Beets, et al. without the frpr, npr, dmsr, trk-1 family
rest.recep <- setdiff(all.recep.order, receptors.order)

rest.recep.TPM <- YA.2.MN[intersect(rownames(YA.2.MN), s2i(rest.recep, gids)),]
rest.recep.TPM[1:5,1:5]
dim(rest.recep.TPM)
rest.recep.TPM <- rest.recep.TPM[which(rowSums(rest.recep.TPM != 0) > 0), ]
dim(rest.recep.TPM)
rest.recep.use <- i2s(rownames(rest.recep.TPM), gids)
rest.recep.use <- intersect(rest.recep, rest.recep.use)

# I'm not going to cluster the rows

rest.recep.order <- str_sort(rest.recep.use, numeric = T)

rest.recep.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% rest.recep.order),]
rest.recep.scaled.long$gene_name <- factor(rest.recep.scaled.long$gene_name, levels = rev(rest.recep.order))

head(rest.recep.scaled.long)

ggplot(rest.recep.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# TPM 
rest.recep.TPM.long <- TPM.long[which(TPM.long$gene_name %in% rest.recep.order),]
rest.recep.TPM.long$gene_name <- factor(rest.recep.TPM.long$gene_name, levels = rev(rest.recep.order))

head(rest.recep.TPM.long)

ggplot(rest.recep.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = TPM, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\n(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# log
ggplot(rest.recep.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = log(TPM,2), size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\nlog2(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 


# combined flp and nlp graphs
pep <- c(flp.order, nlp.order)

pep.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% pep),]
pep.scaled.long$gene_name <- factor(pep.scaled.long$gene_name, levels = rev(pep))

head(pep.scaled.long)

ggplot(pep.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# TPM 
pep.TPM.long <- TPM.long[which(TPM.long$gene_name %in% pep),]
pep.TPM.long$gene_name <- factor(pep.TPM.long$gene_name, levels = rev(pep))

head(pep.TPM.long)

ggplot(pep.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = TPM, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\n(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# log
ggplot(pep.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = log(TPM,2), size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\nlog2(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

saveRDS(med.scaled.long, "./Kratsios_scRNAseq/Analysis/060823_t2_scaled_expression_long_format.rds")
saveRDS(med.prop.long, "./Kratsios_scRNAseq/Analysis/060823_t2_proportion_unthresholded_long_format.rds")
saveRDS(TPM.long, "./Kratsios_scRNAseq/Analysis/060823_t2_TPM_long_format.rds")

# VD vs DD wt comparison

VD.DD.wt <- YA.neuron[,colData(VD.DD)$anatomical_class %in% c("VD", "DD")]
VD.DD.wt <- detect_genes(VD.DD.wt)
VD.DD.wt <- VD.DD.wt[rowData(VD.DD.wt)$num_cells_expressed > 5,]
VD.DD.wt
# 7587 features in 4193 cells

VD.DD.wt <- preprocess_cds(VD.DD.wt, method = "PCA", num_dim = 25)
plot_pc_variance_explained(VD.DD.wt)

VD.DD.wt <- reduce_dimension(VD.DD.wt, 
                             reduction_method = "UMAP",
                             umap.min_dist = 0.3,
                             umap.n_neighbors = 75)

VD.DD.wt <- cluster_cells(VD.DD.wt, res = 3e-3)

VD.DD.wt <- VD.DD.wt[,colData(VD.DD.wt)$Cell.type != "VD13"]

plot_cells(VD.DD.wt, 
           color_cells_by = "Cell.type",
           label_cell_groups = F,
           cell_size = 0.5)

plot_cells(VD.DD.wt, 
           color_cells_by = "anatomical_class",
           label_cell_groups = F,
           cell_size = 0.5)

colData(VD.DD.wt)$UMAP_1 <- reducedDims(VD.DD.wt)[["UMAP"]][,1]
colData(VD.DD.wt)$UMAP_2 <- reducedDims(VD.DD.wt)[["UMAP"]][,2]

colData(VD.DD.wt)$cluster <- monocle3::clusters(VD.DD.wt)

VD.DD.s <- CreateSeuratObject(counts = exprs(VD.DD.wt), meta.data = as.data.frame(colData(VD.DD.wt)))
Idents(VD.DD.s) <- VD.DD.s@meta.data$anatomical_class
table(VD.DD.s$anatomical_class)

DE <- FindMarkers(VD.DD.s, ident.1 = "DD", ident.2 = "VD")
DE$id <- rownames(DE)
DE$gene_short_name <- i2s(rownames(DE), gids)
DE.filt <- DE %>% filter(p_val_adj < 0.05)
head(DE.filt)

ggplot(DE, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point(color = "grey75") +
  geom_point(data = DE[which(DE$p_val_adj < 0.05),], color = "black") + 
  theme_classic() + 
  labs(x = "Avg log2FC", y = "-log10(p-value)") 

DE.filt %>% arrange(desc(abs(avg_log2FC))) %>% head(10)
plot.expr.UMAP(VD.DD.wt, "ZK380.6", size = 0.5)

write.table(DE.filt, "./Kratsios_scRNAseq/Analysis/032223_VD_DD_wild_type_DE_results.csv", sep = ",", col.names = T, row.names = F, quote = F)

#interesting genes between VD and DD (check across all VNC)
genes <- c("flp-13", "flp-14", "flp-24", "far-1", "del-1", "C09G1.5", "irx-1", "gopc-1", "odr-2", "W05E10.5", "T19D12.6", "glb-8", "C56G2.4", "cdh-4",
           "oig-1", "lntl-1", "msp-55", "hot-9", "R10E12.2", "unc-55", "plc-1", "lev-11","C34B2.9", "Y41C4A.17", "C05C12.6", "cblc-1", "F22B7.9", "R10E11.6")

genes.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% genes),]
genes.scaled.long$gene_name <- factor(genes.scaled.long$gene_name, levels = rev(genes))

head(genes.scaled.long)

ggplot(genes.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# TPM 
genes.TPM.long <- TPM.long[which(TPM.long$gene_name %in% genes),]
genes.TPM.long$gene_name <- factor(genes.TPM.long$gene_name, levels = rev(genes))

head(genes.TPM.long)

ggplot(genes.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = TPM, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\n(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# log
ggplot(genes.TPM.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = log(TPM,2), size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Expression\nlog2(TPM)", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

# List of genes that are expressed in all subclasses
colnames(YA.2.MN)

YA.2.MN[1:5,1:5]
dim(YA.2.MN)
colnames(t2)
t2 <- as.data.frame(as.matrix(YA.2.tpm))
t2 <- cbind(id = rownames(t2), gene_short_name = i2s(rownames(t2), gids), t2)
t2.MN <- t2 %>% select(c("gene_short_name", colnames(YA.2.MN)))
identical(rownames(YA.2.MN), rownames(t2.MN))
t2.MN$n_subclasses <- rowSums(YA.2.MN!=0)
head(t2.MN)

table(t2.MN$n_subclasses == 29)
# 1519 genes

identical(colnames(t2)[2:49], colnames(YA.TPM.medium.mat))

t2$n_classes <- rowSums(YA.TPM.medium.mat!=0)
table(t2$n_classes == 48)
# 383 genes that are in all the neurons tested (though some of these did not have as deep 
# of coverage as the MN classes).

# To add in the number of total classes to the number of MN subclasses

identical(rownames(t2), rownames(t2.MN))
t2.MN$n_classes <- rowSums(YA.TPM.medium.mat != 0)

pan.MN.genes <- t2.MN[which(t2.MN$n_subclasses == 29), c(1,31,32)]
head(pan.MN.genes)
pan.MN.genes <- pan.MN.genes %>% arrange(n_classes)
head(pan.MN.genes, 25)

ggplot(pan.MN.genes, aes(x = n_classes)) + geom_histogram(binwidth = 1) +
  theme_classic()

table(rowSums(YA.2.MN) == 0)
# I messed this up, and there are 1357 genes which aren't expressed in any of MNs.

colnames(YA.2.MN)

# All genes that are detected in thresholded data among MNs

MN.gene.list <- t2.MN %>% select(gene_short_name, n_subclasses, n_classes)
head(MN.gene.list)

colnames(MN.gene.list) <- c("gene_short_name", "n_MN_subclasses", "n_neuron_classes")

write.table(MN.gene.list, "./Kratsios_scRNAseq/Analysis/060923_thresholded_genes_expressed_in_MNs.csv",
            sep = ",", row.names = T, col.names = NA, quote = F)

wt.MN <- readRDS("./Kratsios_scRNAseq/060723_wt_VNC_cds.rds")
table(colData(wt.MN)$Cell.type)

MN.TFs <- read.table("./Kratsios_scRNAseq/Analysis/MN_TFs_FromGoAnalysis.csv", sep = ",",
                     header = T, stringsAsFactors = F)

new.tf <- read.table("./Kratsios_scRNAseq/Analysis/TFs in MN subclasses.csv", sep = ",",
                     header = T, stringsAsFactors = F)

head(new.tf)

new.tf$class <- new.tf$Category.2 %>% strsplit(": ") %>% sapply("[", 2)

table(new.tf$class)

identical(new.tf$Wormbase.ID, MN.TFs.use$Wormbase.ID)

head(MN.TFs)
MN.TFs$class <- MN.TFs$Category.2 %>% strsplit( ": ") %>% sapply("[", 2)
table(MN.TFs$class)

YA.2.tpm <- readRDS("./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t2_expression_matrix.rds")
YA.2.tpm <- as.matrix(YA.2.tpm)
MN.2tpm <- YA.2.tpm[, unique(colData(wt.MN)$Cell.type)]

dim(MN.2tpm)
n_subclasses <- data.frame(id = rownames(MN.2tpm), gene_short_name = i2s(rownames(MN.2tpm), gids),
                           n_classes = rowSums(MN.2tpm!=0))

table(n_subclasses$n_classes == 0)

MN.TFs$num_subclasses <- rowSums(MN.2tpm[MN.TFs$Wormbase.ID,]!=0)
head(MN.TFs)

table(MN.TFs$num_subclasses)
#  88 proteins that are showing up as not expressed in any subclasses

MN.tfs.use <- MN.TFs %>% filter(num_subclasses > 0) %>% 
  select(Wormbase.ID, class, num_subclasses) %>%
  arrange(class)

table(MN.tfs.use$class)

identical(sort(new.tf$Wormbase.ID), sort(MN.tfs.use$Wormbase.ID))

MN.tfs.sub <- MN.tfs.use %>% filter(num_subclasses < 29) %>%
  arrange(class)

table(MN.tfs.sub$class)

TF.order <- MN.tfs.sub$Wormbase.ID
TF.order <- i2s(TF.order, gids)

# Plotting this
med.scaled.long <- readRDS("./Kratsios_scRNAseq/Analysis/060823_t2_scaled_expression_long_format.rds")

TF.TPM <- YA.2.MN[intersect(rownames(YA.2.MN), s2i(TF.order, gids)),]
TF.TPM[1:5,1:5]
dim(TF.TPM)
TF.TPM <- TF.TPM[which(rowSums(TF.TPM != 0) > 0), ]
dim(TF.TPM)
TF.use <- i2s(rownames(TF.TPM), gids)

# I'm not going to cluster the rows
library(tidyverse)

TF.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% TF.order),]
TF.scaled.long$gene_name <- factor(TF.scaled.long$gene_name, levels = rev(TF.order))

head(TF.scaled.long)

ggplot(TF.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 

write.table(MN.tfs.sub, "./Kratsios_scRNAseq/Analysis/061523_MN_TF_list_from_Jayson_organized_by_class.rds")

# Possibly Longer list of TFs

tf.all <- read.table("./CeNGEN/L4_scRNAseq/final_gene_model/wTF3.csv", header = T, sep = ",", 
                     stringsAsFactors = F)

head(tf.all)
dim(YA.2.MN)

YA.2.MN.only <- YA.2.MN[rowSums(YA.2.MN)>0,]
dim(YA.2.MN.only)

tf.exprs <- intersect(tf.all$WBGeneID, rownames(YA.2.MN.only))

tf.exprs <- data.frame(id = tf.exprs, 
                       num_subclasses = rowSums(YA.2.MN.only[tf.exprs,]!= 0))

head(tf.exprs)
tf.all.keep <- tf.exprs %>% filter(num_subclasses < 29) %>% pull(id)

TF.all.TPM <- YA.2.MN[intersect(rownames(YA.2.MN), tf.all.keep),]
TF.all.TPM[1:5,1:5]
dim(TF.all.TPM)
TF.all.TPM <- TF.all.TPM[which(rowSums(TF.all.TPM != 0) > 0), ]
dim(TF.all.TPM)
TF.all.use <- i2s(rownames(TF.all.TPM), gids)

TF.all.order <- tf.all %>% filter(WBGeneID %in% s2i(TF.all.use, gids)) %>% pull(Public_name)
# I'm not going to cluster the rows
library(tidyverse)

TF.all.scaled.long <- med.scaled.long[which(med.scaled.long$gene_name %in% TF.all.order),]
TF.all.scaled.long$gene_name <- factor(TF.all.scaled.long$gene_name, levels = rev(TF.all.order))

head(TF.all.scaled.long)

ggplot(TF.all.scaled.long, aes(y = gene_name, x = cell.type)) + 
  geom_point(aes(color = scaled.expr, size = prop)) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust= 0.5, size = 8, hjust = 1), axis.text.y.left = element_text(size = 8),
        panel.background = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA)) + 
  scale_color_gradientn("Scaled \nExpression", colors = c("goldenrod1" ,"darkorange2", "maroon4", "navy")) + 
  scale_size_continuous(name = "Proportion", limit = c(0.02, 100), range = c(0,3)) + 
  labs(y = "", x= "") + theme(panel.grid = element_line(size = 0.1, color = "grey85")) 


# More files for the paper

# Excel sheet with all TFs expressed by class

# Most restricted genes
subclass.DE <- read.table("./Kratsios_scRNAseq/Analysis/060923_wt_VNC_within_class_subclass_DE.csv",
                          sep = ",", header = T, stringsAsFactors = F)

head(subclass.DE)
subclass.DE %>% arrange(pct_of_other_subclasses, desc(pct_of_subclass)) %>% head(20)

table(subclass.DE$pct_of_other_subclasses < 0.01, subclass.DE$subclass, exclude = NULL)
subclass.DE %>% filter(pct_of_other_subclasses < 0.01) %>% head(10)  

YA.t2 <- readRDS("./Kratsios_scRNAseq/Analysis/060823_YA_TPM_t2_expression_matrix.rds")
colnames(YA.t2)
MN.t2 <- YA.t2[,unique(subclass.DE$subclass)]

MN.t2[1:5,1:5]

YA.prop <- readRDS("./Kratsios_scRNAseq/Analysis/060823_YA_proportion_matrix.rds")
MN.prop <- YA.prop[,unique(subclass.DE$subclass)]

MN.prop.keep <- MN.prop[rownames(MN.t2),]
MN.t2.df <- data.frame(id = rownames(MN.t2),
                       gene_short_name = i2s(rownames(MN.t2), gids),
                       n.cells = rowSums(MN.t2 != 0))
head(MN.t2.df)

table(MN.t2.df$n.cells == 1)
MN.t2.df$max.TPM <- apply(MN.t2, 1, max)

max.TPM.loc <- apply(MN.t2, 1, which.max)
max.prop.loc <- apply(MN.prop.keep, 1, which.max)
identical(max.TPM.loc, max.prop.loc)

dim(MN.t2)
max.type <- data.frame(max.TPM.loc = max.TPM.loc,
                       max.TPM.type = colnames(MN.t2)[max.TPM.loc],
                       max.prop.loc = max.prop.loc,
                       max.prop.type = colnames(MN.prop.keep)[max.prop.loc])

head(max.type)

MN.t2.df$max.TPM.type <- max.type$max.TPM.type

MN.t2.df$max.prop <- apply(MN.prop.keep, 1, max)
MN.t2.df$max.prop.type <- max.type$max.prop.type

MN.t2.df %>% filter(n.cells > 0) %>% arrange(n.cells) %>% head(20)

MN.t2.df %>% filter(n.cells > 0) %>% arrange(n.cells, desc(max.TPM)) %>% head(20)

table(MN.t2.df$n.cells == 1, MN.t2.df$max.type)

MN.t2.df <- MN.t2.df %>% filter(n.cells > 0) %>% arrange(n.cells, desc(max.prop))
saveRDS(MN.t2.df, "./Kratsios_scRNAseq/Analysis/070723_gene_list_all_MN_subclasses_n_subclasses_detected.rds")

MN.t2.df.filt <- MN.t2.df %>% filter(n.cells == 1 & max.prop > 0.1) %>% arrange(max.prop.type, desc(max.prop))
table(MN.t2.df.filt$max.prop.type)

head(MN.t2.df.filt)

# Doing this for each neuron class

n.classes <- function(neuron.type) {
  subclasses <- subclass.DE %>% filter(anatomical_class == neuron.type) %>% pull(subclass)
  subclasses <- unique(subclasses)
  sub.TPM <- MN.t2[,subclasses]
  sub.prop <- MN.prop.keep[,subclasses]
  sub.df <- data.frame(id = rownames(sub.TPM),
                       gene_short_name = i2s(rownames(sub.TPM), gids),
                       n.cells = rowSums(sub.TPM != 0),
                       max.TPM = apply(sub.TPM, 1, max),
                       max.TPM.type = colnames(sub.TPM)[apply(sub.TPM,1,which.max)],
                       max.prop = apply(sub.prop, 1, max),
                       max.prop.type = colnames(sub.prop)[apply(sub.prop, 1, which.max)])
  sub.df <- sub.df %>% filter(n.cells > 0) %>% arrange(n.cells, desc(max.prop))
  sub.df
}

AS.df <- n.classes("AS")
head(AS.df)

DA.df <- n.classes("DA")

DB.df <- n.classes("DB")

DD.df <- n.classes("DD")
VA.df <- n.classes("VA")
VB.df <- n.classes("VB")
VC.df <- n.classes("VC")
VD.df <- n.classes("VD")
head(VD.df)

combined.df <- rbind(AS.df, DA.df, DB.df, DD.df, VA.df, VB.df, VC.df, VD.df)
combined.df <- combined.df %>% filter(n.cells == 1) %>% arrange(max.prop.type, desc(max.prop)) 
table(combined.df$max.prop.type)

summary(combined.df$max.prop)
table(combined.df$max.prop > 0.1)

combined.df.filt <- combined.df %>% filter(max.prop > 0.1)
table(combined.df.filt$max.prop.type)
head(combined.df.filt, 20)
MN.t2[s2i("egl-44", gids),]
MN.prop.keep[s2i("egl-44", gids),]
table(MN.t2.df.filt$max.prop.type)

identical(combined.df.filt$max.TPM.type, combined.df.filt$max.prop.type)

combined.df.filt$anatomical_class <- ifelse(
  combined.df.filt$max.prop.type %in% c("AS1_3", "AS11", "AS9_10"),
  "AS",
  ifelse(
    combined.df.filt$max.prop.type %in% c("DA1", "DA2_5", "DA6_8", "DA9"),
    "DA",
    ifelse(
      combined.df.filt$max.prop.type %in% c("DB1", "DB2", "DB3_7"),
      "DB",
      ifelse(
        combined.df.filt$max.prop.type %in% c("DD2_3", "DD4_5"),
        "DD",
        ifelse(
          combined.df.filt$max.prop.type %in% c("VA1", "VA11", "VA12", "VA2", "VA3_8", "VA9_10"),
          "VA",
          ifelse(
            combined.df.filt$max.prop.type %in% c("VB1", "VB10_11", "VB2", "VB3", "VB4_9"),
            "VB",
            ifelse(
              combined.df.filt$max.prop.type %in% c("VC1_3_6", "VC4_5"),
              "VC",
              "VD"
            )
          )
        )
      )
    )
  )
)

table(combined.df.filt$anatomical_class, combined.df.filt$max.prop.type, exclude = NULL)

saveRDS(combined.df, "./Kratsios_scRNAseq/Analysis/070723_gene_list_n_subclass_per_anatomical_class_combined.rds")
saveRDS(combined.df.filt ,"./Kratsios_scRNAseq/Analysis/070723_gene_list_n_subclass_per_anatomical_class_filtered.rds")
saveRDS(MN.t2.df.filt, "./Kratsios_scRNAseq/Analysis/070723_single_subclass_driver_list.rds")

colnames(combined.df.filt)
colnames(MN.t2.df.filt)

write.table(combined.df.filt, "./Kratsios_scRNAseq/Analysis/071123_gene_list_single_subclass_within_anatomical_class.csv",
            sep = ",", quote = F, row.names = F, col.names = T)

write.table(MN.t2.df.filt, "./Kratsios_scRNAseq/Analysis/071123_gene_list_single_subclass_whole_VNC.csv",
            sep = ",", quote = F, row.names = F, col.names = T)

MN.TFs <- read.table("./Kratsios_scRNAseq/Analysis/MN_TFs_FromGoAnalysis.csv", sep = ",",
                     header = T, stringsAsFactors = F)

MN.TFs$class <- MN.TFs$Category.2 %>% strsplit( ": ") %>% sapply("[", 2)
table(MN.TFs$class)
MN.TFs$num_subclasses <- rowSums(MN.t2[MN.TFs$Wormbase.ID,]!=0)
head(MN.TFs)

table(MN.TFs$num_subclasses)
# Why are there 88 proteins that are showing up as not expressed in any subclasses?

MN.tfs.use <- MN.TFs %>% filter(num_subclasses > 0) %>% 
  select(Wormbase.ID, class, num_subclasses) %>%
  arrange(class)

table(MN.tfs.use$class)

head(MN.tfs.use)
length(intersect(MN.tfs.use$Wormbase.ID, rownames(MN.t2)))

MN.TF.t2 <- MN.t2[MN.tfs.use$Wormbase.ID, ]

MN.TF.df <- data.frame(id = rownames(MN.TF.t2),
                       gene_symbol = i2s(rownames(MN.TF.t2), gids),
                       class = MN.tfs.use$class,
                       as.matrix(MN.TF.t2))

head(MN.TF.df)

write.table(MN.TF.df, "./Kratsios_scRNAseq/Analysis/070723_MN_TF_expression_table.csv",
            row.names = F, col.names = T, quote = F, sep = ",")

# now with all TFs from Walhout group
tfs <- read.table("./CeNGEN/L4_scRNAseq/final_gene_model/wTF3.csv", header = T, stringsAsFactors = F,
                  sep = ",")
head(tfs)
table(tfs$class)

tfs.use <- intersect(tfs$WBGeneID, rownames(MN.t2))
tfs.df <- tfs %>% filter(WBGeneID %in% tfs.use)
tfs.df$num_subclasses <- rowSums(MN.t2[tfs.use,]!=0)
head(tfs.df)

table(tfs.df$num_subclasses)
# There are 106 TFs that are in 0

tfs.new.use <- tfs.df %>% filter(num_subclasses > 0) %>% pull(WBGeneID)

length(intersect(tfs.new.use, rownames(MN.t2)))

MN.allTF.t2 <- MN.t2[tfs.new.use, ]

MN.allTF.df <- data.frame(id = rownames(MN.allTF.t2),
                          gene_symbol = i2s(rownames(MN.allTF.t2), gids),
                          class = tfs.df[which(tfs.df$WBGeneID %in% tfs.new.use),]$Domain_Family,
                          as.matrix(MN.allTF.t2))

head(MN.allTF.df)
table(MN.allTF.df$class)

write.table(MN.allTF.df, "./Kratsios_scRNAseq/Analysis/070723_MN_allTF_expression_table.csv",
            row.names = F, col.names = T, quote = F, sep = ",")

L4.updated <- readRDS("./CeNGEN/L1_scRNAseq/8672-ST/SoupX/112122_L4_CeNGEN_updated_annotations_UPR_separate_cds.rds")

L4.updated.data <- as.data.frame(colData(L4.updated)) %>% group_by(Cell.type)
L4.updated.stats <- L4.updated.data %>% summarize(
  n = n(),
  median.umi = median(total_counts),
  median.ngenes = median(num_genes_expressed)
)

L4.updated.stats %>% filter(Cell.type %in% c("AS", "AS11", "DA1", "DA", "DA9", "DB", "DB01", "DB02",
                                             "DD", "VA1", "VA", "VA12", "VB01", "VB02", "VB", "VC", 
                                             "VC_4_5", "VD"))

L4.MN <- L4.updated[,colData(L4.updated)$Cell.type %in% c("AS", "AS11", "DA1", "DA", "DA9", "DB", "DB01", "DB02",
                                                          "DD", "VA1", "VA", "VA12", "VB01", "VB02", "VB", "VC", 
                                                          "VC_4_5", "VD")]

L4.MN <- detect_genes(L4.MN)
L4.MN <- L4.MN[rowData(L4.MN)$num_cells_expressed > 5,]
L4.MN <- preprocess_cds(L4.MN, num_dim = 22)
plot_pc_variance_explained(L4.MN)

L4.MN <- align_cds(L4.MN, alignment_group = 'Experiment', alignment_k = 5)

L4.MN <- reduce_dimension(L4.MN, 
                          reduction_method = "UMAP",
                          preprocess_method = "Aligned",
                          umap.min_dist = 0.1, 
                          umap.n_neighbors = 15)

L4.MN <- cluster_cells(L4.MN, res = 3e-3)

plot_cells(L4.MN, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

source("./Monocle3Functions.txt")

colData(L4.MN)$UMAP_1 <- reducedDims(L4.MN)[["UMAP"]][,1]
colData(L4.MN)$UMAP_2 <- reducedDims(L4.MN)[["UMAP"]][,2]

# Just subsets

AS.MN <- L4.updated[,colData(L4.updated)$Cell.type %in% c("DA1", "DA", "DA9", "VA1", "VA", "VA12")]

AS.MN <- detect_genes(AS.MN)
AS.MN <- AS.MN[rowData(AS.MN)$num_cells_expressed > 5,]
AS.MN <- preprocess_cds(AS.MN, num_dim = 20)
plot_pc_variance_explained(AS.MN)

AS.MN <- align_cds(AS.MN, alignment_group = 'Experiment', alignment_k = 5)

AS.MN <- reduce_dimension(AS.MN, 
                          reduction_method = "UMAP",
                          preprocess_method = "PCA",
                          umap.min_dist = 0.1, 
                          umap.n_neighbors = 15)

plot_cells(AS.MN, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

colData(AS.MN)$UMAP_1 <- reducedDims(AS.MN)[["UMAP"]][,1]
colData(AS.MN)$UMAP_2 <- reducedDims(AS.MN)[["UMAP"]][,2]

plot.expr.UMAP(AS.MN, "mab-5", size = 0.5)
plot.expr.UMAP(AS.MN, "lin-39", size = 0.5)
plot.expr.UMAP(AS.MN, "egl-5", size = 0.5)
plot.expr.UMAP(AS.MN, "ceh-13", size = 0.5)
plot.expr.UMAP(AS.MN, "bnc-1", size = 0.5)
plot.expr.UMAP(AS.MN, "unc-129", size = 0.5)
plot.expr.UMAP(AS.MN, "unc-4", size = 0.5)
plot.expr.UMAP(AS.MN, "elt-1", size = 0.5)
plot.expr.UMAP(AS.MN, "nlp-38", size = 0.5)
plot.expr.UMAP(AS.MN, "cav-1", size = 0.5)
plot.expr.UMAP(AS.MN, "npr-2", size = 0.5)
plot.expr.UMAP(AS.MN, "EGAP4.1", size = 0.5)
plot.expr.UMAP(AS.MN, "F26F12.8", size = 0.5)
plot.expr.UMAP(AS.MN, "W05E10.5", size = 0.5)
plot.expr.UMAP(AS.MN, "nlp-7", size = 0.5)
plot.expr.UMAP(AS.MN, "T02C12.5", size = 0.5)
plot.expr.UMAP(AS.MN, "F02E9.10", size = 0.5)
plot.expr.UMAP(AS.MN, "lntl-1", size = 0.5)

plot.expr.UMAP(AS.MN, "vab-3", size = 0.5)
plot.expr.UMAP(AS.MN, "acr-20", size = 0.5)
plot.expr.UMAP(AS.MN, "nlp-11", size = 0.5)
plot.expr.UMAP(AS.MN, "egl-44", size = 0.5)

AS.MN <- cluster_cells(AS.MN, res = 3e-2)

plot_cells(AS.MN, 
           color_cells_by = "cluster",
           group_label_size = 4,
           label_groups_by_cluster = T,
           cell_size = 0.5)

colData(AS.MN)$cluster <- monocle3::clusters(AS.MN)

plot.expr.UMAP(AS.MN, "unc-129", size = 0.5) + 
  geom_hline(yintercept = -2.2) + 
  geom_vline(xintercept = -0.2)

colData(AS.MN)$Cell.type <- ifelse(
  colData(AS.MN)$cluster == 15 &
    colData(AS.MN)$UMAP_1 < -0.2 & 
    colData(AS.MN)$UMAP_2 < -2.2,
  "VA1",
  as.character(colData(AS.MN)$Cell.type)
)

colData(AS.MN)$Cell.type <- ifelse(
  colData(AS.MN)$cluster == 15 &
    colData(AS.MN)$Cell.type != "VA1",
  "DA1",
  as.character(colData(AS.MN)$Cell.type)
)

plot_cells(AS.MN, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

# Can I distinguish posterior classes at all?
# VA11
plot.expr.UMAP(AS.MN, "vab-3", size = 0.5)
plot.expr.UMAP(AS.MN, "nlp-1", size = 0.5)
plot.expr.UMAP(AS.MN, "mrpl-24", size = 0.5)
plot.expr.UMAP(AS.MN, "nep-26", size = 0.5)
plot.expr.UMAP(AS.MN, "nhr-199", size = 0.5)
plot.expr.UMAP(AS.MN, "avr-15", size = 0.5)
# A few cells here or there expressing some of these genes, but not enough to cluster separately.
# So either VA11 isn't that different at L4 or we don't have enough of them.

# What about more VA9-10, DA6-8?
plot.expr.UMAP(AS.MN, "acr-20", size = 0.5)
plot.expr.UMAP(AS.MN, "nlp-11", size = 0.5)
plot.expr.UMAP(AS.MN, "mab-5", size = 0.5)
# there are clear mab-5/nlp-11 subsets.

plot.expr.UMAP(AS.MN, "lin-39", size = 0.5)
plot.expr.UMAP(AS.MN, "egl-5", size = 0.5)
plot.expr.UMAP(AS.MN, "ceh-13", size = 0.5)
plot.expr.UMAP(AS.MN, "bnc-1", size = 0.5)
plot.expr.UMAP(AS.MN, "unc-129", size = 0.5)
plot.expr.UMAP(AS.MN, "mec-1", size = 0.5)
plot.expr.UMAP(AS.MN, "T22E5.1", size = 0.5)
plot.expr.UMAP(AS.MN, "cwn-1", size = 0.5)
plot.expr.UMAP(AS.MN, "cyk-7", size = 0.5)
plot.expr.UMAP(AS.MN, "nep-21", size = 0.5)
plot.expr.UMAP(AS.MN, "lgc-38", size = 0.5)
plot.expr.UMAP(AS.MN, "lgc-35", size = 0.5)

# There are distinctions between the mab-5 and the lin-39 

colData(AS.MN)$Cell.type <- ifelse(
  colData(AS.MN)$cluster %in% c(20, 14),
  "VA9_11",
  as.character(colData(AS.MN)$Cell.type)
)

AS.sub <- AS.MN[,colData(AS.MN)$cluster %in% c(9,22,14,20,12,2,13,10,15,18,4,8)]

AS.sub <- detect_genes(AS.sub)
AS.sub <- AS.sub[rowData(AS.sub)$num_cells_expressed > 5,]
AS.sub <- preprocess_cds(AS.sub, num_dim = 15)
plot_pc_variance_explained(AS.sub)

AS.sub <- align_cds(AS.sub, alignment_group = 'Experiment', alignment_k = 5)

AS.sub <- reduce_dimension(AS.sub, 
                           reduction_method = "UMAP",
                           preprocess_method = "PCA",
                           umap.min_dist = 0.1, 
                           umap.n_neighbors = 15)

plot_cells(AS.sub, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

colData(AS.sub)$UMAP_1 <- reducedDims(AS.sub)[["UMAP"]][,1]
colData(AS.sub)$UMAP_2 <- reducedDims(AS.sub)[["UMAP"]][,2]

plot.expr.UMAP(AS.sub, "mab-5", size = 0.5)
plot.expr.UMAP(AS.sub, "lin-39", size = 0.5)
plot.expr.UMAP(AS.sub, "ceh-13", size = 0.5)
plot.expr.UMAP(AS.sub, "bnc-1", size = 0.5)
plot.expr.UMAP(AS.sub, "unc-129", size = 0.5)
plot.expr.UMAP(AS.sub, "unc-4", size = 0.5)
plot.expr.UMAP(AS.sub, "elt-1", size = 0.5)
plot.expr.UMAP(AS.sub, "nlp-38", size = 0.5)
plot.expr.UMAP(AS.sub, "cav-1", size = 0.5)
plot.expr.UMAP(AS.sub, "npr-2", size = 0.5)
plot.expr.UMAP(AS.sub, "EGAP4.1", size = 0.5)
plot.expr.UMAP(AS.sub, "F26F12.8", size = 0.5)
plot.expr.UMAP(AS.sub, "W05E10.5", size = 0.5)
plot.expr.UMAP(AS.sub, "nlp-7", size = 0.5)
plot.expr.UMAP(AS.sub, "T02C12.5", size = 0.5)
plot.expr.UMAP(AS.sub, "F02E9.10", size = 0.5)
plot.expr.UMAP(AS.sub, "lntl-1", size = 0.5)

plot.expr.UMAP(AS.sub, "vab-3", size = 0.5)
plot.expr.UMAP(AS.sub, "acr-20", size = 0.5)
plot.expr.UMAP(AS.sub, "nlp-11", size = 0.5)

AS.sub <- cluster_cells(AS.sub, res = 3e-1)

plot_cells(AS.sub, 
           color_cells_by = "cluster",
           group_label_size = 4,
           label_groups_by_cluster = T,
           cell_size = 0.5)

colData(AS.sub)$cluster <- monocle3::clusters(AS.sub)

colData(AS.sub)$Cell.type <- ifelse(
  colData(AS.sub)$cluster %in% c(5),
  "VA1",
  as.character(colData(AS.sub)$Cell.type)
)

colData(AS.sub)$Cell.type <- ifelse(
  colData(AS.sub)$cluster %in% c(23) & colData(AS.sub)$UMAP_2 > -1.7,
  "AS",
  as.character(colData(AS.sub)$Cell.type)
)

colData(AS.sub)$Cell.type <- ifelse(
  colData(AS.sub)$cluster %in% c(23) & colData(AS.sub)$UMAP_2 < -1.7,
  "DA",
  as.character(colData(AS.sub)$Cell.type)
)

# Can I distinguish posterior classes at all?
# VA11
plot.expr.UMAP(AS.sub, "nlp-1", size = 0.5)
plot.expr.UMAP(AS.sub, "mrpl-24", size = 0.5)
plot.expr.UMAP(AS.sub, "nep-26", size = 0.5)

# A few cells here or there expressing some of these genes, but not enough to cluster separately.
# So either VA11 isn't that different at L4 or we don't have enough of them.

# What about more VA9-10, DA6-8?
plot.expr.UMAP(AS.sub, "acr-20", size = 0.5)
plot.expr.UMAP(AS.sub, "nlp-11", size = 0.5)
plot.expr.UMAP(AS.sub, "mab-5", size = 0.5)
# there are clear mab-5/nlp-11 subsets.

plot.expr.UMAP(AS.sub, "lin-39", size = 0.5)
plot.expr.UMAP(AS.sub, "ceh-13", size = 0.5)
plot.expr.UMAP(AS.sub, "bnc-1", size = 0.5)
plot.expr.UMAP(AS.sub, "unc-129", size = 0.5)
plot.expr.UMAP(AS.sub, "mec-1", size = 0.5)
plot.expr.UMAP(AS.sub, "cwn-1", size = 0.5)
plot.expr.UMAP(AS.sub, "cyk-7", size = 0.5)
plot.expr.UMAP(AS.sub, "nep-21", size = 0.5)
plot.expr.UMAP(AS.sub, "lgc-38", size = 0.5)
plot.expr.UMAP(AS.sub, "lgc-35", size = 0.5)

# There are distinctions between the mab-5 and the lin-39 

plot_cells(AS.sub, 
           color_cells_by = "cluster",
           group_label_size = 4,
           label_groups_by_cluster = T,
           cell_size = 0.5)

colData(AS.sub)$Cell.type <- ifelse(
  colData(AS.sub)$cluster %in% c(11,16),
  "VA9_11",
  as.character(colData(AS.sub)$Cell.type)
)

colData(AS.sub)$Cell.type <- ifelse(
  colData(AS.sub)$cluster %in% c(13,22,25,35,6,20,9),
  "VA",
  as.character(colData(AS.sub)$Cell.type)
)

colData(AS.sub)$Cell.type <- ifelse(
  colData(AS.sub)$cluster %in% c(9) & colData(AS.sub)$UMAP_1 < 1.2,
  "VA9_11",
  as.character(colData(AS.sub)$Cell.type)
)

colData(AS.sub)$Cell.type <- ifelse(
  colData(AS.sub)$cluster %in% c(26,10,19,1),
  "DA6_8",
  as.character(colData(AS.sub)$Cell.type)
)

# VA2?
plot.expr.UMAP(AS.sub, "flp-16", size = 0.5)
plot.expr.UMAP(AS.sub, "C24H12.1", size = 0.5)
plot.expr.UMAP(AS.sub, "pat-9", size = 0.5)
plot.expr.UMAP(AS.sub, "acs-20", size = 0.5)
plot.expr.UMAP(AS.sub, "lec-3", size = 0.5)
plot.expr.UMAP(AS.sub, "dma-1", size = 0.5)
plot.expr.UMAP(AS.sub, "tnt-2", size = 0.5)

# possibly a slightly more anterior group, but not clearly separate.
# I could rename DA to be DA2_5, since that is all that is left. 

colData(AS.sub)$Cell.type <- dplyr::recode(colData(AS.sub)$Cell.type, "DA" = "DA2_5")
colData(AS.MN)$Cell.type <- dplyr::recode(colData(AS.MN)$Cell.type, "DA" = "DA2_5")

colData(AS.MN)[colnames(AS.sub),]$Cell.type <- colData(AS.sub)$Cell.type
colData(L4.MN)[colnames(AS.MN),]$Cell.type <- colData(AS.MN)$Cell.type
colData(L4.updated)[colnames(L4.MN),]$Cell.type <- colData(L4.MN)$Cell.type

table(colData(L4.MN)$Cell.type, exclude = NULL)

# AS (2/3 or 9/10?), DB (mab-5 containing? DB7?), DD (4/5?), VB (VB3?, VB10?), VD (mab-5?)

AS.DB.DD <- L4.updated[,colData(L4.updated)$Cell.type %in% c("AS", "DB", "DD")]
AS.DB.DD <- detect_genes(AS.DB.DD)
AS.DB.DD <- AS.DB.DD[rowData(AS.DB.DD)$num_cells_expressed > 5,]
AS.DB.DD <- preprocess_cds(AS.DB.DD, num_dim = 8)
plot_pc_variance_explained(AS.DB.DD)

AS.DB.DD <- align_cds(AS.DB.DD, alignment_group = 'Experiment', alignment_k = 5)

AS.DB.DD <- reduce_dimension(AS.DB.DD, 
                             reduction_method = "UMAP",
                             preprocess_method = "PCA",
                             umap.min_dist = 0.1, 
                             umap.n_neighbors = 15)

plot_cells(AS.DB.DD, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

colData(AS.DB.DD)$UMAP_1 <- reducedDims(AS.DB.DD)[["UMAP"]][,1]
colData(AS.DB.DD)$UMAP_2 <- reducedDims(AS.DB.DD)[["UMAP"]][,2]

plot.expr.UMAP(AS.DB.DD, "mab-5", size = 0.5)
plot.expr.UMAP(AS.DB.DD, "lin-39", size = 0.5)
plot.expr.UMAP(AS.DB.DD, "ceh-13", size = 0.5)
plot.expr.UMAP(AS.DB.DD, "npr-2", size = 0.5)
plot.expr.UMAP(AS.DB.DD, "unc-25", size = 0.5)
plot.expr.UMAP(AS.DB.DD, "vab-7", size = 0.5)
plot.expr.UMAP(AS.DB.DD, "egl-44", size = 0.5)

# only one AS cluster seems to have a strong mab-5/lin-39 difference
# mab-5 in DB and DD is quite spread out. Now just AS

AS <- L4.updated[,colData(L4.updated)$Cell.type %in% c("AS", "AS11")]
AS <- detect_genes(AS)
AS <- AS[rowData(AS)$num_cells_expressed > 5,]
AS <- preprocess_cds(AS, num_dim = 8)
plot_pc_variance_explained(AS)

AS <- align_cds(AS, alignment_group = 'Experiment', alignment_k = 5)

AS <- reduce_dimension(AS, 
                       reduction_method = "UMAP",
                       preprocess_method = "PCA",
                       umap.min_dist = 0.1, 
                       umap.n_neighbors = 15)

plot_cells(AS, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

colData(AS)$UMAP_1 <- reducedDims(AS)[["UMAP"]][,1]
colData(AS)$UMAP_2 <- reducedDims(AS)[["UMAP"]][,2]

plot.expr.UMAP(AS, "mab-5", size = 0.5)
plot.expr.UMAP(AS, "lin-39", size = 0.5)
plot.expr.UMAP(AS, "ceh-13", size = 0.5)
plot.expr.UMAP(AS, "vab-3", size = 0.5)
# There is one subset that has mab-5, some of which also has vab-3. Other AS11 markers?

plot.expr.UMAP(AS, "flp-9", size = 0.5)
plot.expr.UMAP(AS, "flp-7", size = 0.5)
plot.expr.UMAP(AS, "glb-32", size = 0.5)
plot.expr.UMAP(AS, "lin-1", size = 0.5)
plot.expr.UMAP(AS, "avr-15", size = 0.5)
plot.expr.UMAP(AS, "sol-1", size = 0.5)
plot.expr.UMAP(AS, "unc-129", size = 0.5)
plot.expr.UMAP(AS, "kcnl-2", size = 0.5)
plot.expr.UMAP(AS, "C50D2.6", size = 0.5)
plot.expr.UMAP(AS, "ctbp-1", size = 0.5)

# AS9-10 markers
plot.expr.UMAP(AS, "acr-20", size = 0.5)
plot.expr.UMAP(AS, "ser-6", size = 0.5)
plot.expr.UMAP(AS, "F54F12.2", size = 0.5)
plot.expr.UMAP(AS, "ZK1073.2", size = 0.5)
plot.expr.UMAP(AS, "insc-1", size = 0.5)
plot.expr.UMAP(AS, "smp-1", size = 0.5)
plot.expr.UMAP(AS, "Y47A7.2", size = 0.5)
plot.expr.UMAP(AS, "nlp-11", size = 0.5)

AS <- cluster_cells(AS, res = 3e-2)

plot_cells(AS, 
           color_cells_by = "cluster",
           group_label_size = 4,
           label_groups_by_cluster = T,
           cell_size = 0.5)

colData(AS)$cluster <- monocle3::clusters(AS)

colData(AS)$Cell.type <- ifelse(
  colData(AS)$cluster == 5,
  "AS11",
  as.character(colData(AS)$Cell.type)
)

# AS2-3 markers
plot.expr.UMAP(AS, "EGAP4.1", size = 0.5)
plot.expr.UMAP(AS, "egl-44", size = 0.5)
plot.expr.UMAP(AS, "C56G2.4", size = 0.5)
plot.expr.UMAP(AS, "tsp-8", size = 0.5)
plot.expr.UMAP(AS, "F26F12.8", size = 0.5)
plot.expr.UMAP(AS, "msa-1", size = 0.5)
plot.expr.UMAP(AS, "Y105C5B.25", size = 0.5)
plot.expr.UMAP(AS, "F26F12.8", size = 0.5)

# No clear AS2-3 separation. Markers are expressed in cells, but they are scattered around.
# So the "AS" cluster may contain AS1-10.

colData(L4.MN)[colnames(AS),]$Cell.type <- colData(AS)$Cell.type
colData(L4.updated)[colnames(AS),]$Cell.type <- colData(AS)$Cell.type

# VB and VD

v <- L4.MN[,colData(L4.MN)$Cell.type %in% c("VB", "VD")]
v <- detect_genes(v)
v <- v[rowData(v)$num_cells_expressed > 5,]
v <- preprocess_cds(v, num_dim = 10)
plot_pc_variance_explained(v)

v <- align_cds(v, alignment_group = 'Experiment', alignment_k = 5)

v <- reduce_dimension(v, 
                      reduction_method = "UMAP",
                      preprocess_method = "PCA",
                      umap.min_dist = 0.1, 
                      umap.n_neighbors = 15)

plot_cells(v, 
           color_cells_by = "Cell.type",
           group_label_size = 4,
           label_groups_by_cluster = F,
           cell_size = 0.5)

colData(v)$UMAP_1 <- reducedDims(v)[["UMAP"]][,1]
colData(v)$UMAP_2 <- reducedDims(v)[["UMAP"]][,2]

# VB10-11

plot.expr.UMAP(v, "mab-5", size = 0.5)
plot.expr.UMAP(v, "lin-39", size = 0.5)
plot.expr.UMAP(v, "ceh-13", size = 0.5)
plot.expr.UMAP(v, "vab-3", size = 0.5)
plot.expr.UMAP(v, "F38E9.6", size = 0.5)
plot.expr.UMAP(v, "rmh-1", size = 0.5)
plot.expr.UMAP(v, "cwn-1", size = 0.5)
plot.expr.UMAP(v, "bath-15", size = 0.5)
plot.expr.UMAP(v, "cyk-7", size = 0.5)
plot.expr.UMAP(v, "F26F12.8", size = 0.5)
plot.expr.UMAP(v, "nlp-59", size = 0.5)
plot.expr.UMAP(v, "glb-23", size = 0.5)
plot.expr.UMAP(v, "acr-20", size = 0.5)
plot.expr.UMAP(v, "linc-22", size = 0.5)
plot.expr.UMAP(v, "srj-29", size = 0.5)
# few cells expressing - either late maturing or didn't capture

# VB3
plot.expr.UMAP(v, "C24H12.1", size = 0.5)
plot.expr.UMAP(v, "F22F1.2", size = 0.5)
plot.expr.UMAP(v, "F53C11.9", size = 0.5)
plot.expr.UMAP(v, "ham-1", size = 0.5)
plot.expr.UMAP(v, "pat-9", size = 0.5)
plot.expr.UMAP(v, "M02G9.2", size = 0.5)
plot.expr.UMAP(v, "ceh-13", size = 0.5)

# VB3 present, but not a distinct subset

# VD8-11
plot.expr.UMAP(v, "lin-39", size = 0.5)
plot.expr.UMAP(v, "mab-5", size = 0.5)
plot.expr.UMAP(v, "ceh-13", size = 0.5)
plot.expr.UMAP(v, "nlp-11", size = 0.5)
plot.expr.UMAP(v, "T13H5.1", size = 0.5)
plot.expr.UMAP(v, "unc-77", size = 0.5)
plot.expr.UMAP(v, "ptr-6", size = 0.5)
plot.expr.UMAP(v, "T23F2.3", size = 0.5)

# possible anterior markers
plot.expr.UMAP(v, "exp-2", size = 0.5)
plot.expr.UMAP(v, "str-155", size = 0.5)
plot.expr.UMAP(v, "tyra-2", size = 0.5)
plot.expr.UMAP(v, "F26A1.19", size = 0.5)
# No clear VD subsets that I can distinguish

table(colData(L4.MN)$Cell.type)
table(colData(MN.cds)$Cell.type)

# Matching names where possible
colData(L4.MN)$Cell.type <- dplyr::recode(colData(L4.MN)$Cell.type, "DB01" = "DB1",
                                          "DB02" = "DB2",
                                          "DB" = "DB3_7",
                                          "VB01" = "VB1", 
                                          "VB02" = "VB2",
                                          "VC" = "VC1_3_6",
                                          "VC_4_5" = "VC4_5")

table(colData(L4.MN)$Cell.type)
table(colData(MN.cds)$Cell.type)

colData(L4.updated)[colnames(L4.MN),]$Cell.type <- colData(L4.MN)$Cell.type

L4.MN.data <- as.data.frame(colData(L4.MN)) %>% group_by(Cell.type)
L4.MN.stats <- L4.MN.data %>% summarize(
  n = n(),
  median.umi = median(total_counts),
  median.ngenes = median(num_genes_expressed)
)

L4.MN.stats

colData(L4.MN)$anatomical_class <- colData(L4.MN)$Cell.type
colData(L4.MN)$anatomical_class <- dplyr::recode(colData(L4.MN)$anatomical_class,
                                                 "AS11" = "AS",
                                                 "DA1" = "DA",
                                                 "DA2_5" = "DA",
                                                 "DA6_8" = "DA",
                                                 "DA9" = "DA",
                                                 "DB1" = "DB",
                                                 "DB2" = "DB", 
                                                 "DB3_7" = "DB",
                                                 "DD" = "DD",
                                                 "VA1" = "VA",
                                                 "VA9_11" = "VA",
                                                 "VA12" = "VA",
                                                 "VB1" = "VB",
                                                 "VB2" = "VB",
                                                 "VC4_5" = "VC",
                                                 "VC1_3_6" = "VC")

table(colData(L4.MN)$Cell.type, colData(L4.MN)$anatomical_class)

common.MN <- intersect(unique(colData(L4.MN)$Cell.type), unique(colData(MN.cds)$Cell.type))
# 14 subtypes in common. 

# I can easily now run Seurat on the cell types that are commonly identified between the 
# L4 and adult datasets.
colData(L4.MN)$stage <- "L4"
colData(MN.cds)$stage <- "YA"
MN.l4.ad <- combine_cds(list(L4.MN, MN.cds), keep_all_genes = T, cell_names_unique = T)
MN.l4.ad <- detect_genes(MN.l4.ad)
MN.l4.ad <- MN.l4.ad[rowData(MN.l4.ad)$num_cells_expressed > 5,]
table(colData(MN.l4.ad)$Cell.type, colData(MN.l4.ad)$stage)

library(Seurat)
MN.s <- CreateSeuratObject(counts = exprs(MN.l4.ad),
                           meta.data = as.data.frame(colData(MN.l4.ad)))

Idents(MN.s) <- MN.s$Cell.type

l4.ya.de <- list()
for(i in 1:length(common.MN)){
  class <- common.MN[[i]]
  sub.s <- subset(MN.s, subset = Cell.type == class)
  sub.s <- NormalizeData(sub.s)
  Idents(sub.s) <- sub.s$stage
  res <- FindMarkers(sub.s, ident.1 = "YA", ident.2 = "L4")
  res$id <- rownames(res)
  res$gene_short_name <- i2s(rownames(res), gids)
  res$neuron <- class
  res$comp <- "YA-L4"
  avg.exp <- AverageExpression(sub.s)
  avg.exp <- avg.exp[["RNA"]]
  avg.exp <- data.frame(id = rownames(avg.exp),
                        gene_short_name = i2s(rownames(avg.exp), gids),
                        avg.exp)
  res <- left_join(res, avg.exp, by = c("id", "gene_short_name"))
  l4.ya.de[[i]] <- res
}

l4.ya.de.df <- do.call("rbind", l4.ya.de)
head(l4.ya.de.df)
table(l4.ya.de.df$p_val_adj < 0.05)

# possible artifacts - pha-1 higher in YA (pha-1 rescue in lin-39 strain), 

artifacts <- c("eat-4", "lin-15B", "lin-15A", "cex-1", "dpy-20", "cho-1", "C30A5.16", "saeg-2", "unc-119", "F38B6.2", "pha-1",
               "C30F5.8", "gcy-35", "unc-54", "rol-6", "unc-53", "fip-3", "fipr-16", "C06E1.7", "unc-47", "ceh-2", "C30F8.3",
               "srg-64")

l4.ya.de.df$artifact <- ifelse(
  l4.ya.de.df$gene_short_name %in% artifacts,
  TRUE,
  FALSE
)
l4.ya.de.df.filt <- l4.ya.de.df %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) >= 1 & artifact == FALSE)

table(l4.ya.de.df.filt$neuron)
table(l4.ya.de.df.filt$avg_log2FC > 0, l4.ya.de.df.filt$neuron)
table(l4.ya.de.df.filt$avg_log2FC > 0)
# Most instances are higher in YA data

l4.ya.de.df.filt %>% arrange(desc(avg_log2FC)) %>% head(20)

l4.ya.de.df.filt %>% filter(!neuron %in% c("VC4_5", "VC1_3_6")) %>% 
  arrange(desc(avg_log2FC)) %>% head(20)

length(unique(l4.ya.de.df.filt$id))
# 14 subclasses, 1367 instances of differential expression among 783 genes.

l4.ya.de.df.filt %>% arrange(avg_log2FC) %>% head(20)

# lots of stress-related genes (heat-shock proteins, mitochondrial genes) higher in L4 than 
# in YA. Possibly due to isolation at colder temperatures in YA.

gene.l4.ya.de.table <- table(l4.ya.de.df.filt$gene_short_name)
head(gene.l4.ya.de.table)
gene.l4.ya.de.table <- as.data.frame(gene.l4.ya.de.table)
colnames(gene.l4.ya.de.table) <- c("gene_short_name", "n.cells.l4.ya.de")
head(gene.l4.ya.de.table)
rownames(gene.l4.ya.de.table) <- s2i(gene.l4.ya.de.table$gene_short_name, gids)
gene.l4.ya.de.table <- data.frame(id = rownames(gene.l4.ya.de.table),
                                  gene.l4.ya.de.table)
head(gene.l4.ya.de.table)
summary(gene.l4.ya.de.table$n.cells.l4.ya.de)

ggplot(gene.l4.ya.de.table, aes(x = n.cells.l4.ya.de)) + 
  geom_histogram(binwidth = 1) + theme_classic()

gene.l4.ya.de.table <- gene.l4.ya.de.table %>% arrange(desc(n.cells.l4.ya.de))

dir.gene <- as.data.frame(table(l4.ya.de.df.filt$avg_log2FC > 0, l4.ya.de.df.filt$gene_short_name))

head(dir.gene)

colnames(dir.gene) <- c("up_YA", "gene_short_name", "n.cells")
head(dir.gene)
dir.gene$up_YA <- as.character(dir.gene$up_YA)
dir.gene$up_YA <- dplyr::recode(dir.gene$up_YA, "FALSE" = "up_L4", "TRUE" = "up_YA")

head(dir.gene)
library(reshape2)
dir.gene <- dcast(dir.gene, gene_short_name ~ up_YA, value.var = "n.cells")
head(dir.gene)
dir.gene$n.cells.l4.ya.de <- dir.gene$up_YA + dir.gene$up_L4

head(dir.gene)

gene.l4.ya.de.table <- full_join(gene.l4.ya.de.table, dir.gene, by = c("gene_short_name", "n.cells.l4.ya.de"))
head(gene.l4.ya.de.table)

gene.l4.ya.de.table$pct_up_L4 <- 100*(gene.l4.ya.de.table$up_L4/gene.l4.ya.de.table$n.cells.l4.ya.de)

gene.l4.ya.de.table %>% filter(pct_up_L4 > 30 & pct_up_L4 < 70) %>% head(20)

l4.ya.de.df.filt %>% filter(gene_short_name == "rrn-3.1") %>% head(20)
l4.ya.de.df <- l4.ya.de.df[,c(8,6,7,1,5,2,3,4,11,10,9)]
l4.ya.de.df.filt <- l4.ya.de.df.filt[,c(8,6,7,1,5,2,3,4,11,10,9)]

l4.ya.de.df.filt %>% filter(gene_short_name %in% c("nlp-13", "nlp-15", "nlp-11"))
npep <- read.table("./neuropeptides.csv", sep = ",", header = T, stringsAsFactors = F)

l4.ya.de.df.filt %>% filter(gene_short_name %in% c(npep$gene_short_name))
write.table(l4.ya.de.df.filt, "./Kratsios_scRNAseq/Analysis/072023_l4_ya_MN_subclass_DE.csv",
            sep = ",", quote = F, row.names = F, col.names = T)

saveRDS(L4.MN, "./Kratsios_scRNAseq/Analysis/072023_updated_L4_CeNGEN_MN_subset.cds")
saveRDS(L4.updated, "./Kratsios_scRNAseq/Analysis/072023_L4_CeNGEN_neuronal_MN_updated.cds")

write.table(L4.MN.stats, "./Kratsios_scRNAseq/Analysis/072023_L4_updated_MN_per_subclass_stats.csv",
            sep = ",", row.names = F, col.names = T, quote = F)

# per neuron data
neuron.la.table <- table(l4.ya.de.df.filt$neuron)
head(neuron.la.table)
neuron.la.table <- as.data.frame(neuron.la.table)
colnames(neuron.la.table) <- c("neuron", "n.genes.DE")
head(neuron.la.table)

L4.mn.table <- as.data.frame(table(colData(L4.MN[,colData(L4.MN)$Cell.type %in% common.MN])$Cell.type))
YA.mn.table <- as.data.frame(table(colData(MN.cds[,colData(MN.cds)$Cell.type %in% common.MN])$Cell.type))
head(L4.mn.table)

neuron.la.table$n.cells.L4 <- L4.mn.table$Freq
neuron.la.table$n.cells.YA <- YA.mn.table$Freq

head(neuron.la.table)

neuron.la.table <- neuron.la.table %>% arrange(desc(n.genes.DE))

head(neuron.la.table)

dir.gene <- as.data.frame(table(l4.ya.de.df.filt$avg_log2FC > 0, l4.ya.de.df.filt$neuron))

head(dir.gene)

colnames(dir.gene) <- c("up_YA", "neuron", "n.genes")
head(dir.gene)
dir.gene$up_YA <- as.character(dir.gene$up_YA)
dir.gene$up_YA <- dplyr::recode(dir.gene$up_YA, "FALSE" = "up_L4", "TRUE" = "up_YA")

head(dir.gene)
library(reshape2)
dir.gene <- dcast(dir.gene, neuron ~ up_YA, value.var = "n.genes")
dir.gene$n.genes.DE <- dir.gene$up_L4 + dir.gene$up_YA

head(dir.gene)

neuron.la.table <- full_join(neuron.la.table, dir.gene, by = c("neuron", "n.genes.DE"))
head(neuron.la.table)

neuron.la.table$neuron <- factor(neuron.la.table$neuron, levels = neuron.la.table$neuron)

ggplot(neuron.la.table, aes(x = n.cells.L4+n.cells.YA, y = n.genes.DE)) + 
  geom_point()
head(L4.MN.stats)
colnames(L4.MN.stats) <- c("neuron", "n.cells.L4", "Median.UMI.L4", "Median.ngenes.L4")

YA.MN.data <- as.data.frame(colData(MN.cds)) %>% group_by(Cell.type)
YA.MN.stats <- YA.MN.data %>% summarize(
  n.cells.YA = n(),
  Median.UMI.YA = median(total),
  Median.ngenes.YA = median(num_genes_expressed)
)

colnames(YA.MN.stats)[1] <- "neuron"

neuron.la.table <- left_join(neuron.la.table, L4.MN.stats, by = c("neuron", "n.cells.L4"))
neuron.la.table <- left_join(neuron.la.table, YA.MN.stats, by = c("neuron", "n.cells.YA"))

head(neuron.la.table)

neuron.la.table <- neuron.la.table[,c(1,3,4,2,5,6,7,9,8,10)]
neuron.la.table

ggplot(neuron.la.table, aes(x = Median.UMI.L4, y = Median.UMI.YA)) +
  geom_point() + coord_cartesian(xlim = c(0,3000), ylim = c(0,3000)) + 
  geom_abline(intercept = 0, slope = 1, color = "darkred")

ggplot(neuron.la.table, aes(x = Median.ngenes.L4, y = Median.ngenes.YA)) +
  geom_point() + coord_cartesian(xlim = c(0,1000), ylim = c(0,1000)) + 
  geom_abline(intercept = 0, slope = 1, color = "darkred")

# Getting data from the anatomical class level (median genes/cell, UMIs/cell)

L4.orig <- readRDS("./CeNGEN/L4_scRNAseq/final_gene_model/frozen_data_sets/scRNA_publication_datasets/092320_L4_neuron_only_monocle_cds_for_paper.rds")
L4.o.mn <- L4.orig[,colData(L4.orig)$Cell.type %in% c("AS", "DA", "DA9", "DB01", "DB",
                                                      "VA", "VA12", "VB01", "VB02", "VB",
                                                      "VC", "VC_4_5", "VD_DD")]
table(colData(L4.o.mn)$Cell.type)
colData(L4.o.mn)$anatomical_class <- colData(L4.o.mn)$Cell.type
colData(L4.o.mn)$anatomical_class <- dplyr::recode(colData(L4.o.mn)$Cell.type,
                                                   "DA9" = "DA",
                                                   "DB01" = "DB",
                                                   "VA12" = "VA",
                                                   "VB01" = "VB",
                                                   "VB02" = "VB", 
                                                   "VC_4_5" = "VC")


YA.ac.data <- as.data.frame(colData(MN.cds)) %>% group_by(anatomical_class)
YA.ac.stats <- YA.ac.data %>% summarize(
  n.cells.YA = n(),
  Median.UMI.YA = median(total),
  Median.ngenes.YA = median(num_genes_expressed)
)

L4.ac.data <- as.data.frame(colData(L4.MN)) %>% group_by(anatomical_class)
L4.ac.stats <- L4.ac.data %>% summarize(
  n.cells.L4 = n(),
  Median.UMI.L4 = median(total_counts),
  Median.ngenes.L4 = median(num_genes_expressed)
)

L4o.ac.data <- as.data.frame(colData(L4.o.mn)) %>% group_by(anatomical_class)
L4o.ac.stats <- L4o.ac.data %>% summarize(
  n.cells.L4o = n(),
  Median.UMI.L4o = median(total_counts),
  Median.ngenes.L4o = median(num_genes_expressed)
)

class.stats <- full_join(YA.ac.stats, L4.ac.stats, by = "anatomical_class")
class.stats <- full_join(class.stats, L4o.ac.stats, by = "anatomical_class")

head(class.stats)
write.table(class.stats, "./Kratsios_scRNAseq/Analysis/072023_YA_L4updated_L4_original_MN_anatomical_class_stats.csv",
            sep = ",", quote = F, row.names = F, col.names = T)

# Here's a possible way to make a comparison. I can create YA annotations that match the L4 updated
# annotations, then compare the numbers and stats

colData(MN.cds)$L4_match_type <- colData(MN.cds)$Cell.type
colData(MN.cds)$L4_match_type <- dplyr::recode(colData(MN.cds)$L4_match_type, 
                                               "AS2_3" = "AS",
                                               "AS4_8" = "AS",
                                               "AS9_10" = "AS",
                                               "DD2_3" = "DD",
                                               "DD4_5" = "DD",
                                               "VA2" = "VA",
                                               "VA3_8" = "VA",
                                               "VA9_10" = "VA9_11",
                                               "VA11" = "VA9_11",
                                               "VB3" = "VB",
                                               "VB4_9" = "VB",
                                               "VB10_11" = "VB",
                                               "VD3_7" = "VD",
                                               "VD8_11" = "VD",
                                               "VD12" = "VD")

# classes to compare - AS, DD, VA, VA9_11, VB,  VD

YA.l4.match.data <- as.data.frame(colData(MN.cds)) %>% group_by(L4_match_type)
YA.l4.stats <- YA.l4.match.data %>% summarize(
  n.cells.YA = n(),
  Median.UMI.YA = median(total),
  Median.ngenes.YA = median(num_genes_expressed)
)

colnames(YA.l4.stats)[1] <- "neuron"

ya.l4.match <- full_join(YA.l4.stats, L4.MN.stats, by = "neuron")
head(ya.l4.match)

ya.l4.match <- ya.l4.match[,c(1,2,5,3,6,4,7)]

ya.l4.match %>% filter(neuron %in% c("AS", "DD", "VA", "VA9_11", "VB", "VD"))

write.table(ya.l4.match, "./Kratsios_scRNAseq/Analysis/072023_YA_L4_updated_matching_stats.csv",
            sep = ",", quote = F, row.names = F, col.names = T)

# comparing the L4-YA DE data with the within class DE data to see if genes that are enriched in 
# specific subclasses are also increasing in expression between L4 and YA.

# limiting to subclasses that have L4-YA data
common.subclass.DE <- subclass.DE %>% filter(subclass %in% common.MN)

length(intersect(common.subclass.DE$WBGeneID, l4.ya.de.df.filt$id))
# 571 genes, out of 1717 subclass-enriched genes and out of 783 L4-YA differential genes 

length(unique(common.subclass.DE$WBGeneID))
# 1717

length(unique(l4.ya.de.df.filt$id))
# 783

sub.stage.genes <- intersect(common.subclass.DE$WBGeneID, l4.ya.de.df.filt$id)

l4.ya.de.df.filt %>% 
  filter(id %in% sub.stage.genes) %>% 
  arrange(desc(avg_log2FC)) %>%
  head(20)

neuron.la.table$neuron <- factor(neuron.la.table$neuron, levels = neuron.la.table$neuron)

ggplot(neuron.la.table, aes(x = neuron, y = n.genes.DE)) + 
  geom_col(fill = "lightgrey", color = "black") +
  theme_classic() +
  theme(axis.line = element_line(color = "black"),
        axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black")) +
  labs(x = '', y = "Number of DE genes\nbetween adult and L4")

ggplot(l4.ya.de.df.filt, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point() +
  theme_classic() +
  theme(axis.line = element_line(color = "black")) +
  labs(x = "Avg log2FC", y = "-log10(p-value)")

colnames(colData(MN.cds))

colData(L4.MN)$L4_match_type <- colData(L4.MN)$Cell.type

# With L4 matching type

MN.l4.ad <- combine_cds(list(L4.MN, MN.cds), keep_all_genes = T, cell_names_unique = T)
MN.l4.ad <- detect_genes(MN.l4.ad)
MN.l4.ad <- MN.l4.ad[rowData(MN.l4.ad)$num_cells_expressed > 5,]
table(colData(MN.l4.ad)$Cell.type, colData(MN.l4.ad)$stage)
table(colData(MN.l4.ad)$L4_match_type, colData(MN.l4.ad)$stage)

library(Seurat)
MN.s <- CreateSeuratObject(counts = exprs(MN.l4.ad),
                           meta.data = as.data.frame(colData(MN.l4.ad)))

Idents(MN.s) <- MN.s$L4_match_type
matching.common <- sort(unique(colData(MN.l4.ad)$L4_match_type))

l4.match.de <- list()
for(i in 1:length(matching.common)){
  class <- matching.common[[i]]
  sub.s <- subset(MN.s, subset = L4_match_type == class)
  sub.s <- NormalizeData(sub.s)
  Idents(sub.s) <- sub.s$stage
  res <- FindMarkers(sub.s, ident.1 = "YA", ident.2 = "L4")
  res$id <- rownames(res)
  res$gene_short_name <- i2s(rownames(res), gids)
  res$neuron <- class
  res$comp <- "YA-L4"
  avg.exp <- AverageExpression(sub.s)
  avg.exp <- avg.exp[["RNA"]]
  avg.exp <- data.frame(id = rownames(avg.exp),
                        gene_short_name = i2s(rownames(avg.exp), gids),
                        avg.exp)
  res <- left_join(res, avg.exp, by = c("id", "gene_short_name"))
  l4.match.de[[i]] <- res
}

l4.match.de.df <- do.call("rbind", l4.match.de)
head(l4.match.de.df)
table(l4.match.de.df$p_val_adj < 0.05)

# possible artifacts - pha-1 higher in YA (pha-1 rescue in lin-39 strain), 

artifacts <- c("eat-4", "lin-15B", "lin-15A", "cex-1", "dpy-20", "cho-1", "C30A5.16", "saeg-2", "unc-119", "F38B6.2", "pha-1",
               "C30F5.8", "gcy-35", "unc-54", "rol-6", "unc-53", "fip-3", "fipr-16", "C06E1.7", "unc-47", "ceh-2", "C30F8.3",
               "srg-64")

l4.match.de.df$artifact <- ifelse(
  l4.match.de.df$gene_short_name %in% artifacts,
  TRUE,
  FALSE
)
l4.match.de.df.filt <- l4.match.de.df %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) >= 1 & artifact == FALSE)

table(l4.match.de.df.filt$neuron)
table(l4.match.de.df.filt$avg_log2FC > 0, l4.match.de.df.filt$neuron)
table(l4.match.de.df.filt$avg_log2FC > 0)
# Most instances are higher in YA data

l4.match.de.df.filt %>% arrange(desc(avg_log2FC)) %>% head(20)

l4.match.de.df.filt %>% filter(!neuron %in% c("VC4_5", "VC1_3_6")) %>% 
  arrange(desc(avg_log2FC)) %>% head(20)

length(unique(l4.match.de.df.filt$id))
# 14 subclasses, 2247 instances of differential expression among 1037 genes.

l4.match.de.df.filt %>% arrange(avg_log2FC) %>% head(20)

# lots of stress-related genes (heat-shock proteins, mitochondrial genes) higher in L4 than 
# in YA. Possibly due to isolation at colder temperatures in YA.

gene.l4.match.de.table <- table(l4.match.de.df.filt$gene_short_name)
head(gene.l4.match.de.table)
gene.l4.match.de.table <- as.data.frame(gene.l4.match.de.table)
colnames(gene.l4.match.de.table) <- c("gene_short_name", "n.cells.l4.match.de")
head(gene.l4.match.de.table)
rownames(gene.l4.match.de.table) <- s2i(gene.l4.match.de.table$gene_short_name, gids)
gene.l4.match.de.table <- data.frame(id = rownames(gene.l4.match.de.table),
                                     gene.l4.match.de.table)
head(gene.l4.match.de.table)
summary(gene.l4.match.de.table$n.cells.l4.match.de)

ggplot(gene.l4.match.de.table, aes(x = n.cells.l4.match.de)) + 
  geom_histogram(binwidth = 1) + theme_classic()

gene.l4.match.de.table <- gene.l4.match.de.table %>% arrange(desc(n.cells.l4.match.de))

dir.gene <- as.data.frame(table(l4.match.de.df.filt$avg_log2FC > 0, l4.match.de.df.filt$gene_short_name))

head(dir.gene)

colnames(dir.gene) <- c("up_YA", "gene_short_name", "n.cells")
head(dir.gene)
dir.gene$up_YA <- as.character(dir.gene$up_YA)
dir.gene$up_YA <- dplyr::recode(dir.gene$up_YA, "FALSE" = "up_L4", "TRUE" = "up_YA")

head(dir.gene)
library(reshape2)
dir.gene <- dcast(dir.gene, gene_short_name ~ up_YA, value.var = "n.cells")
head(dir.gene)
dir.gene$n.cells.l4.match.de <- dir.gene$up_YA + dir.gene$up_L4

head(dir.gene)

gene.l4.match.de.table <- full_join(gene.l4.match.de.table, dir.gene, by = c("gene_short_name", "n.cells.l4.match.de"))
head(gene.l4.match.de.table)

gene.l4.match.de.table$pct_up_L4 <- 100*(gene.l4.match.de.table$up_L4/gene.l4.match.de.table$n.cells.l4.match.de)

gene.l4.match.de.table %>% filter(pct_up_L4 > 30 & pct_up_L4 < 70) %>% head(20)

l4.match.de.df.filt %>% filter(gene_short_name %in% c("nlp-13", "nlp-15", "nlp-11"))
npep <- read.table("./neuropeptides.csv", sep = ",", header = T, stringsAsFactors = F)

l4.match.de.df.filt %>% filter(gene_short_name %in% c(npep$gene_short_name))
write.table(l4.match.de.df.filt, "./Kratsios_scRNAseq/Analysis/072023_l4_ya_MN_l4_matching_subclass_DE.csv",
            sep = ",", quote = F, row.names = F, col.names = T)

saveRDS(L4.MN, "./Kratsios_scRNAseq/Analysis/072023_updated_L4_CeNGEN_MN_subset.cds")
saveRDS(L4.updated, "./Kratsios_scRNAseq/Analysis/072023_L4_CeNGEN_neuronal_MN_updated.cds")

write.table(L4.MN.stats, "./Kratsios_scRNAseq/Analysis/072023_L4_updated_MN_per_subclass_stats.csv",
            sep = ",", row.names = F, col.names = T, quote = F)

# per neuron data
neuron.la.match.table <- table(l4.match.de.df.filt$neuron)
neuron.la.match.table <- as.data.frame(neuron.la.match.table)
colnames(neuron.la.match.table) <- c("neuron", "n.genes.DE")
head(neuron.la.match.table)

L4.mn.table <- as.data.frame(table(colData(L4.MN)$L4_match_type))
YA.mn.table <- as.data.frame(table(colData(MN.cds)$L4_match_type))
head(L4.mn.table)

neuron.la.match.table$n.cells.L4 <- L4.mn.table$Freq
neuron.la.match.table$n.cells.YA <- YA.mn.table$Freq

head(neuron.la.match.table)

neuron.la.match.table <- neuron.la.match.table %>% arrange(desc(n.genes.DE))

head(neuron.la.match.table)

dir.gene <- as.data.frame(table(l4.match.de.df.filt$avg_log2FC > 0, l4.match.de.df.filt$neuron))

head(dir.gene)

colnames(dir.gene) <- c("up_YA", "neuron", "n.genes")
head(dir.gene)
dir.gene$up_YA <- as.character(dir.gene$up_YA)
dir.gene$up_YA <- dplyr::recode(dir.gene$up_YA, "FALSE" = "up_L4", "TRUE" = "up_YA")

head(dir.gene)
library(reshape2)
dir.gene <- dcast(dir.gene, neuron ~ up_YA, value.var = "n.genes")
dir.gene$n.genes.DE <- dir.gene$up_L4 + dir.gene$up_YA

head(dir.gene)

neuron.la.match.table <- full_join(neuron.la.match.table, dir.gene, by = c("neuron", "n.genes.DE"))
head(neuron.la.match.table)

neuron.la.match.table$neuron <- factor(neuron.la.match.table$neuron, levels = neuron.la.match.table$neuron)

ggplot(neuron.la.match.table, aes(x = n.cells.L4+n.cells.YA, y = n.genes.DE)) + 
  geom_point()
head(L4.MN.stats)
colnames(L4.MN.stats) <- c("neuron", "n.cells.L4", "Median.UMI.L4", "Median.ngenes.L4")

YA.MN.data <- as.data.frame(colData(MN.cds)) %>% group_by(L4_match_type)
YA.MN.stats <- YA.MN.data %>% summarize(
  n.cells.YA = n(),
  Median.UMI.YA = median(total),
  Median.ngenes.YA = median(num_genes_expressed)
)

colnames(YA.MN.stats)[1] <- "neuron"

neuron.la.match.table <- left_join(neuron.la.match.table, L4.MN.stats, by = c("neuron", "n.cells.L4"))
neuron.la.match.table <- left_join(neuron.la.match.table, YA.MN.stats, by = c("neuron", "n.cells.YA"))

head(neuron.la.match.table)

neuron.la.match.table <- neuron.la.match.table[,c(1,3,4,2,5,6,7,9,8,10)]
neuron.la.match.table

# Saving all files from these analyses.
saveRDS(neuron.la.match.table, "./Kratsios_scRNAseq/Analysis/072523_L4_matching_L4vYA_DE_neuron_level_table.rds")
saveRDS(neuron.la.table, "./Kratsios_scRNAseq/Analysis/072523_L4vYA_14_subclasses_DE_neuron_level_table.rds")
saveRDS(gene.l4.ya.de.table, "./Kratsios_scRNAseq/Analysis/072523_L4vYA_14_subclasses_DE_gene_level_table.rds")
saveRDS(gene.l4.match.de.table, "./Kratsios_scRNAseq/Analysis/072523_L4_matching_L4vYA_DE_gene_level_table.rds")
saveRDS(l4.match.de.df.filt, "./Kratsios_scRNAseq/Analysis/072523_L4_matching_L4vYA_DE_results_filtered.rds")
saveRDS(l4.ya.de.df.filt, "./Kratsios_scRNAseq/Analysis/072523_L4vYA_14_subclasses_DE_results_filtered.rds")

# Generating csv file of thresholded data with gene symbols

library(wbData)
gids <- wb_load_gene_ids("WS273")

YA.t2 <- readRDS("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/Analysis/060823_YA_TPM_t2_expression_matrix.rds")
head(rownames(YA.t2))
colnames(YA.t2)

MN.t2 <- YA.t2[,c("AS1_3", "AS4_8", "AS9_10", "AS11", "DA1", "DA2_5", "DA6_8", "DA9", "DB1",
                  "DB2", "DB3_7", "DD2_3", "DD4_5", "VA1", "VA2", "VA3_8", "VA9_10", "VA11",
                  "VA12", "VB1", "VB2", "VB3", "VB4_9", "VB10_11", "VC1_3_6", "VC4_5", "VD3_7",
                  "VD8_11", "VD12")]

rownames(MN.t2) <- i2s(rownames(MN.t2), gids)
head(rownames(MN.t2))
dim(MN.t2)
# 10710 genes in 29 cell types

MN.t2 <- as.matrix(MN.t2)
MN.t2[1:10,1:10]
MN.t2[500:510, 1:10]
colnames(MN.t2)[1] <- "AS2_3"
saveRDS(MN.t2, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/Analysis/072723_MN_threshold2_expression_gene_symbols.rds")
write.table(MN.t2, "~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/Files_for_website/072723_MN_subclass_thresholded_expression.csv",
            sep = ",", quote = F, row.names = T, col.names = NA)

# Checking if Seurat object really uses gene symbols.
library(Seurat)
MN.s <- readRDS("~/Dropbox (VU Basic Sciences)/Miller lab/10X Genomics/Kratsios_scRNAseq/Files_for_website/072523_wt_MN_VNC_Seurat_updated_gene_symbols.rds")
MN.s
MN.s@assays$RNA[1:5,1:5]

