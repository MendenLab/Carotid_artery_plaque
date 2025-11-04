#package installation 
pkgs <- c("org.Hs.eg.db", "snowfall")
new.pkg <- pkgs[!(pkgs %in% utils::installed.packages()[, "Package"])]
if (length(new.pkg)) BiocManager::install(new.pkg)
if(!requireNamespace("BayesPrism")) devtools::install_github("Danko-Lab/BayesPrism/BayesPrism", force = TRUE)
if(!requireNamespace("InstaPrism")) devtools::install_github("humengying0907/InstaPrism", force = TRUE)

library(InstaPrism)
library(scran)
library(data.table)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# parameters
organism <- "human"
ensemblSpecies <- "hsapiens_gene_ensembl"
celldex.method <- "HumanPrimaryCellAtlasData"
rowNamesOfCounts <- "external_gene_name"

mtCutOff <- 15 # mitocondrial percent
nTopVarGenes <- 10
project <- "decon"

################################################################################
################################ Deconvolution #################################
################################################################################
# bulk input
sa <- fread("sample_maifest.csv")

counts <- readRDS("Counts_batch_corrected.rds")

# Ref
sce <- readRDS("GSE159677_filtered_batch_corrected_clustered_annotated.rds")
scExpr <- SummarizedExperiment::assay(sce)

bulkChange2Ensembl <- TRUE
if (bulkChange2Ensembl == TRUE) {
  # to ensembl gene ids
  library(org.Hs.eg.db)
  ensmbl_ids <- mapIds(org.Hs.eg.db, keys = rownames(counts),
                       column="ENSEMBL", keytype="SYMBOL")
  if(identical(rownames(counts), names(ensmbl_ids))) {
    rownames(counts) <- ensmbl_ids
  }
  mean(is.na(rownames(counts)))
  dim(counts); bulk_expr <- counts[!is.na(rownames(counts)), ]; dim(counts)
} else if (sceChange2Symbol == TRUE) {
  rownames(scExpr) <- SummarizedExperiment::rowData(sce)$gene_name
}

cell_type_labels <- SummarizedExperiment::colData(sce)$final_high_level_celltypes



cell.state.labels <- get_subcluster(scExpr = scExpr,
                                    cell_type_labels = cell_type_labels,
                                    subcluster_method = 'scran',
                                    min.subcluster.size = 200)


ref <- refPrepare(sc_Expr = scExpr, cell.type.labels = cell_type_labels,
                  cell.state.labels = cell.state.labels)

deconv_res <- InstaPrism::InstaPrism(bulk_Expr = counts, refPhi_cs = ref)

# Obtain res
estimated_frac <- t(deconv_res@Post.ini.ct@theta)

################################################################################
################################ Visualization #################################
################################################################################
estimated_frac <- as.data.table(estimated_frac, keep.rownames = "Data_ID")
estimated_frac <- merge(estimated_frac, sa[, .(Data_ID, status)], by= "Data_ID")


col_cells <- c(
  'B-cell' = '#1f77b4',
  'Dendritic-cell' = '#ff7f0e',
  'Endothelial' = '#279e68',
  'FB' = '#d62728',
  'Macrophage' = '#aa40fc',
  'Mast-cell' = '#8c564b',
  'NK-cell' = '#e377c2',
  'Plasma-cell' = '#b5bd61',
  'T-cell' = '#17becf',
  'VSMC' = '#aec7e8')

df_long <- estimated_frac %>%
  pivot_longer(
    cols = c(Endothelial, VSMC, FB, Macrophage, `NK-cell`,
             `Dendritic-cell`, `T-cell`, `Mast-cell`, `B-cell`, `Plasma-cell`),
    names_to = "CellType",
    values_to = "Fraction"
  )


cell_order <- c(
  "VSMC", "T-cell", "Macrophage", "Endothelial", "NK-cell",
  "Dendritic-cell", "Mast-cell", "Plasma-cell", "B-cell", "FB"
)

# Factor CellType with correct order
df_long$CellType <- factor(df_long$CellType, levels = cell_order)

# Compute mean fraction per group
df_summary <- df_long %>%
  group_by(status, CellType) %>%
  summarise(mean_fraction = mean(Fraction), .groups = "drop")


## Save source dataframes for plots
library(openxlsx)

# Create directory if it doesn't exist
output_dir = '../../source_data/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)}

## Control dataframe
control_df=df_long[df_long$status == "Control",]
colnames(control_df)=c('original_sample','condition','high_level_celltype','cell_fraction')
output_file <- file.path(output_dir, "Fig2H.xlsx")
write.xlsx(control_df,output_file)

## Plaque dataframe
plaque_df=df_long[df_long$status == "Plaque",]
colnames(plaque_df)=c('original_sample','condition','high_level_celltype','cell_fraction')
output_file <- file.path(output_dir, "Fig2I.xlsx")
write.xlsx(plaque_df,output_file)




#making graph for control
control_plot=ggplot(df_long[df_long$status == "Control",], aes(x = CellType, y = Fraction, fill = CellType)) +
                    geom_col(data = df_summary[df_summary$status == "Control",], aes(y = mean_fraction), width = 0.7) +
                    geom_jitter(shape = 21, color = "black", width = 0.2, size = 1.5, alpha = 0.7) +
                    # facet_wrap(~status, nrow = 2) +
                    scale_fill_manual(values = col_cells) +
                    ylab("Cell fraction per\n sample (deconv.)") +
                    xlab("") +
                    theme_classic() +
                    theme(
                      axis.text.x = element_text(angle = 45, hjust = 1, size = 18,color = "black"),
                      # strip.text = element_text(size = 14),
                      axis.text.y = element_text(size = 14, color = "black"), 
                      axis.title = element_text(size = 16, face = "bold"),
                      legend.position = "none")
dev.off()



## Saving plot
output_dir <- "/root/capsule/results/figure_plots/Fig2/cell_fraction_deconvolution_controls"
output_file <- file.path(output_dir, "cell_fractions.png")

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save
ggsave(
  filename = output_file,
  plot = control_plot,
  dpi = 300
)



#graph for plaque
plaque_plot=ggplot(df_long[df_long$status == "Plaque",], aes(x = CellType, y = Fraction, fill = CellType)) +
                  geom_col(data = df_summary[df_summary$status == "Plaque",], aes(y = mean_fraction), width = 0.7) +
                  geom_jitter(shape = 21, color = "black", width = 0.2, size = 1.5, alpha = 0.7) +
                  # facet_wrap(~status, nrow = 2) +
                  scale_fill_manual(values = col_cells) +
                  ylab("Cell fraction per\n sample (deconv.)") +
                  xlab("") +
                  theme_classic() +
                  theme(
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 18,color = "black"),
                    # strip.text = element_text(size = 14),
                    axis.text.y = element_text(size = 14, color = "black"), 
                    axis.title = element_text(size = 16, face = "bold"),
                    legend.position = "none")

## Saving plot
output_dir <- "/root/capsule/results/figure_plots/Fig2/cell_fraction_deconvolution_plaques"
output_file <- file.path(output_dir, "cell_fractions.png")

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save
ggsave(
  filename = output_file,
  plot = plaque_plot,
  dpi = 300
)



# df_summary <- df_long %>%
#   group_by(status, CellType) %>%
#   summarize(mean_fraction = mean(Fraction), .groups = "drop")
# 
# cell_con <- ggplot(df_summary, aes(x = status, y = mean_fraction, fill = CellType)) +
#   geom_bar(stat = "identity", position = "fill") +
#   scale_fill_manual(values =col_cells) +
#   ylab("Fraction") +
#   xlab("Status") +
#   theme_classic() + 
#   theme(legend.key.size = unit(0.9, "cm"), 
#         strip.text.x = element_text(size = 12),
#         axis.text = element_text(size = 13, color = "black"), 
#         axis.title = element_text(size = 15, face = "bold"), 
#         legend.text = element_text(size =13), 
#         legend.title = element_text(size =15, face = "bold")) 






