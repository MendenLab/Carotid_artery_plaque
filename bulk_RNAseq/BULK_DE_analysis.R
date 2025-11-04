#################################################
## 1. Data Import, Quality Control & Filtering
#################################################

# Load required libraries
library(data.table)
library(ggplot2)
library(magrittr)
library(scales)   
library(DESeq2)
library(BiocParallel)
library(ggrepel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read data files
sa <- fread("sample_maifest.csv")
counts <- readRDS("Counts_batch_corrected.rds")


###########################################
## Differential Expression Analysis
###########################################

# Ensure sample annotation matches filtered counts
sa <- sa[Data_ID %in% colnames(counts)]
sa <- sa[match(colnames(counts), sa$Data_ID), ]

sa[, status := as.factor(status)] #plaques phenotype as factor
sa$status <- relevel(sa$status, ref = "Control")
sa[, Sex := as.factor(Sex)] 
sa[, Patient := as.factor(Patient)] 
sa[, Plaque_Phenotype := as.factor(Plaque_Phenotype)]

# Create a DESeqDataSet with design for unpaired analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sa,
                              design = ~ Age + Sex + status)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# Extract DE results for comparison (Diseased vs Control)
res_status <- results(dds, alpha = 0.05, pAdjustMethod = "BH")
res_status <- as.data.table(res_status, keep.rownames = TRUE)
res_status[, comparison := "Unpaired \n Plaques vs Controls"]


## (B) Unpaired DE Analysis- Stables versus unstables 
design(dds) <- formula(~ Age + Sex + Plaque_Phenotype)
dds <- DESeq(dds)
res_plq_phe <- results(dds,  
                       alpha = 0.05, 
                       contrast =  c("Plaque_Phenotype", "Plaque_stable","Plaque_unstable"), 
                       pAdjustMethod =  "BH") 

res_plq_phe <- as.data.table(res_plq_phe, keep.rownames = TRUE )
res_plq_phe[, comparison := "Unpaired \nStable vs Unstable"]


## (C) Paired DE Analysis Using Patient ID
# Select paired samples (patients with duplicate entries and with a control plaque)
sa_paired <- sa[Patient %in% sa[duplicated(Patient), Patient]]
table(sa_paired$status)

# Create a DESeq2 dataset for paired analysis
counts_paired <- counts[, colnames(counts) %in% sa_paired$Data_ID]
sa_paired <- sa_paired[match(colnames(counts_paired), sa_paired$Data_ID), ]
dds_paired <- DESeqDataSetFromMatrix(countData = counts_paired,
                                     colData = sa_paired,
                                     design = ~ Patient + status)
#BiocParallel::register(MulticoreParam(4))
#dds_paired <- DESeq(dds_paired, parallel = TRUE)
dds_paired <- DESeq(dds_paired, parallel = FALSE)

# Extract DE results for the paired design (status comparison)
res_paired_status <- results(dds_paired, alpha = 0.05, pAdjustMethod = "BH")
res_paired_status <- as.data.table(res_paired_status, keep.rownames = TRUE)
res_paired_status[, comparison := "Paired \n Plaques vs Controls"]

#combining all the result in one table
res_all <- rbind(res_plq_phe, res_status, res_paired_status)
res_all[, significant := ifelse(padj < 0.05 & abs(log2FoldChange) >= 1, "TRUE", "FALSE")]

###########################################
## 4. Visualization
###########################################
res_all[padj < 0.05 & (log2FoldChange) >= 1, significant_2 := "Upregulated"]
res_all[padj < 0.05 & (log2FoldChange) <= -1 , significant_2 := "Downregulated"]
res_all[is.na(significant_2), significant_2 := "n.s"]
table(res_all$comparison, res_all$significant_2)
table(res_all$comparison, res_all$significant)

## Save source dataframes for plots
library(openxlsx)

# Create directory if it doesn't exist
output_dir = '../../source_data_raw/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)}

source_df = res_all[,c("rn","log2FoldChange","padj","comparison","significant_2")]
colnames(source_df) = c("Gene","log2FoldChange","padj","comparison","significance")
output_file <- file.path(output_dir, "Fig2A_B_C.xlsx")
write.xlsx(source_df,output_file)


mycol <-  c("n.s" = "grey", 
            "Upregulated" = 	"#e31a1c",
            "Downregulated" = "#2171b5")


xen=read.csv("Codes_BulkRNAseq/XENIUM panels_060923.csv",sep=';')
xen_gene <- xen[xen[[1]] == "PANEL 1" | xen$LISTS == "CORE", "Gene"]
res_all_fig=res_all

de_plot<-ggplot(res_all[!is.na(padj) , ], 
                aes(log2FoldChange, -log10(padj), 
                    color = significant_2)) + 
  geom_point() + 
  facet_wrap(~comparison) + 
  scale_color_manual(values = mycol) + 
  labs(col = "Significance") +  
  theme_classic() +  
  geom_label_repel(data = res_all_fig[comparison == "Unpaired \n Plaques vs Controls" & rn %in% xen_gene & abs(log2FoldChange) > 2 & padj < 0.00001, ], aes(label = rn), max.overlaps = 20,
                   show.legend = FALSE) +
  geom_label_repel(data = res_all_fig[comparison == "Unpaired \nStable vs Unstable"& significant == TRUE , ], aes(label = rn), max.overlaps = 20, 
                   show.legend = FALSE) +
  geom_label_repel(data = res_all_fig[comparison == "Paired \n Plaques vs Controls"  & rn %in% xen_gene & abs(log2FoldChange) > 3 & padj < 0.00001 , ], aes(label = rn), max.overlaps = 20, 
                   show.legend = FALSE) +
  theme( legend.key.size = unit(0.9, "cm"), 
         strip.text.x = element_text(size = 20, face = "bold"), 
         axis.text = element_text(size = 20, color = "black"), 
         axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size =13), 
         legend.title = element_text(size =20, face = "bold")) + 
  guides(color = guide_legend(override.aes = list(size = 5)))


## Saving plot
output_dir <- "/root/capsule/results/figure_plots/Fig2/Bulk_DGE_results"
output_file <- file.path(output_dir, "bulk_DGE_results.png")

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save
ggsave(
  filename = output_file,
  plot = de_plot,
  width = 14,
  height = 5,
  units='in',
  dpi = 300
)

print('a')

output_dir = '../../Supplementary_tables/'
output_file <- file.path(output_dir, "Supplementary_Table6.xlsx")
write.xlsx(counts,output_file)
