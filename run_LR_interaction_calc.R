#################################################################################
## Modified Date : 190216, by YR
## Original : by Woosung Chung
## Goal : anlaysis of R-L pairs (Tumor, Normal, each patient...)
#################################################################################
library(dplyr)
library(Matrix)
library(rgl)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(circlize)
library(tibble)
library(igraph)

##### Control ############################################################################################
directory <- "D:/HYR/10x/HYR_result/cellranger_2.1/CRC/Final_SMC/After_3rdQC_meta/Interaction"
cell_type_order <- c(tep.order, nep.order, stromal.order, myeloids.order, tcell.order, bcell.order)
cell_type_order

## parameters ##
#ORIGIN = "Tumor" #or Normal
PerSample = "TRUE" #or FALSE of Matched

## Tumor / Normal data
Tumor.all.lnTPMdata = "After_3rdQC_SMC_Tumor_lnTPM.rds"
Normal.all.lnTPMdata = "After_3rdQC_SMC_Normal_lnTPM.rds"
metadata = "Merged_5_major_types_meta_190721_Tumor_EP_Corr_CMS_changed_and_MF_merged.txt"
RL_pair_path = "Ramilowski_receptor_ligand_pairs.txt"
clinical_info_path = "CRC_SMC_Clinical_NYedited_and_YR.txt"

setwd(directory)
#######################################################################
## Loading cell info
total_cell_info <- read.table(metadata, sep="\t", header=T, stringsAsFactors = F)
str(total_cell_info)
table(total_cell_info$cell_type)
table(total_cell_info$cell_subtype)

str(total_cell_info)

total_cell_info <- data.frame(lapply(total_cell_info, as.character), stringsAsFactors = FALSE)
rownames(total_cell_info) <- total_cell_info$cell_barcode
str(total_cell_info)
total_cell_info$origin <- gsub("N", "Normal", total_cell_info$origin)
total_cell_info$origin <- gsub("T", "Tumor", total_cell_info$origin)

## load the receptor-ligand pairs (Ramilowski et al)
RL_pairs <- read.table(RL_pair_path, sep="\t", header=T, stringsAsFactors = F)

## Set the data (TOTAL)
clinic_info <- read.table(clinical_info_path, sep="\t", header=T, stringsAsFactors = F) 

################################################################
## 2. color code 
source(file = "color_code_new_final.R")

alpha50_grey <- adjustcolor("grey", alpha.f = 0.5)
total_cols <- c(heatmap.v3)
################################################################

## load the receptor-ligand pairs (Ramilowski et al)
RL_pairs <- read.table(RL_pair_path, sep="\t", header=T, stringsAsFactors = F)
## Set the data (TOTAL)
clinic_info <- read.table(clinical_info_path, sep="\t", header=T) 

tumor.samples <- unique(total_cell_info[total_cell_info$origin == "Tumor",]$patient)
tumor.samples; length(tumor.samples)

directory <- "Total_Tumor"
setwd(directory)
calc_RL_count_group(ORIGIN = "Tumor", dat_path = Tumor.all.lnTPMdata,
                    SAMPLES = tumor.samples, RL_pairs = RL_pairs, clinical_info = clinical_info, 
                    edge_trim_cutoff = 1000000, edge_label = "10M",
                    resize_node = 0.05, resize_width = 1000000, total_cell_info = total_cell_info,
                    Sample_label = "Total_Tumor")


