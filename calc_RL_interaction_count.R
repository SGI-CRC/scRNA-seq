####### functions ##########
calc_RL_count_group <- function(ORIGIN, dat_path, SAMPLES, Sample_label, RL_pairs, clinical_info, edge_trim_cutoff, edge_label,resize_node, resize_width,
                                total_cell_info, perSamples = NULL){
  
  if(ORIGIN == "Tumor"){
    load(dat_path)
    dat <- smc.t.lnTPM
    dim(dat)  
  }else if(ORIGIN == "Normal"){
    load(dat_path)
    dat <- smc.n.lnTPM
    dim(dat)
  }
  
  ##### 190606 add : using only RL pairs genes ######
  ligand.name <- as.character(RL_pairs$Ligand.ApprovedSymbol)
  receptor.name <- as.character(RL_pairs$Receptor.ApprovedSymbol)
  lr.name <- union(ligand.name, receptor.name)
  lr.name.inter <- intersect(rownames(dat), lr.name)
  dat <- dat[lr.name.inter,]
  dat <- as.matrix(dat)
  ##
  list_WS <- calc_RL_count_per_each(dat, SAMPLES, RL_pairs, total_cell_info, sample_index = Sample_label)
  save(list_WS, file = paste0(Sample_label,"_list_of_Edge_Node_WS.rds"))
}


calc_RL_count_per_each <- function(dat, sample, RL_pairs, total_cell_info, sample_index){
  ##
  
  ##
  ## Data filtering
  cell_barcode <- filter(total_cell_info, patient %in% sample)$cell_barcode
  length(cell_barcode)
  dat2 <- dat[,cell_barcode]; dim(dat2)
  
  #colnames(dat) <- paste0(sample_index, colnames(dat))
  chk_cell_num <- dim(filter(total_cell_info, patient %in% sample)) # check cell number 
  chk_cell_num
  
  df <- merge(total_cell_info, t(dat2), by=0, sort=FALSE); dim(df) #merge datasets
  df <- column_to_rownames(df, "Row.names")
  #df
  dim(df)
  df[1:5,1:10]
  cat(paste0("Check out the cell number of ", sample_index, ": ", chk_cell_num[1] == dim(df)[1]),"\n")
  
  ## rm all zero gene
  all_zero_gene <- grep("TRUE", apply(df[,-c(1:5)], 2, sum) == 0, value=TRUE)
  length(all_zero_gene)
  all_zero_gene
  cat("Remove all zero genes\n")
  
  df_rm_all_zero <- select(df, -one_of(names(all_zero_gene)))
  
  ## remove Unspecified plasma & SMC20-specific macrophage
  df_rm_all_zero <- subset(df_rm_all_zero, !(df_rm_all_zero$cell_subtype %in% c("Unspecified Plasma", "Anti-inflammatory-SMC20")))
  
  df_rm_all_zero$cell_subtype <- factor(df_rm_all_zero$cell_subtype, levels=intersect(names(total_cols),names(table(df_rm_all_zero$cell_subtype))))
  
  ## total receptor-ligand genes
  RL_genes <- unique(c(as.character(RL_pairs$Ligand.ApprovedSymbol), as.character(RL_pairs$Receptor.ApprovedSymbol)))
  length(RL_genes)
  head(RL_genes) # check R-L genes
  
  ## overlaped with expressed genes
  overlap_genes <- intersect(colnames(df_rm_all_zero), RL_genes)
  length(overlap_genes)
  
  
  ## extract only total receptor-ligand genes from dataset (by cell_subtype)
  trimmed_df2 <- df_rm_all_zero[, c("cell_subtype",overlap_genes)]
  trimmed_df2[1:5,1:5]
  print(table(trimmed_df2$cell_subtype))
  write.table(table(trimmed_df2$cell_subtype), paste0("TOTAL_",sample_index, "_Ramilowski_Receptor_Ligand_Analysis_cell_subtype_cell_number.txt"), sep="\t", quote=F, row.names=F)
  
  ## Calculate count of each gene over cutoff (default = 0)
  cutoff = 0
  
  ## check the cell subtype in Normal / Tumor each.
  trimmed_df_count2 <- data.frame(cell_subtype = intersect(names(total_cols),names(table(df_rm_all_zero$cell_subtype))))
  
  ## per each genes, check expressed cell # (for each cell types)
  for(w in 1:length(overlap_genes)){
    count_tmp <- data.frame(table(trimmed_df2[trimmed_df2[,overlap_genes[w]] > cutoff,]$cell_subtype))
    colnames(count_tmp) <- c("cell_subtype", overlap_genes[w])
    
    trimmed_df_count2 <- merge(trimmed_df_count2, count_tmp, by="cell_subtype", all.x=TRUE)
  }
  rownames(trimmed_df_count2) <- trimmed_df_count2[,"cell_subtype"]
  #trimmed_df_count2[is.na(trimmed_df_count2)] <- 0
  
  
  ## Set nodes & links
  nodes2 <- data.frame(id=sprintf("s%02d",1:length(table(df_rm_all_zero$cell_subtype))),table(df_rm_all_zero$cell_subtype)); 
  colnames(nodes2) <- c("id","cell_subtype","Num")
  
  links2 <- cbind(expand.grid(as.character(nodes2$id), as.character(nodes2$id)), 
                  expand.grid(names(table(df_rm_all_zero$cell_subtype)), names(table(df_rm_all_zero$cell_subtype))))
  colnames(links2) <- c("from","to","ligand","receptor")  ## replacement or permutation O A-A O
  links2$num_of_pairs <- rep(0, nrow(links2))
  str(links2)
  
  
  ## remove self interaction
  links2 <- links2[links2$from != links2$to,]
  
  ## YR add : remove intersections per each Tumor (ex, CMS1 - CMS2...)
  cms = c("CMS1", "CMS2", "CMS3", "CMS4")
  link2.new <- c()
  
  ##converting factor to chr
  links2$from <- as.character(links2$from); links2$to <- as.character(links2$to)
  links2$ligand <- as.character(links2$ligand); links2$receptor <- as.character(links2$receptor)
  
  for(i in 1:nrow(links2)){
    if(links2$ligand[i] %in% cms & links2$receptor[i] %in% cms){
      cat(links2$ligand[i], "_", links2$receptor[i], "_rm\n")
    }else{
      link2.new <- rbind(link2.new, c(as.character(links2$from[i]), as.character(links2$to[i]), links2$ligand[i], links2$receptor[i],
                                      links2$num_of_pairs[i]))
    }
  }
  colnames(link2.new) <- colnames(links2)
  ### from      to 
  ### ligand    receptor
  head(link2.new)
  
  link2.total <- as.data.frame(links2, stringsAsFactors = F) # save as link2.total
  links2 <- as.data.frame(link2.new, stringsAsFactors = F)
  
  ## convert as df
  links2 <- as.data.frame(links2, stringsAsFactors = F)
  
  ## Calculate the number of pairs of each cell type pair
  cc_log <- c() # no.of cells in each ligand-receptor pairs
  
  links2$ligand <- as.character(links2$ligand)
  links2$receptor <- as.character(links2$receptor)
  for(i in 1:nrow(links2)){
    ##links2 : per each cell type interection (ligand - receptor)
    ## ex
    ## s02 s01    Mature Enterocytes type 1     Stem            0
    ligand.genes <- trimmed_df_count2[trimmed_df_count2$cell_subtype == links2$ligand[i],]
    receptor.genes <- trimmed_df_count2[trimmed_df_count2$cell_subtype == links2$receptor[i],]
    #
    ligand.genes$cell_subtype <- NULL; receptor.genes$cell_subtype <- NULL
    l_r.df <- data.frame(ligand = t(ligand.genes), receptor = t(receptor.genes))
    colnames(l_r.df) <- c("ligand_NOcell", "receptor_NOcell")
    #l_r.df.new <- subset(l_r.df, rownames(l_r.df) != "cell_subtype")
    genes.new <- data.frame(ligand.genes = rownames(l_r.df), receptor.genes = rownames(l_r.df))
    #l_r.df.new$genes <- rownames(l_r.df.new)
    expand <- as.matrix(expand.grid(l_r.df, stringsAsFactors = F))
    expand.genes <- as.matrix(expand.grid(genes.new, stringsAsFactors = F))
    
    ##
    expand.merged <- cbind(expand, expand.genes)
    expand.merged <- as.data.frame(expand.merged, stringsAsFactors = F)
    expand.merged$ligand_receptor <- as.numeric(expand.merged[,1]) * as.numeric(expand.merged[,2])
    expand.merged$ligand_receptor_name <- paste0(expand.merged$ligand.genes, "_", expand.merged$receptor.genes)
    ##
    expand.merged.rmzero <- subset(expand.merged, expand.merged$ligand_receptor > 0)
    #   rm(expand.merged)
    expand.merged <- NULL
    ##
    expand.merged.rmzero$ligand.subtype <- links2$ligand[i]
    expand.merged.rmzero$receptor.subtype <- links2$receptor[i]
    ##
    
    ## extract only ligand_receptor pairs
    dim(expand.merged.rmzero)
    expand.merged.rmzero <- expand.merged.rmzero[expand.merged.rmzero$ligand_receptor_name %in% RL_pairs$Pair.Name, ]
    dim(expand.merged.rmzero)
    
    ## ligand subtype, receptor subtype, ligand genes, ligand NO cells, receptor genes, receptor NO cells
    cc_log.tmp <- cbind(expand.merged.rmzero$ligand.subtype, expand.merged.rmzero$receptor.subtype,
                        expand.merged.rmzero$ligand.genes, expand.merged.rmzero$ligand_NOcell,
                        expand.merged.rmzero$receptor.genes, expand.merged.rmzero$receptor_NOcell,
                        expand.merged.rmzero$ligand_receptor)
    colnames(cc_log.tmp) <- c("Ligand_subtype", "Receptor_subtype", "Ligand_genes",
                              "Ligand_NO_cells", "Receptor_genes", "Receptor_NO_cells", "Ligand_x_Receptor")
    sum_of_pairs <- sum(as.numeric(as.character(cc_log.tmp[,7])))
    ##
    cc_log <- rbind(cc_log, cc_log.tmp)
    links2$num_of_pairs[i] <- sum_of_pairs
    
    cc_log.tmp <- NULL
    expand.merged <- NULL
    expand.merged.rmzero <- NULL
  }
  
  return(list(nodes2,links2))
}


