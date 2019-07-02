# obtain overlaps and set sizes for the airway analysis
return_cde_stats <- function(input_pairs, p_cut = 0.05){
  merged_data <- input_pairs$merged[[1]]
  # edgeR
  ER <- merged_data[merged_data$edger_adj_p <= p_cut,]
  # DESeq2
  DE <- merged_data[merged_data$deseq_adj_p <= p_cut,]
  # Voom/limma
  VM <- merged_data[merged_data$voom_adj_p <= p_cut,]
  
  ER_ID <- as.character(ER$ID)
  DE_ID <- as.character(DE$ID)
  VM_ID <- as.character(VM$ID)
  ER.DE.VM_ID <- (intersect(intersect(ER_ID, DE_ID), VM_ID))
  
  # proportions of intersect for each method
  ER_prop <- length(intersect(ER_ID, ER.DE.VM_ID))/length(ER_ID)
  DE_prop <- length(intersect(DE_ID, ER.DE.VM_ID))/length(DE_ID)
  VM_prop <- length(intersect(VM_ID, ER.DE.VM_ID))/length(VM_ID)
  
  props <- list("ER" = ER_prop,
                "DE" = DE_prop,
                "VM" = VM_prop)
  
  venn.inputs <- list(ER_ID, DE_ID, VM_ID)
  names(venn.inputs) <- c("ER", "DE", "VM")
  
return(list("venn_input" = venn.inputs,
            "intersect" = ER.DE.VM_ID,
            "props" = props,
            "ER" = ER,
            "DE" = DE,
            "VM" = VM))
}

get_merged_logfc <- function(cde_object){
  # compare the LogFC values.
  DE.short <- cde_object$deseq$short_results[[1]]
  DE.short <- data.frame("logFC_DE"=DE.short$logFC,
                         "ID"= rownames(DE.short))
  
  VM.short <- cde_object$voom$short_results[[1]]
  VM.short <- data.frame("logFC_VM"=VM.short$logFC,
                         "ID"= rownames(VM.short))
  
  ER.short <- cde_object$edger$short_results[[1]]
  ER.short <- data.frame("logFC_ER"=ER.short$logFC,
                         "ID"= rownames(ER.short))
  # compare the modelled FC.
  all_FC <- merge(merge(DE.short, VM.short, by="ID"), ER.short, by="ID")
return(all_FC)
}
