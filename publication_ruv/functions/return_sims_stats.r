return_logFC_sd <- function(input_data, 
                            p_threshold = 0.05){
  merged_data <- input_data$merged
  # establish
  stats_cont <- c()
  stats_cont$inter <- mean(abs(merged_data[merged_data$p_intersect <= p_threshold,]$LogFC_sd))
  stats_cont$union <- mean(abs(merged_data[merged_data$p_union <= p_threshold,]$LogFC_sd))
  stats_cont$eR <- mean(abs(merged_data[merged_data$edger_adj_p <= p_threshold,]$LogFC_sd))
  stats_cont$V <- mean(abs(merged_data[merged_data$voom_adj_p <= p_threshold,]$LogFC_sd))
  stats_cont$DE <- mean(abs(merged_data[merged_data$deseq_adj_p <= p_threshold,]$LogFC_sd))
  # what isnt in the intersect, but DE for each method
  stats_cont$eR_rest <- mean(abs(merged_data[merged_data$edger_adj_p <= p_threshold & merged_data$p_intersect > p_threshold,]$LogFC_sd)) #3
  stats_cont$V_rest <- mean(abs(merged_data[merged_data$voom_adj_p <= p_threshold & merged_data$p_intersect > p_threshold,]$LogFC_sd)) #4
  stats_cont$DE_rest <- mean(abs(merged_data[merged_data$deseq_adj_p <= p_threshold & merged_data$p_intersect > p_threshold,]$LogFC_sd)) #5
  # rename
  names(stats_cont) <- c("inter", 
                         "union", 
                         "EdgeR", 
                         "voom", 
                         "DESeq2",
                         "EdgeR_unique",
                         "voom_unique", 
                         "DESeq2_unique"
  )
return(stats_cont)
}

return_stats <- function(input_data, 
                         p_threshold = 0.05){
  merged_data <- input_data$merged
  #summarize results into one table
  stats <- data.frame("ID"= merged_data$ID, 
                      "test" = merged_data$test)
  # establish results for comparison to truth
  stats$inter <- ifelse(merged_data$p_intersect <= p_threshold, 1, 2)
  stats$union <- ifelse(merged_data$p_union <= p_threshold, 1, 2)
  stats$eR <- ifelse(merged_data$edger_adj_p <= p_threshold, 1, 2)
  stats$V <- ifelse(merged_data$voom_adj_p <= p_threshold, 1, 2)
  stats$DE <- ifelse(merged_data$deseq_adj_p <= p_threshold, 1, 2)
  # obtain "truth" statistics
  stats_cont <- cont(stats)
  # names
  names(stats_cont) <- c("inter", 
                         "union", 
                         "EdgeR", 
                         "voom", 
                         "DESeq2")
return(stats_cont)
}

cont_list <- function(input_data, col){
  # counts
  # TP = 1 * 1 = 1
  TP_n <- length(input_data$test[input_data$test == 1 & input_data[col] == 1])
  # TN = -1 * 2 = -2
  TN_n <- length(input_data$test[input_data$test == 2 & input_data[col] == 2])
  #FP = -1 * 1 = -1
  FP_n <- length(input_data$test[input_data$test == 2 & input_data[col] == 1])
  #FN = 1 * 2 = 2
  FN_n <- length(input_data$test[input_data$test == 1 & input_data[col] == 2])
  # stats
  PPV <- TP_n/(TP_n+FP_n); PPV[is.na(PPV)] <- NA
  NPV <- TN_n/(TN_n+FN_n); NPV[is.na(NPV)] <- NA
  FDR <- FP_n/(FP_n+TP_n); FDR[is.na(FDR)] <- NA
  sens <- TP_n/(TP_n+FN_n); sens[is.na(sens)] <- NA
  spec <- TN_n/(TN_n+FP_n); spec[is.na(spec)] <- NA
  F1 <- 2*((PPV*sens)/(PPV+sens))
  ACC <- (TP_n+TN_n)/(TP_n+TN_n+FP_n+FN_n)
  results_return <- c(TP_n, TN_n, FP_n, FN_n, PPV, NPV, FDR, sens, spec, F1, ACC)
  names(results_return) <- c("TP_n", "TN_n", "FP_n", "FN_n", "PPV", "NPV", "FDR", "sens", "spec", "F1", "ACC")
return(results_return)
}

cont <- function(input_data){
  list_of_results <- lapply(3:ncol(input_data), function(x)
    cont_list(input_data, x))
return(list_of_results)
}

get_mean_FC <- function(which_data, which_experiment){
  count_me <- which_data[[which_experiment]]$sims
  # dimensions and names
  sim_no <- length(count_me)
  # reformat
  return_r <- lapply(1:sim_no, function(x) 
    get_r(count_me, x))
  return_r2 <- t(data.frame(return_r))
  rownames(return_r2) <- 1:sim_no
return(return_r2)
}

get_r <- function(count_data, which_sim){
  # obtain LogFC values
  DE.short <- count_data[[which_sim]]$mde_all$deseq$short_results[[1]]
  DE.short <- data.frame("logFC_DE"=DE.short$logFC,
                         "ID"= rownames(DE.short))
  
  VM.short <- count_data[[which_sim]]$mde_all$voom$short_results[[1]]
  VM.short <- data.frame("logFC_VM"=VM.short$logFC,
                         "ID"= rownames(VM.short))
  
  ER.short <- count_data[[which_sim]]$mde_all$edger$short_results[[1]]
  ER.short <- data.frame("logFC_ER"=ER.short$logFC,
                         "ID"= rownames(ER.short))
  # compare modelled FC
  all_FC <- merge(merge(DE.short, VM.short, by="ID"), ER.short, by="ID")
  my_r <- c(summary(lm(all_FC$logFC_VM~all_FC$logFC_DE))$r.squared,
            summary(lm(all_FC$logFC_VM~all_FC$logFC_ER))$r.squared,
            summary(lm(all_FC$logFC_ER~all_FC$logFC_DE))$r.squared)
  names(my_r) <- c("VM.DE", "VM.ER", "ER.DE")
return(my_r)
}

clean_stats <- function(which_data, which_experiment, isLogFC = FALSE){
  # format difference for logFC stats
  if(isLogFC == FALSE){
    count_me <- which_data[[which_experiment]]$stats
  }
  if(isLogFC == TRUE){
    count_me <- which_data[[which_experiment]]$sd_stats
  }
  # dimensions and names
  sim_no <- length(count_me)
  sim_names <- names(count_me[[1]])
  sim_measures <- length(sim_names)
  if(isLogFC == TRUE){
    sim_stats_names <- "logFC_sd"
    sim_stats_size <- 1
    sim_stats_clean <- sapply(1:sim_no, function(x)
      data.frame(DataFrame(count_me[[x]])))
    sim_stats_clean4 <- t(sim_stats_clean)
    
  }
  if(isLogFC == FALSE){
    sim_stats_names <- names(count_me[[1]][[1]])
    sim_stats_size <- length(names(count_me[[1]][[1]]))
    #step 1 - reformat
    sim_stats_clean <- lapply(1:sim_no, function(x)
      data.frame(DataFrame(count_me[[x]])))
    #step 2
    sim_stats_clean2 <- lapply(1:sim_stats_size, function(x)
      t(sapply(1:sim_no, function(i)
        sim_stats_clean[[i]][x,])))
    # step 3 - tidy list of vectors
    sim_stats_clean3 <- lapply(1:sim_stats_size, function(x)
      cbind(sapply(1:sim_measures, function(i)
        unlist(sim_stats_clean2[[x]][,i]))))
    # step 4 - rename
    sim_stats_clean4 <- lapply(1:length(sim_stats_clean3), function(x)
      setNames(data.frame(sim_stats_clean3[[x]]),
               sim_names))
    names(sim_stats_clean4) <- sim_stats_names
  }
return(sim_stats_clean4)
} 
