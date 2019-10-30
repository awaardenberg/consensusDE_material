wrap_all <- function(path_to_data,
                     p_thresh = 0.05,
                     replicates = 3,
                     ruv = FALSE,
                     norm_method = norm_method,
                     sim_number = 10,
                     de_number = 500,
                     de_size = 10000){
  # for reproducibility
  set.seed(1234)
  # run number of simulations
  sims <- lapply(1:sim_number, function(x) 
    run_sims(real_data = path_to_data,
             replicates = replicates,
             den_n = de_size, 
             n_de = de_number, 
             ruv = ruv,
             norm_method = norm_method))
  # return stats
  sims_stats <- lapply(1:length(sims), function(x) 
    return_stats(sims[[x]],
                 p_threshold = p_thresh))
  # return logFC stats
  sd_stats <- lapply(1:length(sims), function(x) 
    return_logFC_sd(sims[[x]],
                    p_threshold = p_thresh))
  
return(list("sims" = sims,
            "stats" = sims_stats,
            "sd_stats" = sd_stats))
}

# dont need to provide sample_table - will build based on balance design
run_sims <- function(real_data = NULL,
                     sample_table = NULL, 
                     replicates = 3,
                     den_n = 10000, 
                     n_de = 500, 
                     ruv = FALSE,
                     norm_method = norm_method){
  # define sample_table by size and de number
  sample_table <- data.frame("file"=paste(rep(c("G1_rep", "G2_rep"), each = replicates), 1:replicates, sep=""),
                             "group"=rep(c("G1", "G2"), each = replicates))
  
  # divide by 2 = balanced experiment
  n_de <- n_de/2
  # obtain parameters in list
  par_list <- estimate.sim.params(real.counts = real_data)
  # this will estimate parameters for analysis
  sim <- make.sim.data.sd(N = den_n, 
                          param = par_list,
                          ndeg = c(n_de, n_de),
                          samples = c(replicates, replicates))
  # truth set
  true_de <- which(sim$truedeg != 0)
  neg_de <- which(sim$truedeg == 0)
  # data
  syn_data <- sim$simdata[9:ncol(sim$simdata)]
  # put into summarized experiment
  se <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(syn_data)))
  colData(se) <- S4Vectors::DataFrame(sample_table)
  colnames(se) <- sample_table$file
  
  # filter results
  keep <- filterByExpr(assays(se)$counts, group=colData(se)$group)
  se <- se[rownames(se)[keep] ,]
  # run multi_de_pairs from consensusDE
  mde <- multi_de_pairs(summarized = se,
                        ruv_correct = ruv,
                        norm_method = norm_method,
                        verbose = TRUE)
  # obtain merged data results
  merged_data <- mde$merged[[1]]
  
  # establish true positive and true negative sets
  TP <- data.frame("ID" = names(true_de),
                   "test" = 1)
  TN <- data.frame("ID" = names(neg_de),
                   "test" = 2)
  TP_TN <- rbind(TP, TN)
  output <- merge(merged_data, TP_TN, by = "ID")
  output <- output[order(output$rank_sum, decreasing = FALSE),]
return(list("merged" = output, 
            "TP" = true_de,
            "TN" = neg_de,
            "syn_data" = syn_data,
            "mde_all" = mde))
}