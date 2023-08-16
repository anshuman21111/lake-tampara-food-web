# Robustness helper functions to call AB functions and plot
# Modified from code from 10.1038/s41467-021-21824-x

# Function for food web robustness
fw_robustness <- function (MATRIX, BASAL, N_RAND, OUTPUT, FW_NAME, SEQ_NAME, TOT_SUSC, GIVEN) {
  if (GIVEN) {
    # run a given extinction cascade (specified by BASAL)
    mat_results <- collapse_wrap_given(N = MATRIX, basal = BASAL) 
  } else {
    # run standard extinction cascades and subset to the results you want 
    mat_output <- collapse_wrap(N = MATRIX, basal = BASAL, n_rand = N_RAND)
    mat_results <- subset(mat_output, mat_output$id == OUTPUT)
  }
  
  # Y AXIS
  mat_results$nodes_susc <- TOT_SUSC
  mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";")
  mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)
  mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc
  
  # X AXIS
  target_remove <- length(BASAL)
  mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove
  
  # Plot
  root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
  (root + geom_point(shape=19)
    + geom_line() + theme_bw(base_size = 14)
    + xlim(0,1) + ylim(0,1)
    + xlab("Proportion of target species removed")
    + ylab("Proportion of susceptible species remaining")
  )
  ggsave(sprintf("%s_%s.jpeg",FW_NAME, SEQ_NAME), width=6, height=5, units="in", dpi=500)
  
  # Calculate and save robustness (AUC)
  auc <- robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y)
  write.csv(mat_results, sprintf("%s_%s.csv",FW_NAME, SEQ_NAME))
  
  return(auc)
}

# Function for ES Robustness
es_robustness <- function(MATRIX, BASAL, N_RAND, OUTPUT, FW_NAME, SEQ_NAME, TOT_SUSC, GIVEN, SERVICE_NAMES) {
  
  if (GIVEN) {
    # run a given extinction cascade (specified by BASAL)
    mat_results <- collapse_wrap_given(N = MATRIX, basal = BASAL) 
  } else {
    # run standard extinction cascades and subset to the results you want 
    mat_output <- collapse_wrap(N = MATRIX, basal = BASAL, n_rand = N_RAND)
    mat_results <- subset(mat_output, mat_output$id == OUTPUT)
  }
  
  # make a column to track number of services lost at each step (new thing for ES)
  es_pattern <- paste(SERVICE_NAMES, collapse= "|")
  mat_results$service_lost <- ifelse(str_count(string = mat_results$name_lost, pattern = es_pattern),str_count(string = mat_results$name_lost, pattern = es_pattern), NA)
  
  # TO DO - maybe don't need to do this twice for the plotting and the AUC calculation, could combine (?)
  
  NUM_SERVICES <- length(SERVICE_NAMES)
  # track number services remain to use in propES_remain column
  mat_results[["service_lost"]][is.na(mat_results[["service_lost"]])] <- 0
  mat_results$cum_nontarget_lost <- cumsum(mat_results$service_lost)
  mat_results$ES_remain <- c(NUM_SERVICES-cumsum(mat_results$service_lost))
  mat_results$propES_remain <- c(mat_results$ES_remain/NUM_SERVICES)
  
  # ADJUST X AXIS
  target_remove <- length(BASAL)
  mat_results$prop_removed <- mat_results$num_removed_tot / target_remove
  
  root <- ggplot(mat_results, aes(x=prop_removed, y=propES_remain))
  (root + geom_point(shape=19) 
    + geom_line() + theme_bw(base_size = 14) 
    + xlim(0,1) + ylim(0,1)
    + xlab("Proportion of target species removed")
    + ylab("Proportion of ecosystem services remaining")
  )
  ggsave(sprintf("%s_%s_ES.jpeg",FW_NAME, SEQ_NAME), width=6, height=5, units="in", dpi=500)
  write.csv(mat_results, sprintf("%s_%s_ES.csv",FW_NAME, SEQ_NAME))
  
  if (GIVEN){
    # revise output to get ES AUC 
    # this tracks the secondary loss of ecosystem service nodes instead of species
    auc_res <- mat_results %>%
      mutate(service_lost = str_count(string = name_lost, pattern = es_pattern)) %>%
      # add column with the cumulative sum of the services lost
      mutate(service_lost_c = cumsum(service_lost)) %>%
      # overwrite original column of proportion species remaining, for now
      mutate(prop_remain = 1 - (service_lost_c / NUM_SERVICES)) %>%
      auc_wrapper_given()
  } else {
    auc_res <- mat_results %>%
      mutate(service_lost = str_count(string = name_lost, pattern = es_pattern)) %>%
      # add column with the cumulative sum of the services lost
      mutate(service_lost_c = cumsum(service_lost)) %>%
      # overwrite original column of proportion species remaining, for now
      mutate(prop_remain = 1 - (service_lost_c / NUM_SERVICES)) %>%
      auc_wrapper()
  }
  
  return(auc_res$auc)
}

fw_robustness_random <- function(MATRIX, BASAL, N_RAND, TOT_SUSC, FW_NAME) {
  rand <- data.frame(matrix(ncol=2,nrow=N_RAND)) # save robustness results 
  colnames(rand) <- c("Randomization","R")
  rand$Randomization <- seq(1:N_RAND)
  for(i in c(1:N_RAND)){
    seq <- sample(BASAL) # random extinction order 
    mat_results <- collapse_wrap_given(N = MATRIX, basal = seq)
    mat_results$nodes_susc <- TOT_SUSC
    mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";")
    mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each) 
    mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost)/mat_results$nodes_susc
    rand$R[i] <- robust_auc(x = mat_results$prop_removed, y = mat_results$Y)
  }
  write.csv(rand, sprintf("%s_random.csv",FW_NAME))
  # For now return min, mean, and max of the robustness results from random extinctions
  return(list("min"=min(rand$R),"mean"=mean(rand$R),"max"=max(rand$R)))
}

es_robustness_random <- function(MATRIX, BASAL, N_RAND, TOT_SUSC, FW_NAME, SERVICE_NAMES){
  
  es_pattern <- paste(SERVICE_NAMES, collapse= "|")
  NUM_SERVICES <- length(SERVICE_NAMES)
  rand <- data.frame(matrix(ncol=2,nrow=N_RAND)) # save robustness results 
  colnames(rand) <- c("Randomization","R")
  rand$Randomization <- seq(1:N_RAND)
  for(i in c(1:N_RAND)){
    seq <- sample(BASAL)
    mat_output <- collapse_wrap_given(N = MATRIX, basal = seq)
    auc_res <- mat_output %>%
      mutate(service_lost = str_count(string = name_lost, pattern = es_pattern)) %>%
      # add column with the cumulative sum of the services lost
      mutate(service_lost_c = cumsum(service_lost)) %>%
      # overwrite original column of proportion species remaining, for now
      mutate(prop_remain = 1 - (service_lost_c / NUM_SERVICES)) %>%
      auc_wrapper_given()
    rand$R[i] <- auc_res$auc
  }
  write.csv(rand, sprintf("%s_random_ES.csv",FW_NAME))
  return(list("min"=min(rand$R),"mean"=mean(rand$R),"max"=max(rand$R)))
  
}