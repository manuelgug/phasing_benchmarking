library(ggplot2)
library(gridExtra)
library(dplyr)

#################################################
# IMPORT DATA

### EXPECTED CONTROL DATA (only used to label de plots; filtering happening in the MAF step)----
expected <- readxl::read_xlsx("inputs/controls_EXPECTED.xlsx")



### FEM -----
FEM <- read.csv("run_FreqEstimationModel/controls_phased_haplotypes_FEM.csv")
colnames(FEM)<- c("SampleID", "genotypes_FEM", "freq_FEM")

# Function to process rows with brackets: because FEM can't handle multiallelic samples, I'm splitting the "all other alleles" into cloned rows with each genotype to test for accuracy
process_row <- function(row) {
  if (grepl("\\[.*\\]", row$genotypes_FEM)) {
    bracket_content <- gsub("[\\[\\]]", "", regmatches(row$genotypes_FEM, regexpr("\\[.*\\]", row$genotypes_FEM)))
    row$genotypes_FEM <- gsub("\\[.*\\]", "", row$genotypes_FEM)
    new_rows <- data.frame(
      SampleID = rep(row$SampleID, nchar(bracket_content)),
      genotypes_FEM = paste0(row$genotypes_FEM, "_", strsplit(bracket_content, '')[[1]]),
      freq_FEM = rep(row$freq_FEM, nchar(bracket_content)),
      stringsAsFactors = FALSE
    )
    #remove shit rows
    new_rows <- new_rows[!grepl("\\]$", new_rows$genotypes_FEM), ]
    new_rows <- new_rows[!grepl("\\[$", new_rows$genotypes_FEM), ]
    new_rows$genotypes_FEM <- gsub("__", "_", new_rows$genotypes_FEM)  # Replace double underscore with a single one
    new_rows <- new_rows[!grepl("_$", new_rows$genotypes_FEM), ]
    return(rbind(new_rows))
  } else {
    return(row)
  }
}

#loop through rows:
ok <- c()

for (row in 1:length(FEM$SampleID)){
  processed_rows <- process_row(FEM[row,])
  
  ok<-rbind(ok, processed_rows)
}

FEM <- ok



### FAPR ----

FAPR <- read.csv("run_FapR/phased_haplos_relabun_adjusted.csv")
FAPR <- FAPR[,c("SampleID", "haplotype", "HAPLO_FREQ_RECALC")]
colnames(FAPR)<- c("SampleID", "genotypes_FAPR", "freq_FAPR")

##########################################
# remove monoclonal samples from all dfs?

biallelic_benchmark <- FALSE #PUT_BOOLEAN_HERE

if (is.null(biallelic_benchmark)){
  print("SPECIFY TRUE OR FALSE FOR BIALLELIC_BENCHMARK. FALSE MEANS THAT ALL POLIALLELIC LOCI ARE GOING TO BE USED")
}

if (biallelic_benchmark == TRUE){
  
  print("BENCHMARKING ONLY BIALLELIC SAMPLES")
  
  #if benchmarking only for biallelic samples, use this filter instead:
  sample_table <- table(expected$SampleID)
  biallelic_controls <- names(sample_table[sample_table == 2])
  biallelic_controls <- as.data.frame(biallelic_controls)
  
  #if benchmarking for biallelic samples only, run the following 3 lines, else skip
  expected_mixes <- expected[expected$SampleID %in% biallelic_controls$biallelic_controls,] 
  FEM_mixes <- FEM[FEM$SampleID %in% biallelic_controls$biallelic_controls,]
  FEM_mixes <- FEM_mixes %>% distinct(SampleID, freq_FEM, .keep_all = TRUE)
  FAPR_mixes <- FAPR[FAPR$SampleID %in% biallelic_controls$biallelic_controls,]

}else{ #benchmark all samples, not only biallelic ones
  
  print("BENCHMARKING ALL POLIALLELIC SAMPLES")
  
  monoclonal_samples <- expected[expected$freq == 1,]["SampleID"] 
  
  expected_mixes <- expected[!expected$SampleID %in% monoclonal_samples$SampleID,] 
  FEM_mixes <- FEM[!FEM$SampleID %in% monoclonal_samples$SampleID,]
  FAPR_mixes <- FAPR[!FAPR$SampleID %in% monoclonal_samples$SampleID,]
  
}


### keep polyclonal controls only (remove rows with haplo freq of 1, meaning mnoclonal) ------
polyclonal_controls  <- unique(expected[expected$freq < 1,]$SampleID)

expected <- expected[expected$SampleID %in% polyclonal_controls,]

FEM_mixes <- FEM_mixes[FEM_mixes$SampleID %in% polyclonal_controls,]
FAPR_mixes <- FAPR_mixes[FAPR_mixes$SampleID %in% polyclonal_controls,]



############################################
### MAIN FUNCTION-------------------------------------------------------------------

compareMethods <- function(FAPR_mixes, FEM_mixes, expected_mixes) {
  result_data_FINAL <- data.frame()
  
  ############################################
  #filter out samples given a freq threshold (results will reflect a range of MAF thresholds)
  
  for (freq_threshold in seq(0, 0.4, 0.01)) {
    
    #freq_threshold <- 0.1
    
    FAPR_mixes_threshold <- FAPR_mixes[FAPR_mixes$freq > freq_threshold, ]
    FEM_mixes_threshold <- FEM_mixes[FEM_mixes$freq > freq_threshold, ]
    expected_mixes_threshold <- expected_mixes[expected_mixes$freq > 0, ] #EMPIRICAL LIMIT OF DETECTION IN LAB is 0.02 !!! CHANGE AS FIT
    
    ############################################
    # COMPARE PRESENCE OF HAPLOS
    
    # Initialize counters
    total_tp_count_FEM <- 0
    total_fp_count_FEM <- 0
    total_fn_count_FEM <- 0
    
    total_tp_count_FAPR <- 0
    total_fp_count_FAPR <- 0
    total_fn_count_FAPR <- 0
    
    merged_fapr_ <- data.frame()
    merged_fem_ <- data.frame()
    
    for (sample in unique(expected_mixes_threshold$SampleID)){ 
      
      #sample <- unique(expected_mixes_threshold$SampleID)[23]
      
      expected_chunk <- expected_mixes_threshold[expected_mixes_threshold$SampleID == sample,]
      
      FEM_chunk <- FEM_mixes_threshold[FEM_mixes_threshold$SampleID == sample,]
      FAPR_chunk <- FAPR_mixes_threshold[FAPR_mixes_threshold$SampleID == sample,]
      
      # Find common genotypes between FAPR and expected
      common_genotypes_FAPR <- intersect(FAPR_chunk$genotypes, expected_chunk$genotypes)
      
      tp_count_FAPR <- sum(common_genotypes_FAPR %in% expected_chunk$genotypes)
      fp_count_FAPR <- sum(!FAPR_chunk$genotypes %in% expected_chunk$genotypes)
      fn_count_FAPR <- sum(!expected_chunk$genotypes %in% common_genotypes_FAPR)
      
      # Update total counts for FAPR
      total_tp_count_FAPR <- total_tp_count_FAPR + tp_count_FAPR
      total_fp_count_FAPR <- total_fp_count_FAPR + fp_count_FAPR
      total_fn_count_FAPR <- total_fn_count_FAPR + fn_count_FAPR
      
      # Find common genotypes between FEM and expected
      common_genotypes_FEM <- intersect(FEM_chunk$genotypes, expected_chunk$genotypes)
      
      tp_count_FEM <- sum(common_genotypes_FEM %in% expected_chunk$genotypes)
      fp_count_FEM <- sum(!FEM_chunk$genotypes %in% expected_chunk$genotypes)
      fn_count_FEM <- sum(!expected_chunk$genotypes %in% common_genotypes_FEM)
      
      # Update total counts for FEM
      total_tp_count_FEM <- total_tp_count_FEM + tp_count_FEM
      total_fp_count_FEM <- total_fp_count_FEM + fp_count_FEM
      total_fn_count_FEM <- total_fn_count_FEM + fn_count_FEM
      
      ###################################################################
      # calculate difference in frequency from true positives
      if (tp_count_FAPR > 0){
        merged_fapr <-merge(expected_chunk, FAPR_chunk, by.x = "genotypes", by.y = "genotypes_FAPR")
        merged_fapr_ <- rbind(merged_fapr_, merged_fapr)
      }
      
      if (tp_count_FEM > 0){
        merged_fem <-merge(expected_chunk, FEM_chunk, by.x = "genotypes", by.y = "genotypes_FEM")
        merged_fem_ <- rbind(merged_fem_, merged_fem)
      }
      ###################################################################
    }
    
    #root mean square error
    rmse_FAPR <- sqrt(mean((merged_fapr_$freq - merged_fapr_$freq_FAPR)^2)) 
    rmse_FEM <- sqrt(mean((merged_fem_$freq - merged_fem_$freq_FEM)^2)) 
    
    # Calculate MAPE (mean average percentage error; takes into account sample size)
    mape_FAPR <- mean(abs((merged_fapr_$freq - merged_fapr_$freq_FAPR) / merged_fapr_$freq_FAPR) * 100, na.rm = TRUE)
    mape_FEM <- mean(abs((merged_fem_$freq - merged_fem_$freq_FEM) / merged_fem_$freq_FEM) * 100, na.rm = TRUE)
    ###################################################################
    
    # Create a final table
    result_data <- data.frame(
      Method = c("FEM", "FapR"),
      TP = c(total_tp_count_FEM, total_tp_count_FAPR),
      FP = c(total_fp_count_FEM, total_fp_count_FAPR),
      FN = c(total_fn_count_FEM, total_fn_count_FAPR),
      MAF = freq_threshold,
      RMSE = c(rmse_FEM, rmse_FAPR),
      MAPE = c(mape_FEM, mape_FAPR)
    )
    
    # Calculate metrics
    result_data$Accuracy <- (result_data$TP + result_data$FN) / (result_data$TP + result_data$FP + result_data$FN + result_data$FP)
    result_data$Precision <- result_data$TP / (result_data$TP + result_data$FP)
    result_data$Recall <- result_data$TP / (result_data$TP + result_data$FN)
    result_data$F1_Score <- 2 * result_data$Precision * result_data$Recall / (result_data$Precision + result_data$Recall)
    
    result_data_FINAL <- rbind(result_data_FINAL, result_data)
  }
  
  return(result_data_FINAL)
}

###----------------------------------------------------------------------------

#########################################################
# Run comparisons with each parasitaemia and all parasitaemias as a whole

#all
result_data_FINAL_all_parasitaemias <- compareMethods(FAPR_mixes, FEM_mixes, expected_mixes)

#parasitamias separatedly
result_data_FINAL_all_parasitaemias_separatedly <- list()

for (parasitaemia in unique(expected_mixes$parasitaemia)){
  
  parasitaemia_subset <- expected_mixes[expected_mixes$parasitaemia == parasitaemia,]$SampleID
  
  expected_mixes_parasitaemia_subset <- expected_mixes[expected_mixes$SampleID %in% parasitaemia_subset, ]
  
  FAPR_mixes_parasitaemia_subset <- FAPR_mixes[FAPR_mixes$SampleID %in% parasitaemia_subset, ]
  FEM_mixes_parasitaemia_subset <- FEM_mixes[FEM_mixes$SampleID %in% parasitaemia_subset, ]
  
  result_data <- compareMethods(FAPR_mixes_parasitaemia_subset, FEM_mixes_parasitaemia_subset, expected_mixes_parasitaemia_subset)
  
  # Store the result_data in the list with the name as parasitaemia
  result_data_FINAL_all_parasitaemias_separatedly[[as.character(parasitaemia)]] <- result_data
}

#add the all parasitamias df to the list
result_data_FINAL_all_parasitaemias_separatedly[["ALL_parasitaemias"]] <- result_data_FINAL_all_parasitaemias


##################################################################3

# PLOTS

#add parasitaemia column
for (df in 1:length(result_data_FINAL_all_parasitaemias_separatedly)){

  result_data_FINAL_all_parasitaemias_separatedly[[df]]$parasitaemia <- names(result_data_FINAL_all_parasitaemias_separatedly[df])
  
}

# Filter data frames in the list for FEM and FapR
filtered_data_FEM <- lapply(result_data_FINAL_all_parasitaemias_separatedly, function(df) df[df$Method == "FEM", ])
filtered_data_FEM <- do.call(rbind, filtered_data_FEM)
filtered_data_FAPR <- lapply(result_data_FINAL_all_parasitaemias_separatedly, function(df) df[df$Method == "FapR", ])
filtered_data_FAPR <- do.call(rbind, filtered_data_FAPR)


metrics <- c("Accuracy", "Precision", "Recall", "F1_Score", "RMSE", "MAPE")


# Create a 6-panel figure
fig <- grid.arrange(
  arrangeGrob(
    grobs = lapply(metrics, function(metric) {
      p <- ggplot(filtered_data_FEM, aes(x = MAF, y = .data[[metric]], color = factor(parasitaemia))) +
        geom_point() +
        geom_line(aes(group = interaction(Method, parasitaemia))) +
        labs(title = metric, x = "MAF", y = metric) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      
      # Add ylim(0, 1) for the first 5 metrics
      if (metric %in% metrics[1:5]) {
        p <- p + ylim(0, 1)
      }
      
      p
    }),
    ncol = 3
  ),
  top = "FEM"
)

ggsave("benchmark_FEM_metrics_plot.png", fig, width = 20, height = 12)


# Create a 6-panel figure
fig <- grid.arrange(
  arrangeGrob(
    grobs = lapply(metrics, function(metric) {
      p <- ggplot(filtered_data_FAPR, aes(x = MAF, y = .data[[metric]], color = factor(parasitaemia))) +
        geom_point() +
        geom_line(aes(group = interaction(Method, parasitaemia))) +
        labs(title = metric, x = "MAF", y = metric) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      
      # Add ylim(0, 1) for the first 5 metrics
      if (metric %in% metrics[1:5]) {
        p <- p + ylim(0, 1)
      }
      
      p
    }),
    ncol = 3
  ),
  top = "FapR"
)


ggsave("benchmark_FAPR_metrics_plot_relabund_adjusted.png", fig, width = 20, height = 12)

