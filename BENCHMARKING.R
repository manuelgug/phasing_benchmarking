library(ggplot2)
library(gridExtra)

#################################################
# IMPORT DATA

### EXPECTED CONTROL DATA
expected <- readxl::read_xlsx("controls_EXPECTED.xlsx")

### FEM
FEM <- read.csv("controls_phased_haplotypes_FEM_ALL_CONTROLS.csv")
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

### FAPR
FAPR <- read.csv("controls_fapr_phased_haplos.csv") #correct only o all? probar ambas
FAPR <- FAPR[,c("SampleID", "haplotype", "HAPLO_FREQ_RECALC")]
colnames(FAPR)<- c("SampleID", "genotypes_FAPR", "freq_FAPR")

##########################################
# remove monoclonal samples from all dfs

monoclonal_samples <- expected[expected$freq == 1,]["SampleID"]
monoclonal_samples <- monoclonal_samples[!is.na(monoclonal_samples$SampleID),]

expected_mixes <- expected[!expected$SampleID %in% monoclonal_samples$SampleID,]
FEM_mixes <- FEM[!FEM$SampleID %in% monoclonal_samples$SampleID,]
FAPR_mixes <- FAPR[!FAPR$SampleID %in% monoclonal_samples$SampleID,]

############################################
### MAIN FUNCTION-------------------------------------------------------------------

compareMethods <- function(FAPR_mixes, FEM_mixes, expected_mixes) {
  result_data_FINAL <- data.frame()
  
  ############################################
  #filter out samples given a freq threshold (results will reflect a range of MAF thresholds)
  
  for (freq_threshold in seq(0, 0.4, 0.01)) {
    FAPR_mixes_threshold <- FAPR_mixes[FAPR_mixes$freq > freq_threshold, ]
    FEM_mixes_threshold <- FEM_mixes[FEM_mixes$freq > freq_threshold, ]
    expected_mixes_threshold <- expected_mixes[expected_mixes$freq > 0.02, ] #EMPIRICAL LIMIT OF DETECTION IN LAB !!! CHANGE AS FIT
    
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
    
    for (sample in unique(expected_mixes_threshold$SampleID)){ #quitar nnúmero entre []
      
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
#### PLOT FUNCTION------------------------------------------------------------------
plotMetricsGrid <- function(result_data_FINAL_all_parasitaemias, plot_title) {
  # Convert MAF to a factor for better plotting
  result_data_FINAL_all_parasitaemias$MAF <- factor(result_data_FINAL_all_parasitaemias$MAF)
  
  # Plot accuracy
  accuracy_plot <- ggplot(result_data_FINAL_all_parasitaemias, aes(x = MAF, y = Accuracy, color = Method)) +
    geom_point() +
    geom_line(aes(group = Method)) +
    labs(title = "Accuracy",
         x = "MAF",
         y = "Accuracy") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(0, 1)
  
  # Plot precision
  precision_plot <- ggplot(result_data_FINAL_all_parasitaemias, aes(x = MAF, y = Precision, color = Method)) +
    geom_point() +
    geom_line(aes(group = Method)) +
    labs(title = "Precision",
         x = "MAF",
         y = "Precision") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(0, 1)
  
  # Plot recall
  recall_plot <- ggplot(result_data_FINAL_all_parasitaemias, aes(x = MAF, y = Recall, color = Method)) +
    geom_point() +
    geom_line(aes(group = Method)) +
    labs(title = "Recall",
         x = "MAF",
         y = "Recall") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(0, 1)
  
  # Plot F1 Score
  f1_plot <- ggplot(result_data_FINAL_all_parasitaemias, aes(x = MAF, y = F1_Score, color = Method)) +
    geom_point() +
    geom_line(aes(group = Method)) +
    labs(title = "F1 Score",
         x = "MAF",
         y = "F1 Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(0, 1)
  
  # RMSE of freq
  rmse_plot <- ggplot(result_data_FINAL_all_parasitaemias, aes(x = MAF, y = RMSE, color = Method)) +
    geom_point() +
    geom_line(aes(group = Method)) +
    labs(title = "RMSE of freq",
         x = "MAF",
         y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(0, 1)
  
  # MAPE of freq
  mape_plot <- ggplot(result_data_FINAL_all_parasitaemias, aes(x = MAF, y = MAPE, color = Method)) +
    geom_point() +
    geom_line(aes(group = Method)) +
    labs(title = "MAPE of freq",
         x = "MAF",
         y = "MAPE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Create a grid
  #grid.arrange(accuracy_plot, precision_plot, recall_plot, f1_plot, rmse_plot, mape_plot, ncol = 3, top = paste("Parasitaemia =", plot_title))
  
  # Save the grid plot
  ggsave(
    paste0("benchmarking_parasitaemia", plot_title, ".png"),
    arrangeGrob(accuracy_plot, precision_plot, recall_plot, f1_plot, rmse_plot, mape_plot, ncol = 3, top = paste("Parasitaemia =", plot_title)),
    width = 20, 
    height = 12,
    dpi = 300 
  )
  
  }
####----------------------------------------------------------------------------------------

###################
# PLOTS
for (df in 1:length(result_data_FINAL_all_parasitaemias_separatedly)){
  
  plotMetricsGrid(result_data_FINAL_all_parasitaemias_separatedly[[df]], names(result_data_FINAL_all_parasitaemias_separatedly[df]))
  
}


