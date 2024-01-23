library(ggplot2)
library(gridExtra)

#################################################
# IMPORT DATA

### EXPECTED CONTROL DATA
expected <- readxl::read_xlsx("controls_EXPECTED.xlsx")

### FEM
FEM <- read.csv("controls_phased_haplotypes_FEM.csv")
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

cat("There are", length(unique(expected_mixes$SampleID)), "expected mixes")
cat("There are", length(unique(FAPR_mixes$SampleID)), "mixes picked by FapR")
cat("There are", length(unique(FEM_mixes$SampleID)), "mixes picked by FEM")

############################################
  #filter out samples from FAPR (OPTIONAL, could do the same for FEM, maybe even for expected to account for lab...)

result_data_FINAL<-data.frame()

for (freq_threshold in seq(0,0.4, 0.05)){
  
  FAPR_mixes_threshold<- FAPR_mixes[FAPR_mixes$freq > freq_threshold,]
  FEM_mixes_threshold<- FEM_mixes[FEM_mixes$freq > freq_threshold,]
  expected_mixes_threshold<- expected_mixes[expected_mixes$freq > 0.02,] #EMPIRICAL LIMIT OF DETECTION IN LAB !!! CHANGE AS FIT
  
  ############################################
  # COMPARE PRESENCE OF HAPLOS
  
  # Initialize counters
  total_tp_count_FEM <- 0
  total_fp_count_FEM <- 0
  total_fn_count_FEM <- 0
  
  total_tp_count_FAPR <- 0
  total_fp_count_FAPR <- 0
  total_fn_count_FAPR <- 0
  
  for (sample in unique(expected_mixes_threshold$SampleID)){
    
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
  }
  
  # Create a final table
  result_data <- data.frame(
    Method = c("FEM", "FAPR"),
    TP = c(total_tp_count_FEM, total_tp_count_FAPR),
    FP = c(total_fp_count_FEM, total_fp_count_FAPR),
    FN = c(total_fn_count_FEM, total_fn_count_FAPR),
    MAF = freq_threshold
  )
  
  # Calculate metrics
  result_data$Accuracy <- (result_data$TP + result_data$FN) / (result_data$TP + result_data$FP + result_data$FN + result_data$FP)
  result_data$Precision <- result_data$TP / (result_data$TP + result_data$FP)
  result_data$Recall <- result_data$TP / (result_data$TP + result_data$FN)
  result_data$F1_Score <- 2 * result_data$Precision * result_data$Recall / (result_data$Precision + result_data$Recall)
  
  result_data_FINAL <- rbind(result_data_FINAL, result_data)
}  


################################################
# PLOT SHIT

# Convert MAF to a factor for better plotting
result_data_FINAL$MAF <- factor(result_data_FINAL$MAF)

# Plot accuracy
accuracy_plot <- ggplot(result_data_FINAL, aes(x = MAF, y = Accuracy, color = Method)) +
  geom_point() +
  geom_line(aes(group = Method)) +
  labs(title = "Comparison of Accuracy between FAPR and FEM",
       x = "MAF",
       y = "Accuracy") +
  theme_minimal()

# Plot precision
precision_plot <- ggplot(result_data_FINAL, aes(x = MAF, y = Precision, color = Method)) +
  geom_point() +
  geom_line(aes(group = Method)) +
  labs(title = "Comparison of Precision between FAPR and FEM",
       x = "MAF",
       y = "Precision") +
  theme_minimal()

# Plot recall
recall_plot <- ggplot(result_data_FINAL, aes(x = MAF, y = Recall, color = Method)) +
  geom_point() +
  geom_line(aes(group = Method)) +
  labs(title = "Comparison of Recall between FAPR and FEM",
       x = "MAF",
       y = "Recall") +
  theme_minimal()

# Plot F1 Score
f1_plot <- ggplot(result_data_FINAL, aes(x = MAF, y = F1_Score, color = Method)) +
  geom_point() +
  geom_line(aes(group = Method)) +
  labs(title = "Comparison of F1 Score between FAPR and FEM",
       x = "MAF",
       y = "F1 Score") +
  theme_minimal()

grid.arrange(accuracy_plot, precision_plot, recall_plot, f1_plot, ncol = 2)
