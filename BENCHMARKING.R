
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
#filter out samples from FAPR (OPTIONAL, could do the same for FEM)

freq_threshold = 0.05
FAPR_mixes<- FAPR_mixes[FAPR_mixes$freq > freq_threshold,]

############################################
# COMPARE PRESENCE OF HAPLOS

result_TP_FP_FAPR <- c()
result_TP_FP_FEM <- c()

for (sample in unique(expected_mixes$SampleID)){
  
  expected_chunk <- expected_mixes[expected_mixes$SampleID == sample,]
    
  FEM_chunk <- FEM_mixes[FEM_mixes$SampleID == sample,]
  FAPR_chunk <- FAPR_mixes[FAPR_mixes$SampleID == sample,]

#####
    
  TP_FP_fapr <-  FAPR_chunk$genotypes_FAPR %in% expected_chunk$genotypes    # what haplos from FAPR are present in expected, per sample?
  
  if (length(TP_FP_fapr) == 0){
    TP_FP_fapr <- rep(FALSE, length(expected_chunk$genotypes))
  }
  
  TP_FP_fem <- FEM_chunk$genotypes_FEM %in% expected_chunk$genotypes 
  
  if (length(TP_FP_fem) == 0) {
    TP_FP_fem <- rep(FALSE, length(expected_chunk$genotypes))
  }
  
  #print(c(sample, TP_FP_fapr))
  
  result_TP_FP_FAPR <- c(result_TP_FP_FAPR, TP_FP_fapr)
  result_TP_FP_FEM <- c(result_TP_FP_FEM, TP_FP_fem)
  
#####
}

print(sum(result_TP_FP_FAPR == TRUE)) # TP
print(sum(result_TP_FP_FAPR == FALSE)) # FP

print(sum(result_TP_FP_FEM == TRUE)) # TP
print(sum(result_TP_FP_FEM == FALSE)) # FP

length(expected_mixes$SampleID)
