
#import and format data
expected <- readxl::read_xlsx("controls_EXPECTED.xlsx")

FEM <- read.csv("controls_phased_haplotypes_FEM.csv")
colnames(FEM)<- c("SampleID", "genotypes_FEM", "freq_FEM")

FAPR <- read.csv("controls_fapr_phased_haplos_correct_only.csv") #correct only o all? probar ambas
FAPR <- FAPR[,c("SampleID", "haplotype", "HAPLO_FREQ_RECALC")]
colnames(FAPR)<- c("SampleID", "genotypes_FAPR", "freq_FAPR")


# remove monoclonal samples from all dfs first
monoclonal_samples <- expected[expected$freq == 1,]["SampleID"]
monoclonal_samples <- monoclonal_samples[!is.na(monoclonal_samples$SampleID),]

expected_mixes <- expected[!expected$SampleID %in% monoclonal_samples$SampleID,]
FEM_mixes <- FEM[!FEM$SampleID %in% monoclonal_samples$SampleID,]
FAPR_mixes <- FAPR[!FAPR$SampleID %in% monoclonal_samples$SampleID,]

cat("There are", length(unique(expected_mixes$SampleID)), "expected mixes")
cat("There are", length(unique(FAPR_mixes$SampleID)), "mixes picked by FapR")
cat("There are", length(unique(FEM_mixes$SampleID)), "mixes picked by FEM")


# COMPARE PRESENCE OF HAPLOS
result_FAPR <- c()
result_FEM <- c()

for (sample in expected$SampleID){
  
  expected_chunk <- expected_mixes[expected_mixes$SampleID == sample,]
    
  FEM_chunk <- FEM_mixes[FEM_mixes$SampleID == sample,]
  FAPR_chunk <- FAPR_mixes[FAPR_mixes$SampleID == sample,]
  
  ok1 <- FAPR_chunk$genotypes %in% expected_chunk$genotypes # what haplos from FAPR are present in expected, per sample?
  ok2 <- FEM_chunk$genotypes %in% expected_chunk$genotypes
  
  result_FAPR <- c(ok1, result_FAPR)
  result_FEM <- c(ok2, result_FEM)
  
}

print(sum(result_FAPR == TRUE)) # TP
print(sum(result_FAPR == FALSE)) # FP

print(sum(result_FEM == TRUE)) # TP
print(sum(result_FEM == FALSE)) # FP

