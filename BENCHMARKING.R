
#import and format data
expected <- readxl::read_xlsx("controls_EXPECTED.xlsx")

FEM <- read.csv("controls_phased_haplotypes_FEM.csv")
colnames(FEM)<- c("SampleID", "genotypes_FEM", "freq_FEM")

FAPR <- read.csv("controls_fapr_phased_haplos_correct_only.csv")
FAPR <- FAPR[,c("SampleID", "haplotype", "HAPLO_FREQ_RECALC")]
colnames(FAPR)<- c("SampleID", "genotypes_FAPR", "freq_FAPR")


# COMPARE PREENCE OF HAPLOS
# remove monoclonal samples from all dfs first NOT DONE YET

result_FAPR <- c()
result_FEM <- c()

for (sample in expected$SampleID){
  
  expected_chunk <- expected[expected$SampleID == sample,]
    
  FEM_chunk <- FEM[FEM$SampleID == sample,]
  FAPR_chunk <- FAPR[FAPR$SampleID == sample,]
  
  ok1 <- FAPR_chunk$genotypes %in% expected_chunk$genotypes # what haplos from FAPR are present in expected, per sample?
  ok2 <- FEM_chunk$genotypes %in% expected_chunk$genotypes
  
  result_FAPR <- c(ok1, result_FAPR)
  result_FEM <- c(ok2, result_FEM)
  
}

print(sum(result_FAPR == TRUE)) # TP
print(sum(result_FAPR == FALSE)) # FP

print(sum(result_FEM == FALSE)) # FP
print(sum(result_FEM == TRUE)) # TP
