
#import and format data
expected <- readxl::read_xlsx("controls_EXPECTED.xlsx")

FEM <- read.csv("controls_phased_haplotypes_FEM.csv")
colnames(FEM)<- c("SampleID", "genotypes_FEM", "freq_FEM")

FAPR <- read.csv("controls_fapr_phased_haplos.csv")
FAPR <- FAPR[,c("SampleID", "haplotype", "HAPLO_FREQ_RECALC")]
colnames(FAPR)<- c("SampleID", "genotypes_FAPR", "freq_FAPR")


#merged_table <- merge(merge(expected, FEM, by = "SampleID", all.x = TRUE), FAPR, by = "SampleID", all.x = TRUE)

#compare presence of haplos 
result <- c()

for (sample in expected$SampleID){
  
  expected_chunk <- expected[expected$SampleID == sample,]
    
  FEM_chunk <- FEM[FEM$SampleID == sample,]
  FAPR_chunk <- FAPR[FAPR$SampleID == sample,]
  
  ok <- FAPR_chunk$genotypes %in% expected_chunk$genotypes
  result <- c(ok, result)
}

print(sum(result) / length(expected$genotypes))

