

# this script converts resmarker_table from mad4hatter >= v0.1.8 to input for FreqEstimationModel (https://github.com/aimeertaylor/FreqEstimationModel/)

library(dplyr)
library(tidyr)

resmarkers_table <- read.csv("combined_mixture_controls_resmarker_table_16_17_21_22_27_manuel.csv")

#####################################################################3
#subset relevant markers
markers_to_phase <- c("dhps_431", "dhps_437", "dhps_540", "dhps_581", "dhfr_51", "dhfr_59", "dhfr_108")

#create resmarker column
resmarkers_table$resmarker <-  paste0(resmarkers_table$Gene, "_", resmarkers_table$CodonID)

resmarkers_table <- resmarkers_table %>%
  filter(grepl(paste(markers_to_phase, collapse = "|"), resmarker))

#for resmarkers with 2 amplicons (as dhps_581), keep the one with the highest amount of reads for each variant
resmarkers_table <- resmarkers_table %>%
  group_by(SampleID, resmarker, AA) %>%
  filter(Reads == max(Reads)) %>%
  ungroup()

#if there is no norm.read.locus column, create it
resmarkers_table<- resmarkers_table %>%
  group_by(SampleID,resmarker, CodonStart) %>%
  mutate(norm.reads.locus = Reads/sum(Reads)) %>%
  mutate(n.alleles = n())
#####################################################################3

# (mean) major allele will be 1, the rest will be 0
mean_allele_freq <- resmarkers_table %>%
  group_by(resmarker, AA) %>%
  summarise(mean_freq = mean(norm.reads.locus))%>%
  group_by(resmarker) %>%
  mutate(FEMcoded = ifelse(mean_freq == max(mean_freq), 1, 0))

#add 0.5 to bi/multiallelic samples
resmarkers_table$FEMcoded <- ifelse(resmarkers_table$n.alleles > 1, 0.5, NA)

#add 1 to samples with major allele only and 0 to samples with minor allele only
merged_table <- left_join(resmarkers_table, mean_allele_freq, by = c("resmarker", "AA"))

# Fill NA values in resmarker_table$FEMcoded with corresponding values from mean_allele_freq$FEMcoded
resmarkers_table$FEMcoded <- ifelse(is.na(resmarkers_table$FEMcoded), merged_table$FEMcoded.y, resmarkers_table$FEMcoded)

#final input file for FEM
FEMcoded <- resmarkers_table %>%
  group_by(SampleID, resmarker) %>%
  summarise(FEMcoded = unique(FEMcoded))

formatted_table <- FEMcoded %>%
  spread(resmarker, FEMcoded)

formatted_table <- as.data.frame(formatted_table)  

row.names(formatted_table) <- formatted_table$SampleID
formatted_table$SampleID <- NULL

formatted_table <- formatted_table %>%
  select(markers_to_phase)

write.csv(formatted_table, file = "FEMcoded_combined_mixture_controls_resmarker_table_16_17_21_22_27_manuel.csv")
write.csv(mean_allele_freq, file = "FEMcoded_CODE_combined_mixture_controls_resmarker_table_16_17_21_22_27_manuel.csv", row.names = F)
