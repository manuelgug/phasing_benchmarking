
library(dplyr)

pool1B2 <- read.csv("truth_from_monoclonals_1B_2.csv")

#subset amplicons of interest
amps_names <- c("Pf3D7_04_v3-748105-748359-1B", "Pf3D7_04_v3-748374-748611-1B") #dhfr
amps_names <- c(amps_names, "Pf3D7_08_v3-549583-549807-1B", "Pf3D7_08_v3-549960-550215-1B") #dhps

condition <- pool1B2$locus %in% amps_names
subsetted_pool1B2 <- subset(pool1B2, condition)

subsetted_pool1B2$resmarkers<-NA
subsetted_pool1B2[subsetted_pool1B2$locus=="Pf3D7_04_v3-748105-748359-1B",][["resmarkers"]]<-"dhfr_51_59"
subsetted_pool1B2[subsetted_pool1B2$locus=="Pf3D7_04_v3-748374-748611-1B",][["resmarkers"]]<-"dhfr_108"
subsetted_pool1B2[subsetted_pool1B2$locus=="Pf3D7_08_v3-549583-549807-1B",][["resmarkers"]]<-"dhps_437"
subsetted_pool1B2[subsetted_pool1B2$locus=="Pf3D7_08_v3-549960-550215-1B",][["resmarkers"]]<-"dhps_540"

subsetted_pool1B2$GENOTYPE<-NA
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="36+10N",][["GENOTYPE"]]<-"NC"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="36+10N110T133C",][["GENOTYPE"]]<-"IR"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="3A18+9N72+9N",][["GENOTYPE"]]<-"N"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="3A18+9N72+9N170T",][["GENOTYPE"]]<-"N"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="71G75C",][["GENOTYPE"]]<-"A"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="75C",][["GENOTYPE"]]<-"A"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="72T",][["GENOTYPE"]]<-"G"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar==".",][["GENOTYPE"]]<-"G"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="38+8N",][["GENOTYPE"]]<-"K"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="18+9N72+9N",][["GENOTYPE"]]<-"S"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="4G38+8N",][["GENOTYPE"]]<-"E"
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="5T36+10N",][["GENOTYPE"]]<-"NC" #DUDA.....
subsetted_pool1B2[subsetted_pool1B2$pseudo_cigar=="3C18+9N72+9N",][["GENOTYPE"]]<-"T" #DUDA.....
subsetted_pool1B2[is.na(subsetted_pool1B2$GENOTYPE),][["GENOTYPE"]]<-"?"

# build genotypes: dhps then dhfr
df <- subsetted_pool1B2 %>%
  group_by(Strain) %>%
  summarise(Full_Genotype = paste(GENOTYPE, collapse = "")) %>%
  ungroup()

df
