library(FreqEstimationModel)
library(dplyr)
library(reshape2)

file <- c("FEMcoded_CONTROLS_ALL.csv")

#load data
control_data <- read.csv(file, row.names = 1) # Load data 
control_data
#colnames(control_data) <- paste0("locus_", colnames(control_data))

control_data <- mutate_all(control_data, as.numeric)

#remove rows with NAs
control_data <- control_data[complete.cases(control_data), ]

control_data_list<-c()
control_data_list$Data <- control_data


thinning_interval <- 1 # Number of iterations per chain that are not saved 
no_traces_preburnin <- 10000 # For more traces but manageable pdf plots, don't exceed 10k and increase thinning interval instead
no_mcmc_chains <- 3 # Number of MCMC chains to run
parallel <- FALSE # Set to true if running code in parallel (this option isn't active yet)
NGS <- FALSE # Set to true if data are in NGS format (this option isn't active yet)
log_like_zero <- FALSE # QC check: sets log(p(yi|ai)) to zero (only impacts NGS)
mcmc_variable_list <- list(no_mcmc_chains = no_mcmc_chains,
                           no_traces_preburnin = no_traces_preburnin,
                           thinning_interval = thinning_interval,
                           NGS = NGS,
                           log_like_zero = log_like_zero)


if(!NGS | log_like_zero){ # Set to true to keep partial observations rather than discard them 
  augment_missing_data <- TRUE 
} else {
  augment_missing_data <- FALSE
}


if (log_like_zero) { # When log_like_zero == TRUE, code should return the prior
  moi_prior <- 'Uniform' # Use a prior that is easy to eye-ball
} else {
  moi_prior <- 'Poisson' # Choose between 'Uniform', 'Poisson', 'Geometric' or 'nBinomial'
}

moi_max <- 5 # Maximum MOI regarded as possible by the model (I haven't tested beyond 20)
moi_hyperparameter <- 2 # Specify population-average MOI (parameter of the prior on the MOI)
moi_size_hyperparameter <- 0.5 # Only applies if moi_prior == 'nBinomial' (hyperparameter for the prior on the MOI if the prior is the negative Binomial)
moi_prior_min2 <- NULL # Specify the lower bound for the MOI per individual
moi_initial <- NULL # If null, the initial vector of MOIs is set internally, otherwise set to input moi_initial
moi_list <- list(moi_hyperparameter = moi_hyperparameter,
                 moi_size_hyperparameter = moi_size_hyperparameter,
                 moi_prior = moi_prior,
                 moi_max = moi_max,
                 moi_prior_min2 = moi_prior_min2,
                 moi_initial = moi_initial)

########################################################### modify preprocess_data function 'cause it's not working as intended
preprocess_data_TEST<-function (data_summary, log_like_zero, NGS, augment_missing_data, 
          moi_prior_min2) 
{
  nlocus <- length(grep("locus", colnames(control_data_list$Data)))
  if (nlocus > 0) {
    no_marker_loci <- nlocus
  }
  else {
    no_marker_loci <- ncol(data_summary$Data)
  }
  if (log_like_zero) {
    control_data_list$Data[, 1:no_marker_loci] <- 99
  }
  if (NGS) {
    yij <- data_summary$yij[, 1:no_marker_loci, drop = FALSE]
    zij <- data_summary$zij[, 1:no_marker_loci, drop = FALSE]
  }
  else {
    yij <- NA
    zij <- NA
  }
  raw_data_pre_augment_step <- control_data_list$Data[, 1:no_marker_loci, 
                                                 drop = FALSE]
  ind_partial_initial <- apply(raw_data_pre_augment_step, 1, 
                               function(X) {
                                 99 %in% X
                               })
  if (!augment_missing_data & sum(ind_partial_initial) > 0) {
    raw_data <- raw_data_pre_augment_step[!ind_partial_initial, 
                                          , drop = FALSE]
  }
  else {
    raw_data <- raw_data_pre_augment_step
  }
  datasampleID <- rownames(raw_data)
  markerID <- colnames(raw_data)
  ind_partial <- apply(raw_data, 1, function(X) {
    99 %in% X
  })
  ind_mixed <- apply(raw_data, 1, function(X) {
    0.5 %in% X
  })
  ind_non_mixed <- as.logical((!ind_partial) * (!ind_mixed))
  if (is.null(moi_prior_min2)) {
    y_no_mxed <- datasampleID[!ind_mixed]
    y_mxed <- datasampleID[ind_mixed]
  }
  else {
    y_mxed <- unique(c(datasampleID[ind_mixed], moi_prior_min2))
    y_no_mxed <- datasampleID[match(datasampleID, y_mxed, 
                                    nomatch = 0L) == 0]
  }
  y_pure_nomissing <- datasampleID[ind_non_mixed]
  y_mxed_andor_missing <- datasampleID[!ind_non_mixed]
  no_total <- nrow(raw_data)
  masked_data <- raw_data
  masked_data[masked_data == 99] <- 0.5
  masked_data_unique <- unique(masked_data)
  masked_data_character <- apply(masked_data, 1, paste, sep = "", 
                                 collapse = ",")
  masked_data_character_unique <- unique(masked_data_character)
  no_masked_data_character_unique <- length(masked_data_character_unique)
  comp_genotypes <- c()
  for (i in 1:no_masked_data_character_unique) {
    x <- compatible_genotypes_TEST(masked_data_unique[i, ])
    comp_genotypes <- rbind(comp_genotypes, x)
  }
  comp_genotypes <- unique(comp_genotypes)
  comp_genotypes_character <- apply(comp_genotypes, 1, paste, 
                                    collapse = "")
  dimnames(comp_genotypes) <- list(comp_genotypes_character, 
                                   markerID)
  no_haplotypes <- nrow(comp_genotypes)
  frequency_truncated <- matrix(0, nrow = no_masked_data_character_unique, 
                                ncol = no_haplotypes, dimnames = list(masked_data_character_unique, 
                                                                      comp_genotypes_character))
  for (i in 1:no_masked_data_character_unique) {
    frequency_truncated[i, compatible_genotypes_TEST(unique(masked_data)[i, 
    ], as_numeric = FALSE)] <- 1
  }
  neg_log_sum_frequency_truncated <- -log(rowSums(frequency_truncated))
  processed_data_list <- list(no_marker_loci = no_marker_loci, 
                              yij = yij, zij = zij, raw_data = raw_data, datasampleID = datasampleID, 
                              markerID = markerID, y_no_mxed = y_no_mxed, y_mxed = y_mxed, 
                              y_pure_nomissing = y_pure_nomissing, y_mxed_andor_missing = y_mxed_andor_missing, 
                              no_total = no_total, masked_data_character = masked_data_character, 
                              comp_genotypes = comp_genotypes, comp_genotypes_character = comp_genotypes_character, 
                              no_haplotypes = no_haplotypes, frequency_truncated = frequency_truncated, 
                              neg_log_sum_frequency_truncated = neg_log_sum_frequency_truncated)
  return(processed_data_list)
}


compatible_genotypes_TEST<-function (observation, as_numeric = TRUE) 
{
  L <- length(observation)
  observation<-as.numeric(observation) ## ADDED THIS!!! WITHOUT THIS IT DOESN'T WORK (PUSH TO AIMEE's GH)
  if (!is.list(observation)) {
    if (!all(unique(observation) %in% c(0, 1, 0.5))) {
      stop("invalid allelic observations")
    }
    else {
      yl = vector("list", length = L)
      for (j in 1:L) {
        if (observation[j] == 0.5) {
          yl[[j]] = 0:1
        }
        else {
          yl[[j]] = observation[j]
        }
      }
    }
  }
  else {
    yl <- observation
  }
  hapmat = as.matrix(expand.grid(yl))
  genotypes <- apply(hapmat, 1, paste, sep = "", collapse = "")
  if (as_numeric == TRUE) {
    return(hapmat)
  }
  else {
    return(genotypes)
  }
}
############################################################33

processed_data_list <- preprocess_data_TEST(control_data_list,
                                       log_like_zero,
                                       NGS,
                                       augment_missing_data,
                                       moi_prior_min2)

processed_data_list

#############################################################3




frequency_hyperparameter <-rep(1, processed_data_list$no_haplotypes) # The Parameter vector for the Dirichlet prior on the frequency vector
frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to input external_frequency_initial
frequency_list <- list(frequency_hyperparameter = frequency_hyperparameter,
                       frequency_initial = frequency_initial)


seed <- 1 # For reproducibility
set.seed(seed)

# Run MCMC
results <- mcmc_sampling_parallel(processed_data_list,
                                  moi_list,
                                  frequency_list,
                                  mcmc_variable_list,
                                  cores_max = 8)


# Install, load, and attach packages that are helpful for unpacking example results
if(!require("plyr")) install.packages("plyr")
if(!require("coda")) install.packages("coda")
if(!require("abind")) install.packages("abind")



burnin <- 1:(0.5*mcmc_variable_list$no_traces_preburnin) # remove first half
if(mcmc_variable_list$no_mcmc_chains > 1){
  alply_genotype_freq_store_chains_burnin <- plyr::alply(results$genotype_freq_store_chains[-burnin,,],3)
}else{
  alply_genotype_freq_store_chains_burnin <- results$genotype_freq_store_chains[-burnin,,]
}

# Use mcmc.list to represent frequency chains across parallel runs (see ?coda::mcmc.list)
mcmc_frequency_chains <- coda::mcmc.list(lapply(alply_genotype_freq_store_chains_burnin, # 3 for splitting by third dimension (nothing to do with no. of chains)
                                                coda::mcmc,
                                                start = (max(burnin)+1)*mcmc_variable_list$thinning_interval,
                                                end = mcmc_variable_list$no_traces_preburnin*mcmc_variable_list$thinning_interval,
                                                thin = mcmc_variable_list$thinning_interval))

mcmc_As <- abind::abind(plyr::alply(results$genotype_count_store_chains[-burnin,,,],4), along = 1) # Haplotype counts for all chains excluding burnin
mcmc_mois <- apply(mcmc_As, c(1,2), sum)



pop_freq <- cbind(mean = summary(mcmc_frequency_chains)$statistics[,"Mean"],
                  median = summary(mcmc_frequency_chains)$quantiles[,3],
                  "CI2.5%" = summary(mcmc_frequency_chains)$quantiles[,1],
                  "CI97.5%" = CI_upper<-summary(mcmc_frequency_chains)$quantiles[,5])


compute_mode <- function(x) {z <- table(x); names(z)[which.max(z)]}
qprobs <- c(0.025, 0.5, 0.975)
inf_MOI <- data.frame(colMeans(mcmc_mois),
                      as.numeric(apply(mcmc_mois, 2, compute_mode)),
                      t(apply(mcmc_mois, 2, quantile, probs = qprobs)))
colnames(inf_MOI) <- c("mean","mode","median","CI2.5%","CI97.5%")


mode_hap_count_chr <- apply(mcmc_As, 2, function(x) {compute_mode(apply(x, 1, paste, collapse = ""))})
inf_As <- t(sapply(mode_hap_count_chr, function(x) as.numeric(strsplit(x, split = "")[[1]])))
colnames(inf_As) <- dimnames(mcmc_As)[[3]]


inf_freq <- t(apply(mcmc_As, 2, function(x) colMeans(x / rowSums(x))))


pop_prev <- 1-(1-pop_freq)^median(mcmc_mois)


inf_prev <- t(apply(mcmc_As, 2, function(x) colMeans(x > 0)))
pop_prev2 <- colMeans(inf_prev)


# Marginal allele frequencies 
haplotypes_str <- rownames(pop_freq) 
haplotypes_chr <- do.call(rbind, strsplit(haplotypes_str, split = ""))
ind <- apply(haplotypes_chr, 1, function(x) x == "1") # Compute the frequency of the default allele
allele_freq <- sapply(1:nrow(ind), function(i) sum(pop_freq[,"mean"][ind[i,]]))
biallele_freq <- rbind(allele_freq, 1-allele_freq)
colnames(biallele_freq) <- processed_data_list$markerID
rownames(biallele_freq) <- c("0", "1")

# pop_freq 
# inf_MOI
# inf_As
# inf_freq
# pop_prev 
# inf_prev
# biallele_freq

################################################################################################
# translate encodings into actual genotypes, use inf_freq
markers_code <- read.csv("FEMcoded_CODE_CONTROLS_ALL.csv")
inf_freq

#dimensions
markers_order<-colnames(control_data)
col_names <- colnames(inf_freq)

# Create a matrix with zeros
matrix_data <- matrix(0, nrow = length(col_names), ncol = nchar(col_names[1]))

# Fill in the matrix based on the occurrences of each digit in each position
for (i in 1:length(col_names)) {
  matrix_data[i, ] <- as.numeric(strsplit(col_names[i], NULL)[[1]])
}

#binary codes
df_genotypes_binary <- data.frame(matrix_data)
colnames(df_genotypes_binary) <- markers_order


#account for multillelic loci 
code_results_list <- list()

for (i in markers_order) {
  subset <- markers_code[markers_code$resmarker == i,]
  
  if (sum(subset$FEMcoded == 0) > 1) {
    AA_new <- paste(subset$AA[subset$FEMcoded == 0], collapse = "_")
    AA_new <- paste0("[", AA_new, "]")
    mean_freq <- paste(subset$mean_freq[subset$FEMcoded == 0], collapse = "_")
    
    new_line <- as.data.frame(t(c(i, AA_new, mean_freq, 0)))
    colnames(new_line) <- colnames(subset)
    subset <- subset[!subset$FEMcoded == 0,]
    subset <- rbind(subset, new_line)
    
    print(subset)
  } else {
    print(subset)
  }
  
  # Add the subset to the list
  code_results_list[[i]] <- subset
}

# Concatenate all data frames in the list into a single data frame
code_results <- do.call(rbind, code_results_list) # LISTO EL CÓDIGO!!
rownames(code_results)<-NULL


#decode binary code to genotype
df_genotypes_binary_ <- df_genotypes_binary
# Iterate over rows in code_results
for (i in 1:nrow(code_results)) {
  resmarker <- code_results$resmarker[i]
  AA <- code_results$AA[i]
  FEMcoded <- code_results$FEMcoded[i]
  
  # Find corresponding column in df_genotypes_binary
  col_index <- which(colnames(df_genotypes_binary_) == resmarker)
  
  # Replace values in df_genotypes_binary based on FEMcoded
  df_genotypes_binary_[, col_index] <- ifelse(df_genotypes_binary_[, col_index] == FEMcoded, AA, df_genotypes_binary_[, col_index])
}

#concat genotypes and binary code
genos <- apply(df_genotypes_binary_, 1, function(row) paste(row, collapse = "_"))
FEMcode <- apply(df_genotypes_binary, 1, function(row) paste(row, collapse = ""))
FINAL_DECODED_GENOTYPES <- as.data.frame(cbind(genotypes=genos, FEMcode=as.numeric(FEMcode)))


# Melt the data frame
melted_inf_freq <- melt(as.matrix(inf_freq))
colnames(melted_inf_freq) <- c("SampleID", "FEMcode", "freq")

melted_inf_prob <- melt(as.matrix(inf_As)) #keep only most likely haplos

melted_inf_freq$value <- melted_inf_prob$value

melted_inf_freq <- melted_inf_freq %>%
  arrange(SampleID) %>%
  filter(value != 0)

#decode
merged_data <- merge(melted_inf_freq, FINAL_DECODED_GENOTYPES, by = "FEMcode", all.x = TRUE)

FINAL_merged_data <- merged_data %>%
  arrange(SampleID) %>%
  select("SampleID", "genotypes", "freq")

write.csv(FINAL_merged_data, "controls_phased_haplotypes_FEM.csv", row.names = F)
