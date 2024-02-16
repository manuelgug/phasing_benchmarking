# FapR vs FreqEstimationModel Phasing Benchmarking

## Overview

This R script benchmarks the performance of [FapR](https://github.com/manuelgug/FapR) and [Aimee Taylor's FreqEstimationModel](https://github.com/aimeertaylor/FreqEstimationModel/) for phasing 7 loci across the `dhfr` and `dhps` genes from *Plasmodium falciparum*. The benchmarking process involves comparing the results of each method against an expected control dataset.

### Steps

The benchmarking is done across minor allele frequency (MAF) thresholds ranging from 0 to 0.4, ivolving 2 main steps:

1. Verify the accurate identification of haplotypes (presence or absence) and assess the outcomes using metrics like accuracy, precision, recall, and F1 score, derived from:

   - __TP (True Positives)__: The algorithm correctly identifies a haplotype that exists in the expected data. 
  
   - __FP (False Positives)__: The algorithm incorrectly identifies a haplotype that does not exist in the expected data.
  
   - __FN (False Negatives)__: The algorithm fails to identify a haplotype that actually exists in the expected data.

2. Ensure the accuracy of sample frequency estimates by computing both the Root Mean Square Error (RMSE) and the Mean Absolute Percentage Error (MAPE), which considers sample size, against the expected frequencies.


### Expected Dataset Description

The [expected control dataset](https://github.com/manuelgug/phasing_benchmarking/blob/main/inputs/controls_EXPECTED.xlsx) is a 4-column Excel file containing 100 single and 95 mixed controls including `3D7`, `HB3`, `DD2`, `D10`, `D6`, `FCR3`, `U659`, `V1S` and `W2` strains with different genotypes, frequencies, and parasitaemia levels:

- __SampleID__: Unique identifier for each control.
- __Genotypes__: A string representing the genotypes of the sample in the followig order: `dhps_431`, `dhps_437`, `dhps_540`, `dhps_581`, `dhfr_51`, `dhfr_59` and `dhfr_108`.
- __Freq__: The frequency of genotypes within its sample. (This refers theoretical proportions of strains in the mixes; in practice, it differs).
- __Parasitaemia__: Parasitaemia level associated with the sample in parasites/uL.


## Usage

1. Set the correct file paths for the expected control data, FreqEstimationModel output, and FapR output.

```R
expected <- readxl::read_xlsx("inputs/controls_EXPECTED.xlsx")
FEM <- read.csv("run_FreqEstimationModel/controls_phased_haplotypes_FEM.csv")
FAPR <- read.csv("run_FapR/controls_fapr_phased_haplos.csv")
```

2. Specify whether to benchmark only biallelic samples (`biallelic_benchmark = TRUE`) or all samples, including poliallelic ones (`biallelic_benchmark = FALSE`).

```R
biallelic_benchmark <- FALSE # Change to TRUE if benchmarking only biallelic samples
```

3. Run the main benchmarking function:

```R
result_data_FINAL_all_parasitaemias <- compareMethods(FAPR_mixes, FEM_mixes, expected_mixes)
```

4. Run the plotting function for visualizing metrics:

```R
plotMetricsGrid(result_data_FINAL_all_parasitaemias, "All_Parasitaemias", save_plot = TRUE)
```

## Results

### FreqEstimationModel Metrics Plot
![EstimationModel Metrics Plot](https://github.com/manuelgug/phasing_benchmarking/blob/main/results/benchmark_FEM_metrics_plot.png)

### FapR Metrics Plot
![FAPR Metrics Plot](https://github.com/manuelgug/phasing_benchmarking/blob/main/results/benchmark_FAPR_metrics_plot.png)

### Overall comparison across parasitaemias
![Overall comparison across parasitaemias](https://github.com/manuelgug/phasing_benchmarking/blob/main/results/benchmarking_parasitaemiaALL_parasitaemias.png)

