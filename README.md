# FapR vs FreqEstimationModel Phasing Benchmark

## Overview

This R script benchmarks the performance of FapR and FreqEstimationModel for phasing haplotypes of *Plasmodium falciparum* (P. falciparum). The benchmarking process involves comparing the results of each method against expected control data, considering metrics such as accuracy, precision, recall, F1 score, root mean square error (RMSE), and mean absolute percentage error (MAPE) at different minor allele frequency (MAF) thresholds.

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

The script also compares metrics between FreqEstimationModel and FapR based on parasitaemia, utilizing various minimum allele frequency (MAF) cutoffs. Accuracy, precision, recall, and F1 score assess the phasing of haplotypes (presence/absence), whereas root mean square error (RMSE) and mean absolute percentage error (MAPE) gauge the error in frequency relative to the expected values*, the later taking into account the sample size.

### Freq EstimationModel Metrics Plot
![EstimationModel Metrics Plot](https://github.com/manuelgug/phasing_benchmarking/blob/main/results/benchmark_FEM_metrics_plot.png)

### FAPR Metrics Plot
![FAPR Metrics Plot](https://github.com/manuelgug/phasing_benchmarking/blob/main/results/benchmark_FAPR_metrics_plot.png)

### Overall comparison across parasitaemias
![Overall comparison across parasitaemias](https://github.com/manuelgug/phasing_benchmarking/blob/main/results/benchmarking_parasitaemiaALL_parasitaemias.png)

## Notes
*This refers theoretical proportions of strains in the mixes. In practice, it tends to differ.
