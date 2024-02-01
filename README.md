# FapR vs FreqEstimationModel Phasing Benchmark

## Overview

This R script benchmarks the performance of FapR and FreqEstimationModel (FEM) for phasing haplotypes of *Plasmodium falciparum* (P. falciparum). The benchmarking process involves comparing the results of each method against expected control data, considering metrics such as accuracy, precision, recall, F1 score, root mean square error (RMSE), and mean absolute percentage error (MAPE) at different minor allele frequency (MAF) thresholds.

## Requirements

Make sure to install the required R packages before running the script:

```R
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("dplyr")
install.packages("readxl")
```

## Usage

1. Set the correct file paths for the expected control data, FEM output, and FapR output.

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

## Additional Analysis

The script further provides additional analyses, such as separating the benchmarking results for different parasitaemias and generating individual plots. The resulting plots are saved as images for further inspection.

```R
for (df in 1:length(result_data_FINAL_all_parasitaemias_separatedly)){
  plotMetricsGrid(result_data_FINAL_all_parasitaemias_separatedly[[df]], names(result_data_FINAL_all_parasitaemias_separatedly[df]), save_plot = TRUE)
}
```

## Metrics Comparison

The script also compares metrics between FEM and FapR, providing separate plots for each metric:

### FEM Metrics Plot
![FEM Metrics Plot](benchmark_FEM_metrics_plot.png)

### FAPR Metrics Plot
![FAPR Metrics Plot](benchmark_FAPR_metrics_plot.png)

## Notes

- Adjust the file paths and boolean values according to your specific dataset and benchmarking preferences.
- Ensure that the required packages are installed before running the script.
