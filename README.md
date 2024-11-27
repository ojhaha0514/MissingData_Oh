
# Surrogate Assisted Coefficients Updating in Logistic Regression considering Partially Missing Covariates

This repository provides an implementation of a framework for correcting logistic regression coefficients in the presence of missing covariates. The core implementation is encapsulated in the `function_oh` script, which allows users to apply one of six adjustment methods to their dataset.

## Overview

The methods included address scenarios where logistic regression suffers from biases due to missing covariates. Key features include:
- **Six Correction Methods**: Options for various data setups and biases.
- **Framingham Data Example**: Demonstrates the methods using the well-documented Framingham Heart Study dataset.

## Files in the Repository

- **`function_oh.R`**: Main R script implementing six correction methods for logistic regression.
- **`framingham_data_example.csv`**: Example dataset derived from the Framingham Heart Study.
- **`framingham_data_explanation.pdf`**: Detailed documentation about the Framingham dataset.

## Prerequisites

- **R Version**: 4.0.0 or later
- Required R packages:
  - `tidyverse`
  - `dplyr`
  - `ggplot2`

## Usage

1. **Setup**:
   - Place `function_oh.R` and `framingham_data_example.csv` in the working directory.

2. **Load the Function**:
   ```R
   source("function_oh.R")
   ```

3. **Apply a Correction Method**:
   - Example: Run method 1 on the Framingham dataset.
   ```R
   result <- function_oh(data = "framingham_data_example.csv", method = 1)
   print(result)
   ```

4. **Visualization**:
   - Use the `ggplot2` package for visualizing results if necessary.

## Framingham Data

The Framingham Heart Study dataset is a subset of a well-known longitudinal study focusing on cardiovascular health. The example file (`framingham_data_example.csv`) includes:
- Demographic variables: Age, Sex, etc.
- Health metrics: BMI, Blood Pressure, etc.
- Outcomes: Myocardial infarction, Stroke, etc.

For more information, see the included `framingham_data_explanation.pdf`.

---

## Citation

If you use this code or data in your work, please cite:
**Jooha Oh (2024)**: "Comparison of Correction Methods for Logistic Regression Coefficients in the Presence of Missing Covariates."
