
# Surrogate Assisted Coefficients Updating in Logistic Regression considering Partially Missing Covariates
![R Programming](https://img.shields.io/badge/R%20Programming-4.0.0%2B-blue)
![Logistic Regression](https://img.shields.io/badge/Logistic%20Regression-Coefficient%20Correction-green)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)


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

## The Six Methods

1. **Golden Method**:
   - Treats the validation dataset as fully observed, combining the main and validation datasets for estimation.
   - Produces the most accurate and unbiased estimates but assumes complete and accurate data for all covariates.

2. **Main Method**:
   - Uses only the main dataset for estimation.
   - Limited to fully observed covariates and does not leverage any information from partially observed data, leading to potential biases in coefficients.

3. **Naïve Method**:
   - Substitutes partially missing covariates with their surrogate covariates (e.g., proxies that are fully observed in the validation dataset).
   - Simplistic approach but prone to significant biases if the surrogate covariates poorly represent the missing data.

4. **Calibration Method**:
   - Uses a calibration model to adjust coefficients based on the relationships between missing covariates and surrogate covariates observed in the validation dataset.
   - Assumes a linear relationship between missing and surrogate covariates and is robust in scenarios where this assumption holds.

5. **Updating Method**:
   - Starts with a reduced logistic model that omits missing covariates.
   - Corrects the coefficients using an adjustment term derived from the main dataset, assuming a large sample size and full column rank of the design matrix.

6. **New Method**:
   - Combines the calibration and updating methods.
   - Utilizes the robustness of the updating method and the surrogate relationships from the calibration method to achieve more accurate and reliable estimations.


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
   result <- function_oh(data = "framingham_data_example.csv", method)
   print(result)
   ```

## Framingham Data

The Framingham Heart Study dataset is a subset of a well-known longitudinal study focusing on cardiovascular health. The example file (`framingham_data_example.csv`) includes:
- Demographic variables: Age, Sex, etc.
- Health metrics: BMI, Blood Pressure, etc.
- Outcomes: Myocardial infarction, Stroke, etc.

For more information, see the included `framingham_data_explanation.pdf`.

---

## Acknowledgments

- Publicly accessible data for research purposes can be requested from the Framingham Heart Study: [https://www.framinghamheartstudy.org/fhs-for-researchers/](https://www.framinghamheartstudy.org/fhs-for-researchers/).

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

