
# A Surrogate-calibrated Updating Method for Logistic Regression with Missing Covariates
![R Programming](https://img.shields.io/badge/R%20Programming-4.0.0%2B-blue)
![Logistic Regression](https://img.shields.io/badge/Logistic%20Regression-Coefficient%20Correction-green)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)


This repository provides an implementation of a framework for correcting logistic regression coefficients in the presence of missing covariates. The core implementation is encapsulated in the `function_oh` script, which allows users to apply one of six adjustment methods to their dataset.

This repository accompanies the paper:

**A Surrogate-Calibrated Updating Method for Logistic Regression with Missing Covariates**

## Overview

In logistic regression, coefficient estimates can be biased when some covariates are missing for part of the sample. This project implements several correction methods for this setting, including the proposed **SCU** method, which combines calibration and updating to improve estimation.

The code is designed for settings with:

- **`X1`**: fully observed covariates
- **`X2`**: partially missing covariates
- **`Z`**: surrogate covariates related to `X2`
- **`D`**: binary outcome

## Files in the Repository

- **`function_oh.R`**: Main R script implementing six correction methods for logistic regression.
- **`framingham_data_example.csv`**: Example dataset derived from the Framingham Heart Study.
- **`framingham_data_explanation.pdf`**: Detailed documentation about the Framingham dataset.

## Requirements

- **R version**: 4.0.0 or later

Required packages:

- `dplyr`
- `tidyr`
- `MASS`
- `ggplot2`
- `forcats`
- `reshape2`
- `Matrix`
- `Rcpp`
- `RcppArmadillo`

Install missing packages with:

```R
install.packages(c(
  "dplyr", "tidyr", "MASS", "ggplot2",
  "forcats", "reshape2", "Matrix", "Rcpp", "RcppArmadillo"
))
```

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

6. **SCU(Surrogate-calibrated Upating) Method**:
   - Combines the calibration and updating methods.
   - Utilizes the robustness of the updating method and the surrogate relationships from the calibration method to achieve more accurate and reliable estimations.


## Usage

1. **Setup**:
   - Place function_oh.R in your working directory and prepare the input objects X1, X2, Z, and D.

2. **Load the Function**:
   ```R
   source("function_oh.R")
   ```

3. **Apply a Correction Method**:
   - Example: Run SCU method on the Framingham dataset.
     
   ```R
    result <- function_oh(
      X1 = X1,
      X2 = X2,
      Z  = Z,
      D  = D,
      method = "new"
    )
    
    print(result)
   ```
---

## Output

The function returns a data frame containing inferential results for each coefficient, including:

- estimate
- standard error
- p-value
- odds ratio
- confidence interval

---

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

## Citation

Oh J, Shin YE. A Surrogate-Calibrated Updating Method for Logistic Regression With Missing Covariates. Stat Med. 2026 Mar;45(6-7):e70489. doi: 10.1002/sim.70489. PMID: 41830128; PMCID: PMC12988319.

---
## License

This project is licensed under the MIT License. See the LICENSE file for details.

