# Monte Carlo Permutation Tests and Nearest Neighbor Methods in R

This project was developed as part of a graduate-level course on **Monte Carlo Methods**. It explores how nonparametric permutation-based hypothesis testing can be implemented using Monte Carlo simulation techniques. The focus is on assessing whether two samples are drawn from the same distribution using custom and built-in statistical methods in R.

---

## Objectives

- Implement Monte Carlo-based permutation tests to evaluate distributional equality.
- Explore the power and interpretation of:
  - Kolmogorov-Smirnov test
  - Cramér-von Mises test
  - r-th Nearest Neighbor test
- Simulate the null distribution of test statistics via repeated resampling (up to 1000 permutations).
- Visualize permutation test distributions and p-values.

---

## Methods

### 1. **Kolmogorov-Smirnov Permutation Test**
- Univariate test on chick weight distributions (casein vs sunflower feed).
- Monte Carlo simulation with 999 permutations.
- Test statistic from original sample compared to null distribution.

### 2. **Cramér-von Mises Test**
- Implemented via `twosamples` R package.
- Measures squared difference between empirical CDFs.
- Applied to the same `chickwts` data with 1000 permutations.

### 3. **r-th Nearest Neighbor Test (Custom Implementation)**
- Designed to test for multivariate distributional equality using nearest neighbor structure.
- Applied to the `iris` dataset to compare *setosa* vs *virginica*.
- Computed via the `boot` package using user-defined function and `yaImpute` for efficient neighbor lookup.

---

## Datasets Used

- `chickwts`: Built-in R dataset comparing chicken weights by feed type.
- `iris`: Multivariate measurements of iris flower species.

---

## Packages Used

- `boot`: for permutation resampling
- `yaImpute`: for nearest neighbor indexing
- `twosamples`: for Cramér-von Mises test
- Base R for data wrangling and visualization

---

## Visualizations

Each test includes a histogram of the simulated (null) distribution of the test statistic with the observed value marked. These plots help visualize where the observed statistic falls under the null hypothesis.

---

## File Overview

- `week7_analysis.R`: Full R script implementing all tests and visualizations
- `Week7_Rorick.pdf`: Written summary of code, outputs, and conclusions (course submission)

---

## Learning Takeaways

- Applied Monte Carlo simulation to derive empirical p-values through permutation.
- Gained practical experience with resampling techniques for nonparametric testing.
- Developed custom statistical functions and interpreted high-dimensional test outcomes.
