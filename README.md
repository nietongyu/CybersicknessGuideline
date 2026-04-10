# Cybersickness Guideline
This is the cybersickness guideline project's code.
It includes the R code for meta-regression and the Python code for power analysis

## Data
The effect sizes we collected for this paper is in `Results.csv`. 
`Results_delta_only.csv` contains only delta-SSQ data, which is a subset of the data in `Results.csv`.

## Mega-Regression
See `Meta.R` for details.

## Simulation Based Power Analysis
For between-subjects and within-subjects (with or without carryover) designs, we simulate their power given 30 participants by generating data and run corresponding t-tests on the data.
The power calculation and visualization is in `PowerCalculation.R`.
