# Cybersickness Guideline
This is the repository for the code and data related to our meta-analysis of the washout period's impact on cybersickness experiments.
It includes the R code for meta-regression and the Python code for power analysis.
The results were published in our IEEE TVCG paper: "[Revenge of the Sick: A Meta-Analysis of Washout Periods in Cybersickness Research](https://ieeexplore.ieee.org/document/11569370)"

## Data
The effect sizes we collected for this paper is in `Results.csv`. 
`Results_delta_only.csv` contains only delta-SSQ data, which is a subset of the data in `Results.csv`.

## Mega-Regression
See `Meta.R` for details.

## Simulation Based Power Analysis
For between-subjects and within-subjects (with or without carryover) designs, we simulate their power given 30 participants by generating data and run corresponding t-tests on the data.
The power calculation and visualization is in `PowerCalculation.R`.

## Bibtex
```
@article{nie_revenge_2026,
 author = {Nie, Tongyu and Cantory, Ville and Hutton Pospick, Courtney and Zhang, Danhua and Tan, Haoyu and Adeniyi, Sam and Wu, Fei and Zielasko, Daniel and Adhanom, Isayas Berhe and Suma Rosenberg, Evan},
 doi = {10.1109/TVCG.2026.3704898},
 journal = {IEEE Transactions on Visualization and Computer Graphics},
 keywords = {journal},
 number = {9},
 pages = {1--20},
 title = {Revenge of the Sick: A Meta-Analysis of Washout Periods in Cybersickness Research},
 volume = {32},
 year = {2026}
}
```

