# Advanced Predictive Analytics: Shrinkage Estimation for Portfolio Construction and Insurance Claim Severity

**MSc Quantitative Finance — SMM069 Advanced Predictive Analytics**
Bayes Business School, City St George's, University of London
Submitted: July 2026 | Grade: pending

Group project by Camilla Ayres, Ramprashanth Elangovan, Mahad Khan, and Raghav Iyengar.

---

## Overview

This project applies shrinkage estimation methods to two problems in quantitative finance and insurance: mean-vector shrinkage for portfolio construction (Coursework 1) and shrinkage generalised linear models for insurance claim severity estimation (Coursework 2). All implementation is in R.

The full write-up is in [`report/SMM069_Group7_report.pdf`](report/SMM069_Group7_report.pdf).

---

## Data

**Coursework 1:** 441 S&P 500 stocks, 1,251 weekly return observations spanning 2000–2023, plus a weekly risk-free rate series. 100 stocks selected using `set.seed(7)`. Data not included due to licence restrictions.

**Coursework 2:** French motor third-party liability data (`freMTPL2freq.csv` and `freMTPL2sev.csv`), publicly available via the `CASdatasets` R package. 1,000 policies sampled from 23,571 eligible single-claim policies using `set.seed(7)`.

---

## Coursework 1 — Shrinkage Estimation for Portfolio Construction

### Setup

100 stocks drawn at random from 441 S&P 500 constituents (2000–2023, fixed by `set.seed(7)`). Weekly excess returns computed by subtracting the risk-free rate. Rolling-window backtest: 156-week training window, 4-week test window, 4-week step, producing 273 out-of-sample windows (2003–2023).

Three portfolios constructed per window: Mean-Variance (γ = 1), Tangency, and Equal-Weight benchmark. Covariance estimated throughout using `cov1Para` (Ledoit–Wolf single-parameter shrinkage toward scaled identity), which guarantees a well-conditioned inverse with n = 156 observations for a 100 × 100 matrix.

### Mean estimators compared

Eight mean-vector estimators evaluated: Sample mean, St-MSh, D-MSh, O-LSh, T-LSh, Wang, BOP, and Jorion. CW excluded because n > p makes it collapse to a scalar multiple of the sample mean. Performance scored by annualised Sharpe ratio across three periods: full sample, Financial Crisis (Oct 2008 – Aug 2011), and Pandemic (Feb 2020 – Dec 2023).

### Key findings

No optimised portfolio beats the equal-weight benchmark over the full sample (EWP: 1.54 vs best optimised: 1.14 for Wang MV), confirming severe estimation error when estimating 100 means from 156 weekly observations. St-MSh reproduces the sample mean exactly in the Sharpe ratio (shrinkage scalar cancels). BOP is the weakest method overall.

Regime behaviour differs sharply: during the Financial Crisis, grand-mean shrinkage estimators (Wang, O-LSh) hold up best; during the Pandemic, sample mean and Jorion deliver higher but more volatile returns. The Tangency and MV portfolios diverge during stress because the Tangency normalisation flips sign when 1ᵀΣ⁻¹µ is negative.

---

## Coursework 2 — Shrinkage GLMs for Insurance Claim Severity

### Setup

1,000 French MTPL policies with exactly one positive claim (`set.seed(7)`). Target: per-policy severity (€4.1 to €114,760, mean €2,051, strongly right-skewed). Gamma GLM with log link. Feature engineering: bounded log transforms for skewed continuous covariates, natural cubic splines for DrivAge and VehAge (df = 3), credibility-weighted LOO target encoding for high-cardinality categoricals (Region, Area, VehBrand, k = 20). 5-fold cross-validation repeated 100 times (500 train/test runs per model).

### Models compared

Nine models: unpenalised GLM2 (baseline), Ridge (α = 0), Elastic Net (α = 0.5), and six shrinkage GLMs (St, DSh, SR, GSR, LW, QIS). Ridge and Elastic Net select λ by inner 5-fold CV. Metrics: RMSE, MAE, winsorised RMSE/MAE at 95% and 99%, Gamma deviance, Gini coefficient, actual-to-expected ratio, computation time, and solver iterations.

### Key findings

Ridge is the best overall performer: lowest RMSE (5,971), lowest RMSE₉₅ (0.918× GLM2), lowest deviance (0.940× GLM2), with positive Gini and AE of 1.05. DSh is the efficiency-optimal choice: recovers ~99% of Ridge's RMSE gain at 4.9× GLM2 runtime versus Ridge's 169×, and converges in 0.40× the GLM2 iterations. LW achieves the highest Gini but is poorly calibrated (AE = 0.886, ~11% aggregate over-shrinkage).

In-sample fit of the selected Ridge model on all 1,000 policies: aggregate AE = 1.038, mean actual/fitted = €2,051/€1,977. Fitted values occupy a narrow band (expected for a shrinkage estimator), so the model over-predicts typical small claims and cannot reach extreme tails. Residual error is driven by a small number of catastrophic claims rather than systematic bias.

---

## Repository Structure

```
├── README.md
├── report/
│   └── SMM069_Group7_report.pdf
├── code/
│   └── smm069_shrinkage_portfolio_and_glm.R
└── data/
    └── README.md
```

---

## Requirements

```r
install.packages(c(
  "tidyverse", "PerformanceAnalytics", "nlshrink", "glmnet",
  "CASdatasets", "splines", "Matrix"
))
```

Set `set.seed(7)` is used throughout for reproducibility. The script is self-contained: run it top to bottom with the data files in the working directory.

---

## Key References

- Ledoit, O. and Wolf, M. (2004). A well-conditioned estimator for large-dimensional covariance matrices. *Journal of Multivariate Analysis*, 88(2).
- Jorion, P. (1986). Bayes–Stein estimation for portfolio analysis. *Journal of Financial and Quantitative Analysis*, 21(3).
- Wang, C. et al. (2020). Shrinkage estimation of high-dimensional mean vectors. *Journal of Multivariate Analysis*.
- Ohlsson, E. and Johansson, B. (2010). *Non-Life Insurance Pricing with Generalised Linear Models*. Springer.
