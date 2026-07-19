# Data

## Coursework 1 — Portfolio construction

Weekly return data for 441 S&P 500 constituents (2000–2023) and a matching weekly risk-free rate series. Not included in this repository due to data licence restrictions.

Files used:
- `DA441_weekly_returns.csv` — 1,251 × 441 matrix of weekly stock returns
- `risk_free_rate_weekly.csv` — weekly risk-free rate in percentage form

100 stocks are selected from the 441 using `set.seed(7)` at the start of the script.

## Coursework 2 — Insurance claim severity

French motor third-party liability (MTPL) data, publicly available via the `CASdatasets` R package.

Files used:
- `freMTPL2freq.csv` — 677,991 policies with 12 risk features
- `freMTPL2sev.csv` — 26,444 individual claims

To obtain the data:
```r
install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/", type = "source")
library(CASdatasets)
data(freMTPL2freq)
data(freMTPL2sev)
```

1,000 policies are sampled from the 23,571 eligible single-claim policies using `set.seed(7)`.
