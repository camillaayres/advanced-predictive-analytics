############################################################
# SMM069 - Advanced Predictive Analytics
# Coursework 1: Shrinkage Estimation for Portfolio Construction
# Tasks A1-A5  |  GROUP 7
############################################################

rm(list = ls())
options(scipen = 999)

## ----------------------------------------------------------
## 0. Parameters
## ----------------------------------------------------------
GROUP_ID       <- 7      # group number -> seed and 100-stock subset
TRAIN_LEN      <- 52 * 3 # 156 weeks (3 years)
TEST_LEN       <- 4      # 4 weeks (1 month)
STEP_LEN       <- 4      # 4 weeks
WEEKS_PER_YEAR <- 52

## ----------------------------------------------------------
## List of mean estimators
## ----------------------------------------------------------

#Main function for St-MSh estimator
St_MSh_estimator <- function(x) {
  Tn <- nrow(x); N_T <- ncol(x)
  x_bar <- colMeans(x)
  a1T_hat <- sum(x_bar^2) / N_T
  trace_Sn <- sum(diag(cov(x)))
  a2T_hat <- a1T_hat + trace_Sn / (Tn * N_T)
  c_star_hat <- a1T_hat / a2T_hat
  as.vector(c_star_hat * x_bar)
}

#Main function for D-MSh estimator
D_MSh_estimator <- function(x) {
  Tn <- nrow(x)
  x_bar <- colMeans(x)
  Sn <- cov(x)
  c_star_hat <- (x_bar^2) / (x_bar^2 + (1 / Tn) * diag(Sn))
  as.vector(c_star_hat * x_bar)
}

#Main function for O-LSh estimator
O_LSh_estimator <- function(x) {
  Tn <- nrow(x); N_T <- ncol(x)
  x_bar <- colMeans(x)
  a1T_hat <- sum(x_bar^2) / N_T
  trace_Sn <- sum(diag(cov(x)))
  a2T_hat <- a1T_hat + trace_Sn / (Tn * N_T)
  s_i <- rowSums(x)
  term1 <- sum(s_i^2)
  term2 <- ((sum(s_i))^2 - term1) / (Tn - 1)
  d2_hat <- (term1 - term2) / (Tn^2 * N_T^2)
  grand_mean <- mean(x_bar)
  d3_hat <- sum((x_bar - grand_mean)^2) / N_T
  tilde_a1_hat <- a2T_hat - a1T_hat - d2_hat
  tilde_a2_hat <- tilde_a1_hat + d3_hat
  delta_star_hat <- tilde_a1_hat / tilde_a2_hat
  as.vector((1 - delta_star_hat) * x_bar + delta_star_hat * rep(grand_mean, N_T))
}

#Main function for T-LSh estimator
T_LSh_estimator <- function(x) {
  Tn <- nrow(x); N_T <- ncol(x)
  x_bar <- colMeans(x)
  a1T_hat <- sum(x_bar^2) / N_T
  trace_Sn <- sum(diag(cov(x)))
  a2T_hat <- a1T_hat + trace_Sn / (Tn * N_T)
  s_i <- rowSums(x)
  term1 <- sum(s_i^2)
  term2 <- ((sum(s_i))^2 - term1) / (Tn - 1)
  d2_hat <- (term1 - term2) / (Tn^2 * N_T^2)
  grand_mean <- mean(x_bar)
  d3_hat <- sum((x_bar - grand_mean)^2) / N_T
  d4_hat <- grand_mean^2
  tilde_a1_hat <- a2T_hat - a1T_hat - d2_hat
  tilde_a2_hat <- tilde_a1_hat + d3_hat
  delta_star_hat <- tilde_a1_hat / tilde_a2_hat
  xi_star_hat <- 1 - (d2_hat / (d2_hat + d4_hat)) / delta_star_hat
  as.vector((1 - delta_star_hat) * x_bar +
              delta_star_hat * xi_star_hat * rep(grand_mean, N_T))
}

# Estimator from (Wang et al., 2014)
Wang_estimator <- function(x) {
  n <- nrow(x); p <- ncol(x)
  x_bar <- colMeans(x)
  col_sum <- colSums(x)
  total_sq_col <- sum(col_sum^2)
  sum_sq_x <- sum(x^2)
  Y1n <- (total_sq_col - sum_sq_x) / (p * (n - 1))
  Y2n <- (sum_sq_x - p * Y1n) / (n * p)
  row_sum <- rowSums(x)
  total_row <- sum(row_sum)
  Y3n <- (total_row^2 - sum(row_sum^2)) / (p^2 * (n - 1))
  Y4n <- total_row / (n * p)
  denom <- Y1n + Y2n - Y3n
  alpha_star <- as.numeric((Y1n - Y3n) / denom)
  beta_star  <- as.numeric((Y2n * Y4n) / denom)
  as.vector(alpha_star * x_bar + beta_star * rep(1, p))
}

# Estimator from (Bodnar et al., 2019): BOP
BOP_estimator <- function(x) {
  n <- nrow(x); p <- ncol(x)
  y_bar <- colMeans(x)
  S_n <- cov(x)
  S_n_inv <- if (p > n) MASS::ginv(S_n) else solve(S_n)
  mu_0 <- rep(1, p)
  yS_inv_y      <- as.numeric(t(y_bar) %*% S_n_inv %*% y_bar)
  mu0_S_inv_mu0 <- as.numeric(t(mu_0) %*% S_n_inv %*% mu_0)
  yS_inv_mu0    <- as.numeric(t(y_bar) %*% S_n_inv %*% mu_0)
  term1 <- if (p > n) yS_inv_y - n / (p - n) else yS_inv_y - p / (n - p)
  term2 <- yS_inv_y * mu0_S_inv_mu0
  alpha_mean <- (term1 * mu0_S_inv_mu0 - yS_inv_mu0^2) / (term2 - yS_inv_mu0^2)
  beta_mean  <- (1 - alpha_mean) * yS_inv_mu0 / mu0_S_inv_mu0
  as.vector(alpha_mean * y_bar + beta_mean * mu_0)
}

# CW estimator (Bodnar et al., 2022) tested but EXCLUDED: at n > p the sample
# covariance is full rank, so S %*% ginv(S) = I and the (I - S S^+) y_bar term
# vanishes. CW then collapses to a scalar multiple of the sample mean and
# reproduces the Sample estimator's Sharpe ratios exactly. Not used.
# CW_estimator <- function(x) {
#   n <- nrow(x)
#   p <- ncol(x)
#   y_bar <- colMeans(x)
#   S_n <- cov(x)
#   S_n_plus <- MASS::ginv(S_n)
#   a <- (n - 3) / (p - n + 4)
#   scalar <- as.numeric(t(y_bar) %*% S_n_plus %*% y_bar)
#   term2 <- max(1 - (a / scalar), 0)
#   theta_hat_CW <- (diag(p) - S_n %*% S_n_plus) %*% y_bar +
#     term2 * (S_n %*% S_n_plus %*% y_bar)
#   return(as.vector(theta_hat_CW))
# }

# Estimator from (Jorion, 1986): Jorion
Jorion_estimator <- function(x) {
  Tn <- nrow(x); N_T <- ncol(x)
  mu_MLE <- colMeans(x)
  Sn <- (Tn - 1) / (Tn - N_T - 2) * cov(x)
  Sn_inv <- solve(Sn)
  ones <- rep(1, N_T)
  num <- as.numeric(t(ones) %*% Sn_inv %*% mu_MLE)
  denom <- as.numeric(t(ones) %*% Sn_inv %*% ones)
  mu0 <- num / denom
  diff <- mu_MLE - mu0 * ones
  d <- as.numeric(Tn * t(diff) %*% Sn_inv %*% diff)
  w <- (N_T + 2) / (N_T + 2 + d)
  as.vector((1 - w) * mu_MLE + w * rep(mu0, N_T))
}

## ----------------------------------------------------------
## Covariance shrinkage estimator: cov1Para (Ledoit-Wolf 2004b)
## ----------------------------------------------------------
cov1Para <- function(Y, k = -1) {
  dim.Y <- dim(Y); N <- dim.Y[1]; p <- dim.Y[2]
  if (k < 0) { Y <- scale(Y, scale = FALSE); k <- 1 }
  n <- N - k
  sample <- (t(Y) %*% Y) / n
  meanvar <- mean(diag(sample))
  target <- meanvar * diag(p)
  Y2 <- Y^2
  sample2 <- (t(Y2) %*% Y2) / n
  piMat <- sample2 - sample^2
  pihat <- sum(piMat)
  gammahat <- norm(c(sample - target), type = "2")^2
  rhohat <- 0
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  shrinkage * target + (1 - shrinkage) * sample
}

## ----------------------------------------------------------
## Annualised Sharpe ratio
## ----------------------------------------------------------
sharpe_annual <- function(r) {
  r <- as.numeric(r)
  (mean(r) / sd(r)) * sqrt(WEEKS_PER_YEAR)
}

############################################################
# A1) DATA PREPROCESSING
############################################################

asset_data <- read.csv("DA441_weekly_returns.csv", stringsAsFactors = FALSE)
asset_data$Date <- as.Date(as.character(asset_data$Date), format = "%Y%m%d")

rf_data <- read.csv("risk_free_rate_weekly.csv", stringsAsFactors = FALSE)
rf_data$Date <- as.Date(as.character(rf_data$Date), format = "%Y%m%d")
rf_data$RF <- rf_data$RF / 100

rf_aligned <- rf_data$RF[match(asset_data$Date, rf_data$Date)]

set.seed(GROUP_ID)
all_stock_cols <- 2:ncol(asset_data)
sel_cols <- sort(sample(all_stock_cols, 100))
selected_names <- colnames(asset_data)[sel_cols]

dates   <- asset_data$Date
returns <- as.matrix(asset_data[, sel_cols])
excess_returns <- returns - rf_aligned

cat("== A1 ==\n")
cat("Total weekly observations:", nrow(excess_returns), "\n")
cat("Number of selected stocks:", ncol(excess_returns), "\n")
cat("First 6 selected tickers:", paste(head(selected_names), collapse = ", "), "\n")
cat("Excess-return matrix dims:", paste(dim(excess_returns), collapse = " x "), "\n\n")

############################################################
# A2) ESTIMATOR SELECTION
############################################################

mean_estimators <- list(
  Sample = function(x) colMeans(x),
  St_MSh = St_MSh_estimator,
  D_MSh  = D_MSh_estimator,
  O_LSh  = O_LSh_estimator,
  T_LSh  = T_LSh_estimator,
  Wang   = Wang_estimator,
  BOP    = BOP_estimator,
  Jorion = Jorion_estimator
)
mean_order <- names(mean_estimators)

############################################################
# A3) ROLLING WINDOW + PORTFOLIO CONSTRUCTION
############################################################

n_total <- nrow(excess_returns)
max_windows <- floor((n_total - TRAIN_LEN - TEST_LEN) / STEP_LEN) + 1
cat("== A3 ==\n")
cat("Total number of rolling windows:", max_windows, "\n\n")

gam1 <- 1
sharpes_mv  <- vector("list", max_windows)
sharpes_tan <- vector("list", max_windows)

for (w in seq_len(max_windows)) {
  start_idx <- 1 + (w - 1) * STEP_LEN
  train_end <- start_idx + TRAIN_LEN - 1
  test_start <- train_end + 1
  test_end   <- train_end + TEST_LEN
  
  train_returns <- excess_returns[start_idx:train_end, , drop = FALSE]
  test_returns  <- excess_returns[test_start:test_end, , drop = FALSE]
  test_date     <- dates[test_start]
  
  Sigma <- cov1Para(train_returns)
  Sigma_inv <- solve(Sigma)
  N_assets <- ncol(train_returns)
  ones <- rep(1, N_assets)
  
  mu_list <- lapply(mean_estimators, function(f) f(train_returns))
  
  sr_mv <- sapply(mu_list, function(mu) {
    w_mv <- (1 / gam1) * (Sigma_inv %*% mu)
    sharpe_annual(test_returns %*% w_mv)
  })
  
  sr_tan <- sapply(mu_list, function(mu) {
    raw_w <- Sigma_inv %*% mu
    w_tan <- as.vector(raw_w / sum(raw_w))
    sharpe_annual(test_returns %*% w_tan)
  })
  
  w_ewp <- rep(1 / N_assets, N_assets)
  sr_ewp <- sharpe_annual(test_returns %*% w_ewp)
  
  sharpes_mv[[w]]  <- data.frame(Date = test_date,
                                 t(as.data.frame(c(sr_mv,  EWP = sr_ewp))),
                                 row.names = NULL)
  sharpes_tan[[w]] <- data.frame(Date = test_date,
                                 t(as.data.frame(c(sr_tan, EWP = sr_ewp))),
                                 row.names = NULL)
}

sharpes_mv  <- do.call(rbind, sharpes_mv)
sharpes_tan <- do.call(rbind, sharpes_tan)

############################################################
# A4) PERFORMANCE SUMMARY TABLES
############################################################
summary_table <- function(df, row_order) {
  df$Date <- as.Date(df$Date)
  idx_full     <- rep(TRUE, nrow(df))
  idx_crisis   <- df$Date >= as.Date("2008-10-06") & df$Date <= as.Date("2011-08-03")
  idx_pandemic <- df$Date >= as.Date("2020-02-26")
  get_stats <- function(idx) {
    sub <- df[idx, row_order, drop = FALSE]
    data.frame(Mean = sapply(sub, mean, na.rm = TRUE),
               SD   = sapply(sub, sd,   na.rm = TRUE))
  }
  out <- cbind(get_stats(idx_full), get_stats(idx_crisis), get_stats(idx_pandemic))
  colnames(out) <- c("Full_Mean","Full_SD","Crisis_Mean","Crisis_SD",
                     "Pandemic_Mean","Pandemic_SD")
  round(out, 4)
}

row_order <- c(mean_order, "EWP")
A4_mv  <- summary_table(sharpes_mv,  row_order)
A4_tan <- summary_table(sharpes_tan, row_order)

cat("== A4: Mean-Variance Portfolio (gamma = 1) - Annualised Sharpe ==\n")
print(A4_mv)
cat("\n== A4: Tangency Portfolio - Annualised Sharpe ==\n")
print(A4_tan)

write.csv(A4_mv,  "A4_MeanVariance_gamma1.csv")
write.csv(A4_tan, "A4_Tangency.csv")
saveRDS(list(mv = sharpes_mv, tan = sharpes_tan,
             selected = selected_names, n_windows = max_windows),
        "cw1_results.rds")

############################################################
# A5 SUPPORT: GICS sector composition of the 100 stocks
############################################################
gics_ok <- requireNamespace("readxl", quietly = TRUE)
if (gics_ok) {
  gics <- readxl::read_excel("DA441_GICS_Sectors.xlsx")
  cat("\n== A5 support: GICS columns ==\n"); print(colnames(gics))
}
cat("\nDONE.\n")



############################################################
############################################################
##                                                        ##
##                    COURSEWORK 2                         ##
##   Shrinkage GLMs for Insurance Claim Severity          ##
##   Tasks B1-B3  |  GROUP 7                               ##
##                                                        ##
##   NOTE: Part 2 is self-contained. It reads the French  ##
##   MTPL datasets (freMTPL2freq.csv, freMTPL2sev.csv).   ##
##   The B1 cross-validation fits 9 models x 500 runs and ##
##   is computationally heavy; see the runner at the end. ##
##   Required packages: glm2, glmnet, MASS, splines,      ##
##   savvyGLM.                                            ##
############################################################
############################################################

suppressMessages({
  library(glm2); library(glmnet); library(MASS)
  library(splines); library(savvyGLM)
})
options(scipen = 999)

## ==========================================================
## B0) DATA PREP (Group 7)
## ==========================================================
freq <- read.csv("freMTPL2freq.csv")
sev  <- read.csv("freMTPL2sev.csv")
sev_agg <- aggregate(ClaimAmount ~ IDpol, data = sev, FUN = sum)
data <- merge(freq, sev_agg, by = "IDpol", sort = FALSE)

# Brief: keep policies with EXACTLY 1 claim and a positive amount
data_sev <- data[data$ClaimNb == 1 & data$ClaimAmount > 0, ]
set.seed(7)                                   # GROUP 7
data_sev <- data_sev[sample(nrow(data_sev), 1000), ]
data_sev$Severity <- data_sev$ClaimAmount / data_sev$ClaimNb

# Bounded log transforms for skewed continuous covariates
data_sev$VehPower_Log   <- log(pmax(pmin(data_sev$VehPower, 9), 1))
data_sev$Density_Log    <- log(pmax(data_sev$Density, 1))
data_sev$BonusMalus_Log <- log(pmin(data_sev$BonusMalus, 150))

## ---- Metric helpers ----
calc_rmse_win <- function(a, p, cap) sqrt(mean((pmin(a, cap) - pmin(p, cap))^2))
calc_mae_win  <- function(a, p, cap) mean(abs(pmin(a, cap) - pmin(p, cap)))

WeightedGini <- function(actual, weights, predicted) {
  df <- data.frame(actual, weights, predicted)
  df <- df[order(df$predicted, decreasing = TRUE), ]
  df$random <- cumsum(df$weights / sum(df$weights))
  tot <- sum(df$actual * df$weights)
  df$cumPosFound <- cumsum(df$actual * df$weights)
  df$Lorentz <- df$cumPosFound / tot
  n <- nrow(df)
  sum(df$Lorentz[-n] * df$random[-1]) - sum(df$Lorentz[-1] * df$random[-n])
}
NormalizedWeightedGini <- function(a, w, p) WeightedGini(a, w, p) / WeightedGini(a, w, a)

n_folds      <- 5
control_list <- list(maxit = 250, epsilon = 1e-6, trace = FALSE)
family_type  <- Gamma(link = "log")

## ==========================================================
## B1) OUT-OF-SAMPLE EVALUATION: 5-fold CV x 100 repeats = 500 runs
##     9 models: GLM2, Ridge (RR), Elastic Net (EN),
##                St, DSh, SR, GSR, LW, QIS
## ==========================================================
run_single_iteration <- function(i) {
  repeat_idx <- floor((i - 1) / n_folds) + 1
  fold_idx   <- ((i - 1) %% n_folds) + 1
  # Same split seed per repeat -> all 9 models see identical folds (paired)
  set.seed(123 + repeat_idx)
  folds <- sample(rep(1:n_folds, length.out = nrow(data_sev)))
  test_idx   <- which(folds == fold_idx)
  train_data <- data_sev[-test_idx, ]
  test_data  <- data_sev[test_idx, ]

  # Natural cubic splines: knots learned on train, projected onto test (no leakage)
  driv <- ns(train_data$DrivAge, df = 3); veh <- ns(train_data$VehAge, df = 3)
  train_data$DrivAge_Sp1 <- driv[,1]; train_data$DrivAge_Sp2 <- driv[,2]; train_data$DrivAge_Sp3 <- driv[,3]
  train_data$VehAge_Sp1  <- veh[,1];  train_data$VehAge_Sp2  <- veh[,2];  train_data$VehAge_Sp3  <- veh[,3]
  tdb <- predict(driv, test_data$DrivAge); tvb <- predict(veh, test_data$VehAge)
  test_data$DrivAge_Sp1 <- tdb[,1]; test_data$DrivAge_Sp2 <- tdb[,2]; test_data$DrivAge_Sp3 <- tdb[,3]
  test_data$VehAge_Sp1  <- tvb[,1]; test_data$VehAge_Sp2  <- tvb[,2]; test_data$VehAge_Sp3  <- tvb[,3]

  # Credibility-weighted leave-one-out target encoding on train; dictionary applied to test
  te_cols <- c("Region", "Area", "VehBrand")
  gmean <- sum(train_data$ClaimAmount) / sum(train_data$ClaimNb); k <- 20
  for (col in te_cols) {
    ga <- tapply(train_data$ClaimAmount, train_data[[col]], sum)
    gc <- tapply(train_data$ClaimNb,     train_data[[col]], sum)
    oa <- train_data$ClaimAmount; oc <- train_data$ClaimNb
    gi <- as.character(train_data[[col]])
    la <- ga[gi] - oa; lc <- gc[gi] - oc          # leave-one-out
    lam <- lc / (lc + k); mu <- ifelse(lc > 0, la / lc, gmean)
    train_data[[paste0("TE_", col)]] <- lam * mu + (1 - lam) * gmean
    tgi <- as.character(test_data[[col]])
    fl <- gc / (gc + k); fm <- ga / gc; dict <- fl * fm + (1 - fl) * gmean
    test_data[[paste0("TE_", col)]] <- dict[tgi]
    test_data[[paste0("TE_", col)]][is.na(test_data[[paste0("TE_", col)]])] <- gmean
  }

  # Standardise engineered predictors using TRAIN moments only
  sc <- c("BonusMalus_Log","VehPower_Log","Density_Log","DrivAge_Sp1","DrivAge_Sp2","DrivAge_Sp3",
          "VehAge_Sp1","VehAge_Sp2","VehAge_Sp3","TE_Region","TE_Area","TE_VehBrand")
  tm <- colMeans(train_data[, sc]); ts <- apply(train_data[, sc], 2, sd); ts[ts == 0] <- 1
  for (col in sc) {
    train_data[[col]] <- (train_data[[col]] - tm[col]) / ts[col]
    test_data[[col]]  <- (test_data[[col]]  - tm[col]) / ts[col]
  }

  glm_formula <- Severity ~ BonusMalus_Log + VehPower_Log + Density_Log +
    DrivAge_Sp1 + DrivAge_Sp2 + DrivAge_Sp3 +
    VehAge_Sp1 + VehAge_Sp2 + VehAge_Sp3 + VehGas +
    TE_Region + TE_Area + TE_VehBrand
  X_train <- model.matrix(glm_formula, train_data)
  X_test  <- model.matrix(glm_formula, test_data)
  y_train <- train_data$Severity; y_test <- test_data$Severity
  w_train <- train_data$ClaimNb;  w_test <- test_data$ClaimNb
  cap_99 <- quantile(y_train, .99, na.rm = TRUE)
  cap_95 <- quantile(y_train, .95, na.rm = TRUE)

  eval_model <- function(coefs) {
    if (all(is.na(coefs))) return(rep(NA_real_, 10))
    s <- ifelse(is.na(coefs), 0, coefs)
    p <- exp(as.numeric(X_test %*% s))          # inverse log-link
    c(sqrt(mean((y_test - p)^2)), mean(abs(y_test - p)),
      calc_rmse_win(y_test, p, cap_99), calc_rmse_win(y_test, p, cap_95),
      calc_mae_win(y_test, p, cap_99),  calc_mae_win(y_test, p, cap_95),
      WeightedGini(y_test, w_test, p), NormalizedWeightedGini(y_test, w_test, p),
      sum(w_test * Gamma()$dev.resids(y_test, p, 1)),
      sum(y_test * w_test) / sum(p * w_test))
  }

  t1 <- system.time({ m1 <- glm.fit2(X_train, y_train, weights = w_train, family = family_type, control = control_list) }); res1 <- eval_model(m1$coefficients); it1 <- m1$iter
  t2 <- system.time({ m2 <- cv.glmnet(X_train[,-1], y_train, weights = w_train, family = family_type, alpha = 0,   nlambda = 50, nfolds = 5) }); res2 <- eval_model(as.vector(coef(m2, s = "lambda.min"))); it2 <- if (is.null(m2$glmnet.fit$npasses)) NA else m2$glmnet.fit$npasses[1]
  t3 <- system.time({ m3 <- cv.glmnet(X_train[,-1], y_train, weights = w_train, family = family_type, alpha = 0.5, nlambda = 50, nfolds = 5) }); res3 <- eval_model(as.vector(coef(m3, s = "lambda.min"))); it3 <- if (is.null(m3$glmnet.fit$npasses)) NA else m3$glmnet.fit$npasses[1]
  fit_sv <- function(mc) savvy_glm.fit2(x = X_train, y = y_train, weights = w_train, family = family_type, model_class = mc, control = control_list)
  t4 <- system.time({ m4 <- fit_sv("St")  }); res4 <- eval_model(m4$coefficients); it4 <- m4$iter
  t5 <- system.time({ m5 <- fit_sv("DSh") }); res5 <- eval_model(m5$coefficients); it5 <- m5$iter
  t6 <- system.time({ m6 <- fit_sv("SR")  }); res6 <- eval_model(m6$coefficients); it6 <- m6$iter
  t7 <- system.time({ m7 <- fit_sv("GSR") }); res7 <- eval_model(m7$coefficients); it7 <- m7$iter
  t8 <- system.time({ m8 <- fit_sv("LW")  }); res8 <- eval_model(m8$coefficients); it8 <- m8$iter
  t9 <- system.time({ m9 <- fit_sv("QIS") }); res9 <- eval_model(m9$coefficients); it9 <- m9$iter

  models  <- c("GLM2","RR","EN","QIS","LW","SR","GSR","St","DSh")
  all_res <- rbind(res1, res2, res3, res9, res8, res6, res7, res4, res5)
  times   <- c(t1[3], t2[3], t3[3], t9[3], t8[3], t6[3], t7[3], t4[3], t5[3])
  iters   <- c(it1, it2, it3, it9, it8, it6, it7, it4, it5)
  data.frame(Iteration = i, Model = models,
    RMSE = all_res[,1], MAE = all_res[,2], RMSE_Win_99 = all_res[,3], RMSE_Win_95 = all_res[,4],
    MAE_Win_99 = all_res[,5], MAE_Win_95 = all_res[,6], Gini_Raw = all_res[,7], Gini_Norm = all_res[,8],
    Dev = all_res[,9], AE_Ratio = all_res[,10], CT = as.numeric(times), NoIt = as.numeric(iters))
}

## ---- Run all 500 evaluations (heavy). For chunked/resumable runs, wrap in a
##      loop over blocks of i and saveRDS each block instead. ----
run_B1 <- function() {
  final_results <- do.call(rbind, lapply(1:500, run_single_iteration))
  saveRDS(final_results, "final_results.rds")

  ord <- c("GLM2","RR","EN","St","DSh","SR","GSR","LW","QIS")   # coursework order
  final_results$Model <- factor(final_results$Model, levels = ord)

  # Raw absolute table (mean over 500 runs)
  avg_raw <- aggregate(
    cbind(RMSE, MAE, RMSE_Win_99, RMSE_Win_95, MAE_Win_99, MAE_Win_95,
          Dev, Gini_Raw, Gini_Norm, AE_Ratio, CT, NoIt) ~ Model,
    data = final_results, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
  write.csv(avg_raw, "tab_raw.csv", row.names = FALSE)

  # Relative table: per-run ratio to GLM2, then averaged
  ratio_list <- lapply(split(final_results, final_results$Iteration), function(df) {
    g <- function(m) df[[m]][df$Model == "GLM2"]
    base <- list(RMSE = g("RMSE"), MAE = g("MAE"), R99 = g("RMSE_Win_99"), R95 = g("RMSE_Win_95"),
                 M99 = g("MAE_Win_99"), M95 = g("MAE_Win_95"), Dev = g("Dev"),
                 CT = ifelse(is.na(g("CT")) | g("CT") == 0, NA, g("CT")),
                 NoIt = ifelse(is.na(g("NoIt")) | g("NoIt") == 0, NA, g("NoIt")))
    df$RMSE_R <- df$RMSE / base$RMSE; df$MAE_R <- df$MAE / base$MAE
    df$R99_R <- df$RMSE_Win_99 / base$R99; df$R95_R <- df$RMSE_Win_95 / base$R95
    df$M99_R <- df$MAE_Win_99 / base$M99; df$M95_R <- df$MAE_Win_95 / base$M95
    df$Dev_R <- df$Dev / base$Dev; df$CT_R <- df$CT / base$CT; df$NoIt_R <- df$NoIt / base$NoIt
    df
  })
  ar <- do.call(rbind, ratio_list); ar$Model <- factor(ar$Model, levels = ord)
  avg_rel <- aggregate(
    cbind(RMSE_R, MAE_R, R99_R, R95_R, M99_R, M95_R, Dev_R,
          Gini_Raw, Gini_Norm, AE_Ratio, CT_R, NoIt_R) ~ Model,
    data = ar, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
  write.csv(avg_rel, "tab_rel.csv", row.names = FALSE)

  cat("B1 done: tab_raw.csv and tab_rel.csv written.\n")
  list(raw = avg_raw, rel = avg_rel)
}

## ==========================================================
## B3) IN-SAMPLE ANALYSIS of the selected model (Ridge)
##     Refit on all 1,000 obs; compare fitted vs actual.
## NOTE: in-sample target encoding here is the plain group mean
##       (not leave-one-out); this matches the B3 in-sample brief
##       but is slightly more optimistic than the CV pipeline.
## ==========================================================
run_B3 <- function() {
  d <- data_sev
  driv <- ns(d$DrivAge, df = 3); veh <- ns(d$VehAge, df = 3)
  d$DrivAge_Sp1 <- driv[,1]; d$DrivAge_Sp2 <- driv[,2]; d$DrivAge_Sp3 <- driv[,3]
  d$VehAge_Sp1  <- veh[,1];  d$VehAge_Sp2  <- veh[,2];  d$VehAge_Sp3  <- veh[,3]
  te <- c("Region","Area","VehBrand"); gm <- sum(d$ClaimAmount) / sum(d$ClaimNb); k <- 20
  for (col in te) {
    ga <- tapply(d$ClaimAmount, d[[col]], sum); gc <- tapply(d$ClaimNb, d[[col]], sum)
    gi <- as.character(d[[col]]); la <- ga[gi] - d$ClaimAmount; lc <- gc[gi] - d$ClaimNb
    lam <- lc / (lc + k); mu <- ifelse(lc > 0, la / lc, gm)
    d[[paste0("TE_", col)]] <- lam * mu + (1 - lam) * gm
  }
  sc <- c("BonusMalus_Log","VehPower_Log","Density_Log","DrivAge_Sp1","DrivAge_Sp2","DrivAge_Sp3",
          "VehAge_Sp1","VehAge_Sp2","VehAge_Sp3","TE_Region","TE_Area","TE_VehBrand")
  mn <- colMeans(d[, sc]); s <- apply(d[, sc], 2, sd); s[s == 0] <- 1
  for (col in sc) d[[col]] <- (d[[col]] - mn[col]) / s[col]

  f <- Severity ~ BonusMalus_Log + VehPower_Log + Density_Log +
    DrivAge_Sp1 + DrivAge_Sp2 + DrivAge_Sp3 +
    VehAge_Sp1 + VehAge_Sp2 + VehAge_Sp3 + VehGas + TE_Region + TE_Area + TE_VehBrand
  X <- model.matrix(f, d); y <- d$Severity; w <- d$ClaimNb
  set.seed(7)
  cvfit <- cv.glmnet(X[,-1], y, weights = w, family = Gamma(link = "log"),
                     alpha = 0, nlambda = 50, nfolds = 5)
  co <- as.vector(coef(cvfit, s = "lambda.min"))
  pred <- exp(as.numeric(X %*% ifelse(is.na(co), 0, co)))
  rmse <- sqrt(mean((y - pred)^2)); mae <- mean(abs(y - pred))
  ae <- sum(y * w) / sum(pred * w); dev <- sum(w * Gamma()$dev.resids(y, pred, 1))
  cat(sprintf("IS RMSE=%.2f MAE=%.2f AE=%.4f Dev=%.2f lambda.min=%.4f\n",
              rmse, mae, ae, dev, cvfit$lambda.min))
  cat(sprintf("Actual mean=%.2f Pred mean=%.2f | Actual median=%.2f Pred median=%.2f\n",
              mean(y), mean(pred), median(y), median(pred)))

  dr <- sign(y - pred) * sqrt(w * Gamma()$dev.resids(y, pred, 1))
  pdf("fig_insample.pdf", width = 10, height = 4.2)
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  plot(pred, y, log = "xy", pch = 16, col = rgb(0.15, 0.35, 0.75, 0.45),
       main = "In-Sample: Actual vs Predicted", xlab = "Predicted severity (EUR)", ylab = "Actual severity (EUR)")
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  plot(pred, dr, pch = 16, col = rgb(0.80, 0.20, 0.20, 0.45),
       main = "Deviance Residuals vs Predicted", xlab = "Predicted severity (EUR)", ylab = "Deviance residual")
  abline(h = 0, col = "black", lwd = 2, lty = 2); lines(lowess(pred, dr), col = "blue", lwd = 2)
  dev.off()
  cat("figure saved: fig_insample.pdf\n")
}

## ---- To reproduce Part 2 end to end, uncomment: ----
# b1 <- run_B1()   # heavy: 9 models x 500 runs
# run_B3()         # in-sample fit + figure
