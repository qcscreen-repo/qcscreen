library(fbi)
library(Rfast)
library(doParallel)
library(foreach)
library(quantdr)
library(QCSIS)

filepath <- "2025-08-MD.csv"
data <- fredmd(
  filepath,
  date_start = as.Date("2016-01-01"),
  date_end   = as.Date("2024-02-01"),
  transform  = TRUE
)
N <- ncol(data)


data_clean <- rm_outliers.fredmd(data)
col_na_prop <- apply(is.na(data_clean), 2, mean)
data_select <- data_clean[, (col_na_prop < 0.05)]
data_bal <- na.omit(data_select)
X_bal <- data_bal[,2:ncol(data_bal)]

############## main functions #################
# (2.6) Conditional copula function
conditional_copula_C2_given_C1 <- function(u, v, s_vec, t_vec, h) {
  n <- length(s_vec)
  z_u <- qnorm(u)
  z_v <- qnorm(v)
  kernel_s <- dnorm((z_u - s_vec) / h)
  kernel_t <- pnorm((z_v - t_vec) / h)
  numerator <- sum(kernel_s * kernel_t)
  denominator <- n * h * dnorm(z_u)
  numerator / denominator
}

# (2.8) Copula quantile dependence (vectorized over u_hat)
estimate_D_star <- function(s_vec, t_vec, x_vec, h = 0.1, alpha = 0.5, omega = function(u) 1) {
  n <- length(s_vec)
  F_Xj <- ecdf(x_vec)
  u_hat <- (n / (n + 1)) * F_Xj(x_vec)
  
  # Evaluate C(u_i, alpha | s_vec, t_vec) for all u_hat at once
  C_hat_vec <- vapply(
    u_hat,
    FUN = function(u_i) conditional_copula_C2_given_C1(u = u_i, v = alpha, s_vec = s_vec, t_vec = t_vec, h = h),
    FUN.VALUE = numeric(1)
  )
  mean((C_hat_vec - alpha)^2 * omega(u_hat))
}

# phi'(Phi^{-1}(alpha)) and related constants / variances
phi_prime_at_alpha <- function(alpha) {
  z <- qnorm(alpha)
  -z * dnorm(z)                 # φ'(z) = -z φ(z)
}

# Returns M_{perp alpha,2}, M_{perp alpha,4}, sigma^2_{perp alpha,1}, sigma^2_{perp alpha,3}
alpha_constants <- function(alpha) {
  stopifnot(alpha > 0, alpha < 1)
  phip <- phi_prime_at_alpha(alpha)
  
  # M_{perp alpha,2}(omega) = (alpha - alpha^2) / (4*pi)
  M2 <- (alpha - alpha^2) / (4 * pi)
  
  # sigma^2_{perp alpha,3}(omega) = [alpha - alpha^2]^2 / (16 pi^2)
  sig3_sq <- (alpha - alpha^2)^2 / (16 * pi^2)
  
  # M_{perp alpha,4}(omega) = (2*alpha^2 - 4*alpha*phi' + 3*(phi')^2) / (24*sqrt(3)*pi)
  M4 <- (2*alpha^2 - 4*alpha*phip + 3*(phip^2)) / (24 * sqrt(3) * pi)
  
  # sigma^2_{perp alpha,1}(omega)
  # = (alpha - alpha^2) * [ (18*alpha^2 - 40*alpha*phip + 25*phip^2)/(400*pi^2*sqrt(5)) - (3*phip - 2*alpha)^2/(432*pi^2) ]
  term1 <- (18*alpha^2 - 40*alpha*phip + 25*(phip^2)) / (400 * pi^2 * sqrt(5))
  term2 <- (3*phip - 2*alpha)^2 / (432 * pi^2)
  sig1_sq <- (alpha - alpha^2) * (term1 - term2)
  
  list(M2 = M2, M4 = M4, sig1_sq = sig1_sq, sig3_sq = sig3_sq)
}

# Calls estimate_D_star() for each predictor column
compute_D_dagger_from_alpha <- function(X, Y, h, alpha, omega=function(u) dnorm(qnorm(u))^2,
                                        delta_n_alpha = 0, eps = 1e-12, n_workers = max(1, parallel::detectCores() - 5)) {
  stopifnot(is.matrix(X), length(Y) == nrow(X))
  n <- length(Y); p <- ncol(X)
  
  # 1) Compute hat D_alpha for each j
  F_Y <- ecdf(Y)
  t_vec <- qnorm((n / (n + 1)) * F_Y(Y))  # Transformed Y
  
  # Parallel cluster
  cl <- parallel::makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  
  D_alpha_vals <- foreach(
    j = 1:p, .combine = c, .packages = "stats",
    .export = c("estimate_D_star", "conditional_copula_C2_given_C1")
  ) %dopar% {
    F_Xj <- ecdf(X[, j])
    s_vec <- qnorm((n / (n + 1)) * F_Xj(X[, j]))
    estimate_D_star(s_vec, t_vec, X[, j], h = h, alpha = alpha, omega = omega)
  }
  
  stopCluster(cl)
  
  # 2) Asymptotic constants (depend on alpha only)
  consts <- alpha_constants(alpha)
  A_num  <- n^(-1) * h^(-1) * consts$M2 + h^4 * consts$M4
  denom  <- sqrt(4 * n^(-1) * h^4 * consts$sig1_sq +
                   2 * n^(-2) * h^(-1) * consts$sig3_sq)
  
  if (any(denom <= 0)) stop("Denominator non-positive; check α/h.")
  if (any(A_num <= 0))  warning("A_num ≤ 0; log-ratio may be unstable.")
  
  log_ratio <- log( pmax(D_alpha_vals, eps) / pmax(A_num, eps) )
  D_dagger  <- (A_num / denom) * log_ratio + as.numeric(delta_n_alpha)
  
  list(
    D_alpha  = D_alpha_vals,
    M2 = consts$M2, M4 = consts$M4,
    sig1_sq = consts$sig1_sq, sig3_sq = consts$sig3_sq,
    A_num = A_num, denom = denom,
    D_dagger = D_dagger
  )
}

# Inputs:
#   D_dagger: length-p vector of hat D^dagger_alpha
#   gammas:   vector of gamma values for the adaptive rule
# Output: selected indices, thresholds tau, and per-gamma details
adaptive_threshold_select <- function(D_dagger, gammas) {
  gammas <- as.numeric(gammas)
  p <- length(D_dagger)
  
  # Precompute sorted |D_dagger| grid (independent of gamma)
  absD     <- abs(D_dagger)
  t_sorted <- sort(absD, decreasing = FALSE)
  t_all    <- c(0, t_sorted, Inf)   # t_0=0, t_{p+1}=Inf
  
  count_leq_neg_t <- function(t) sum(D_dagger <= -t)
  count_geq_pos_t <- function(t) sum(D_dagger >=  t)
  
  # For fixed gamma: find tau, m, and selected predictors
  solve_for_gamma <- function(g) {
    m <- 0L
    while (m <= p) {
      t_m <- t_all[m + 1L]
      num <- 1L + count_leq_neg_t(t_m)
      den <- count_geq_pos_t(t_m)
      ratio <- if (den == 0L) Inf else num / den
      if (ratio > g && m <= p) {
        m <- m + 1L
      } else {
        break
      }
    }
    tau <- t_all[m]
    selected <- which(D_dagger >= tau)
    list(gamma = g, tau = tau, m = m, selected_idx = selected)
  }
  
  # Run the rule for each gamma
  details_list <- lapply(gammas, solve_for_gamma)
  names(details_list) <- paste0("gamma=", format(gammas, trim = TRUE))
  
  # Summary table across gammas
  summary_df <- data.frame(
    gamma = vapply(details_list, function(z) z$gamma, numeric(1)),
    tau   = vapply(details_list, function(z) z$tau,   numeric(1)),
    m     = vapply(details_list, function(z) z$m,     integer(1)),
    selected_size = vapply(details_list, function(z) length(z$selected_idx), integer(1))
  )
  
  list(
    summary = summary_df,
    details = details_list,
    D_dagger = D_dagger,
    t_grid = t_all
  )
}






############## compared methods ################
wdGscreening <- function(Y, X, ppp = 1,
                         n_workers = max(1, parallel::detectCores() - 5)) {
  Y <- as.matrix(Y)
  n <- nrow(X)
  
  Gaussian_temp <- qnorm((1:n) / (n + 1))
  gt <- function(sq) Gaussian_temp[sq]
  
  ry <- Rfast::Rank(Y)
  rx <- colRanks(as.matrix(X))  # parallel=TRUE omitted
  gry <- gt(ry)
  grx <- apply(rx, 2, gt)
  
  X <- as.matrix(grx)
  Y <- gry
  
  G <- Rfast::rmvnorm(n, mu = rep(0, 2), sigma = diag(2),seed=1234)
  
  # Parallel cluster
  cl <- parallel::makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  
  # Parallel Wasserstein distances
  wxy <- foreach::foreach(
    i = seq_len(ncol(X)),
    .combine = "c",
    .packages = c("transport")
  ) %dopar% {
    transport::wasserstein(transport::pp(G),
                           transport::pp(cbind(X[, i], Y)),
                           p = ppp)
  }
  
  stopCluster(cl)
  
  return(wxy)
}




################## Analysis time 1 ################
X_bal1=X_bal[6:47,]

Y=X_bal1$CPIAUCSL
X=as.matrix(X_bal1[, !(names(X_bal1) %in% c("CPIAUCSL", "CPIAPPSL", "CPITRNSL", "CPIMEDSL", "CUSR0000SAC",
                                          "CUSR0000SAD", "CUSR0000SAS", "CPIULFSL", "CUSR0000SA0L2", "CUSR0000SA0L5"))])

n=nrow(X)
p=ncol(X)
alpha <- 0.5
gammas_vec <- c(0.15)

hh <- alpha_constants(alpha)
h  <- (hh$M2/(4*hh$M4))^(1/5)*n^(-1/5)
if (alpha == 0.5){
  delta_n_alpha <- 2.858*n^(-0.1)
}else if (alpha == 0.25){
  delta_n_alpha <- 2.478*n^(-0.1)
}else {
  delta_n_alpha <- 2.455*n^(-0.1)
}

X_norm <- scale(X)
Y_norm <- scale(Y)

# --- Proposed method
comp <- compute_D_dagger_from_alpha(X_norm, Y_norm, h = h, alpha = alpha, delta_n_alpha = delta_n_alpha)
res  <- adaptive_threshold_select(comp$D_dagger, gammas = gammas_vec)
# Example selected_idx (one run): 26 43 44 61 64 95 96 97 98 99 102 104
res$details$`gamma=0.15`$selected_idx
# FRED series codes for selected predictors, ordered as in the vector above
colnames(X[,c(26,64,43,99,95,98,96,97,104,61,102,44)])

# --- SIS method ---
order(abs(cor(X, Y, method = "pearson")), decreasing = TRUE)[1:12]
# Corresponding variable names for that ranking (top 12)
colnames(X[,c(104,102,98,97,96,99,95,100,64,13,11,16)])

# --- WD Method ---
order(wdGscreening(Y, X, ppp = 1), decreasing = TRUE)[1:12]
colnames(X[,c(104,102,98,97,96,99,95,11,93,64,100,112)])

# --- QcSIS method ---
order(QCSIS(X, Y, alpha, n-1)$w, decreasing = TRUE)[1:12]
colnames(X[,c(102,104,96,98,97,99,64,95,36,100,105,94)])


################## Analysis time 2 ################
X_bal2=X_bal[48:89,]

Y=X_bal2$CPIAUCSL
X=as.matrix(X_bal2[, !(names(X_bal2) %in% c("CPIAUCSL", "CPIAPPSL", "CPITRNSL", "CPIMEDSL", "CUSR0000SAC",
                                            "CUSR0000SAD", "CUSR0000SAS", "CPIULFSL", "CUSR0000SA0L2", "CUSR0000SA0L5"))])


n=nrow(X)
p=ncol(X)
alpha <- 0.5
gammas_vec <- c(0.15)

hh <- alpha_constants(alpha)
h  <- (hh$M2/(4*hh$M4))^(1/5)*n^(-1/5)
if (alpha == 0.5){
  delta_n_alpha <- 2.858*n^(-0.1)
}else if (alpha == 0.25){
  delta_n_alpha <- 2.478*n^(-0.1)
}else {
  delta_n_alpha <- 2.455*n^(-0.1)
}

X_norm <- scale(X)
Y_norm <- scale(Y)

# --- Proposed method ---
comp <- compute_D_dagger_from_alpha(X_norm, Y_norm, h = h, alpha = alpha, delta_n_alpha = delta_n_alpha)
res  <- adaptive_threshold_select(comp$D_dagger, gammas = gammas_vec)
# Example selected_idx (one run): 3 10 21 22 36 44 47 49 54 56 57 58 61 74 80 95 96 97 98 100 102 104 105
res$details$`gamma=0.15`$selected_idx
colnames(X[,c(10,95,80,21,54,56,36,105,3,58,22,57,44,49,96,97,47,98,104,61,100,74,102)])

# --- SIS method ---
order(abs(cor(X, Y, method = "pearson")), decreasing = TRUE)[1:12]
colnames(X[,c(102,104,97,96,98,100,56,105,61,68,103,15)])

# --- WD Method ---
order(wdGscreening(Y, X, ppp = 1), decreasing = TRUE)[1:12]
colnames(X[,c(102,104,97,96,100,56,98,64,18,61,105,67)])

# --- QcSIS method ---
order(QCSIS(X, Y, alpha, n-1)$w, decreasing = TRUE)[1:12]
colnames(X[,c(102,104,97,96,100,98,56,105,58,69,15,61)])
