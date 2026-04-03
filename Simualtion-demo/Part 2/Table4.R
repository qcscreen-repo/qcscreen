pkgs <- c(
  "doParallel",
  "foreach",
  "Rfast",
  "dcov",
  "splines",
  "quantreg",
  "QCSIS",
  "transport"
)

lapply(pkgs, library, character.only = TRUE)

############## main functions (our method) #################
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

# (2.8) Copula Quantile Dependence 
estimate_D_star <- function(s_vec, t_vec, x_vec, h = 0.1, alpha = 0.5, omega = function(u) 1) {
  n <- length(s_vec)
  F_Xj <- ecdf(x_vec)
  u_hat <- (n / (n + 1)) * F_Xj(x_vec)
  
  # Compute C(u_i, alpha | s_vec, t_vec) for all u_hat at once
  C_hat_vec <- vapply(
    u_hat,
    FUN = function(u_i) conditional_copula_C2_given_C1(u = u_i, v = alpha, s_vec = s_vec, t_vec = t_vec, h = h),
    FUN.VALUE = numeric(1)
  )
  mean((C_hat_vec - alpha)^2 * omega(u_hat))
}

# Main QC-Screen function
qc_screen <- function(X, Y, alpha = 0.5, h = 0.1, omega = function(u) 1, t_n_alpha = 0.01) {
  n <- length(Y)
  p <- ncol(X)
  
  F_Y <- ecdf(Y)
  t_vec <- qnorm((n / (n + 1)) * F_Y(Y))  # Transformed Y
  
  D_alpha_vals <- numeric(p)
  
  for (j in 1:p) {
    F_Xj <- ecdf(X[, j])
    s_vec <- qnorm((n / (n + 1)) * F_Xj(X[, j]))
    
    # Pass original X_j to compute u_hat correctly
    D_alpha_vals[j] <- estimate_D_star(s_vec, t_vec, X[, j], h = h, alpha = alpha, omega = omega)
  }
  
  selected_vars <- which(D_alpha_vals > t_n_alpha)
  
  return(list(
    selected_index = selected_vars,
    D_alpha = D_alpha_vals
  ))
}

# Main function: outer parallel loop
qc_screen_with_min_model_size <- function(X, Y, alpha = 0.5, h = NULL,
                                          omega = function(u) dnorm(qnorm(u))^2,
                                          active_index = NULL, n_workers = max(1, parallel::detectCores() - 5)) {
  n <- length(Y)
  p <- ncol(X) 
  if (is.null(h)) {
    if (alpha==0.5){
      h <- 1.02 * n^(-1/5)
    }
    else{
      h <- 0.78* n^(-1/5)
    }
  }
  
  F_Y <- ecdf(Y)
  t_vec <- qnorm((n / (n + 1)) * F_Y(Y))  # Transformed Y
  
  # Parallel setting
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
  
  ranked_index <- order(D_alpha_vals, decreasing = TRUE)
  
  if (!is.null(active_index)) {
    rank_positions <- match(active_index, ranked_index)
    min_model_size <- max(rank_positions)
  } else {
    min_model_size <- NA
  }
  
  list(
    D_alpha = D_alpha_vals,
    ranked_index = ranked_index,
    min_model_size = min_model_size
  )
}



############## conditional functions #################
# -----------------------------
# Kernels
# -----------------------------
gaussian_kernel <- function(t) dnorm(t)

Kh <- function(u_diff, h, kernel = gaussian_kernel) {
  kernel(u_diff / h) / h
}

# -----------------------------
# Local linear quantile: Q_tau(Y|U=u0) ≈ a_hat
# -----------------------------
local_linear_qtau <- function(Y, U, tau, u0, h, kernel = gaussian_kernel) {
  w <- Kh(U - u0, h, kernel)
  if (sum(w) <= 1e-12) {
    return(as.numeric(quantile(Y, probs = tau, type = 7, na.rm = TRUE)))
  }
  Xloc <- cbind(1, U - u0)
  fit <- rq.wfit(x = Xloc, y = Y, tau = tau, weights = w)
  as.numeric(fit$coefficients[1])
}

local_linear_qtau_all <- function(Y, U, tau, h, kernel = gaussian_kernel) {
  n <- length(Y)
  qhat <- numeric(n)
  for (i in seq_len(n)) {
    qhat[i] <- local_linear_qtau(Y, U, tau, u0 = U[i], h = h, kernel = kernel)
  }
  qhat
}

# -----------------------------
# Nadaraya–Watson conditional mean at u0
# -----------------------------
nw_cond_mean <- function(g, U, u0, h, kernel = gaussian_kernel) {
  w <- Kh(U - u0, h, kernel)
  sw <- sum(w)
  if (sw <= 1e-12) return(NA_real_)
  sum(w * g) / sw
}

# Vectorized across all u=U_i for many columns at once
# g_mat: n x p matrix. returns n x p matrix of E[g|U=U_i]
nw_cond_mean_all_mat <- function(g_mat, U, h, kernel = gaussian_kernel) {
  n <- nrow(g_mat)
  p <- ncol(g_mat)
  out <- matrix(NA_real_, n, p)
  
  for (i in seq_len(n)) {
    w <- Kh(U - U[i], h, kernel)
    sw <- sum(w)
    if (sw <= 1e-12) next
    # Weighted column means: (w' g_mat)/sum(w)
    out[i, ] <- as.numeric(crossprod(w, g_mat) / sw)
  }
  out
}

# -----------------------------
# Main: X is n x p
# -----------------------------
estimate_qcor_tau_p <- function(Y, X, U, tau = 0.5,
                                h_q = NULL, h_k = NULL,
                                kernel = gaussian_kernel) {
  
  stopifnot(is.numeric(Y), is.numeric(U))
  stopifnot(is.matrix(X) || is.data.frame(X))
  X <- as.matrix(X)
  
  n <- length(Y); p <- ncol(X)
  stopifnot(nrow(X) == n, length(U) == n)
  
  # bandwidth defaults
  if (is.null(h_q)) h_q <- 1.06 * sd(U) * n^(-1/5)
  if (is.null(h_k)) h_k <- h_q
  
  # 1) qhat_i = Q_tau(Y|U=U_i)
  qhat <- local_linear_qtau_all(Y, U, tau = tau, h = h_q, kernel = kernel)
  
  # 2) psi_i = tau - 1(Y_i < qhat_i)
  psi <- tau - as.numeric(Y < qhat)
  
  # 3) Kernel regression estimates:
  #    E[psi * X | U], E[X|U], E[X^2|U]
  psiX <- X * psi                  # broadcast psi over columns
  EXpsiX_u <- nw_cond_mean_all_mat(psiX, U, h = h_k, kernel = kernel)
  EX_u     <- nw_cond_mean_all_mat(X,    U, h = h_k, kernel = kernel)
  EX2_u    <- nw_cond_mean_all_mat(X^2,  U, h = h_k, kernel = kernel)
  
  varX_u <- EX2_u - EX_u^2
  
  # 4) qcor: elementwise
  denom <- sqrt((tau - tau^2) * varX_u)
  qcor_u <- EXpsiX_u / denom
  qcor_u[!is.finite(qcor_u)] <- NA_real_
  
  # 4) omega*_k = (1/n) sum_j qcor^2 at U_j
  omega_star <- colMeans(qcor_u^2, na.rm = TRUE)
  
  list(
    omega_star = omega_star,  # length p
    qcor_u = qcor_u,          # n x p (optional, can be huge)
    qhat = qhat,
    h_q = h_q, h_k = h_k
  )
}

############## joint functions #################
# -----------------------------
# shrinkage: soft-thresholding with threshold thld
# MATLAB:
# if |a| > thld: z = a - sign(a)*thld else 0
# -----------------------------
shrinkage <- function(a, thld) {
  if (abs(a) > thld) {
    return(a - sign(a) * thld)
  } else {
    return(0)
  }
}

# -----------------------------
# ITA: "FISTA-like" solver used by the authors for initialization
# Inputs:
#   X: N x p
#   Y: length N
#   k: sparsity index used to pick threshold s(k) where s=sort(abs(x))
#   maxi: iterations
# Output:
#   x: length p
# -----------------------------
ITA <- function(X, Y, k, maxi) {
  n <- ncol(X)
  x <- rep(0, n)
  x1 <- x
  
  t_k <- 1
  t_km1 <- 1
  
  for (iter in seq_len(maxi)) {
    xs <- x
    # momentum step
    x <- x + ((t_km1 - 1) / t_k) * (x - x1)
    x1 <- xs
    
    # gradient step for least squares: x = x + X' (Y - X x)
    resid <- Y - as.numeric(X %*% x)
    x <- x + as.numeric(crossprod(X, resid))
    
    s <- sort(abs(x), decreasing = TRUE)
    # IMPORTANT: MATLAB uses s(k). If k > n, that would error there too.
    if (k > length(s)) stop("ITA: k is larger than number of features.")
    thld <- s[k]
    
    # coordinate-wise shrinkage
    for (j in seq_len(n)) {
      x[j] <- shrinkage(x[j], thld)
    }
    
    t_kp1 <- 0.5 * (1 + sqrt(1 + 4 * t_k * t_k))
    t_km1 <- t_k
    t_k <- t_kp1
  }
  return(x)
}

# -----------------------------
# Lassoinit: initialization for each component
# MATLAB:
#   normtrain = norm(X,'fro')*1.1;
#   tmp = ITA(X/normtrain, Y/normtrain, kt*S(kt)+1, maxitLasso)
# -----------------------------
Lassoinit <- function(X, Y, K, S, maxitLasso) {
  p <- ncol(X)
  normtrain <- norm(X, type = "F") * 1.1
  beta <- matrix(0, nrow = p, ncol = K)
  
  for (kt in seq_len(K)) {
    kk <- kt * S[kt] + 1
    beta[, kt] <- ITA(X / normtrain, Y / normtrain, kk, maxitLasso)
  }
  return(beta)
}

# -----------------------------
# Evaluation: SSR/PSR metrics when K_T == K
# -----------------------------
Evaluation <- function(beta, Beta_true, K_T) {
  IIndex <- (Beta_true != 0)
  index_beta <- (beta != 0)
  index_Beta <- (Beta_true != 0)
  inc <- abs(index_beta - index_Beta)  # logical->numeric in R gives 0/1
  
  ST <- numeric(K_T)
  sk <- numeric(K_T)
  ssr_each <- numeric(K_T)
  psr_eKK <- numeric(K_T)
  
  for (kt in seq_len(K_T)) {
    # mismatch count on true-nonzero positions
    ST[kt] <- sum(inc[IIndex[, kt], kt])
    sk[kt] <- sum(IIndex[, kt])
    ssr_each[kt] <- as.numeric(ST[kt] == 0)
    psr_eKK[kt] <- (sk[kt] - ST[kt]) / sk[kt]
  }
  
  SSR_eKK <- matrix(ssr_each, ncol = 1)
  PSR_eKK <- matrix(psr_eKK, ncol = 1)
  ssr <- as.numeric(sum(ssr_each) == K_T)
  psr <- as.numeric(sum(sk * psr_eKK) / sum(sk))
  
  return(list(ssr = ssr, psr = psr, SSR_eKK = SSR_eKK, PSR_eKK = PSR_eKK))
}

# -----------------------------
# Evaluation_union: union-based SSR/PSR when K_T != K (mis/over-spec)
# -----------------------------
Evaluation_union <- function(beta, Beta_true, K_T, K) {
  IIndex <- (Beta_true != 0)
  Index <- integer(0)
  
  for (ss in seq_len(K_T)) {
    ttmp <- which(IIndex[, ss] != 0)
    Index <- union(Index, ttmp)
  }
  
  index_beta <- integer(0)
  for (ks in seq_len(K)) {
    Ind <- which(beta[, ks] != 0)
    index_beta <- union(index_beta, Ind)
  }
  
  ssr <- all(Index %in% index_beta)
  psr <- length(intersect(Index, index_beta)) / length(Index)
  
  return(list(ssr = as.numeric(ssr), psr = as.numeric(psr)))
}

# -----------------------------
# predict: mixture prediction (weighted average of linear predictors)
# -----------------------------
predict_sEAM <- function(X, beta, Pi, K) {
  N <- nrow(X)
  pred <- rep(0, N)
  for (k in seq_len(K)) {
    yy <- as.numeric(X %*% beta[, k])
    pred <- pred + Pi[k] * yy
  }
  return(pred)
}

# -----------------------------
# sEAM_main: translated from your MATLAB sEAM_main
# -----------------------------
sEAM_main <- function(X, Y, OPTS) {
  stopifnot(is.matrix(X), is.numeric(Y))
  N <- nrow(X); p <- ncol(X)
  t0 <- proc.time()[["elapsed"]]
  
  K <- OPTS$K
  S <- OPTS$S
  maxit <- OPTS$maxit
  maxitLasso <- OPTS$maxitLasso
  LLHP <- isTRUE(OPTS$LLHP) || identical(OPTS$LLHP, "True")
  
  # initialization
  beta <- Lassoinit(X, Y, K, S, maxitLasso)
  
  # E/A/M init
  nu_square <- stats::var(Y)
  sigma_hat <- sqrt(nu_square) / K
  sigma <- rep(sigma_hat, K)
  u <- 1000
  pi <- rep(1 / K, K)
  r <- matrix(0, nrow = N, ncol = K)
  
  beta_path  <- matrix(0, nrow = p, ncol = K * (maxit + 1))
  sigma_path <- matrix(0, nrow = K, ncol = (maxit + 1))
  pi_path    <- matrix(0, nrow = K, ncol = (maxit + 1))
  u_path     <- numeric(maxit + 1)
  
  LLLH <- NULL
  if (LLHP) {
    LLLH <- numeric(maxit + 1)
    temp <- rep(0, N)
    psi <- 0
    for (k in seq_len(K)) {
      mu_k <- as.numeric(X %*% beta[, k])
      temp <- temp + pi[k] * stats::dnorm(Y, mean = mu_k, sd = sigma[k])
      psi <- psi + (1 / N) * (nu_square / (sigma[k]^2) + log((sigma[k]^2) / (nu_square^2)) - 1)
    }
    LLLH[1] <- sum(log(pmax(temp, .Machine$double.eps))) - N * psi
  }
  
  beta_path[, 1:K] <- beta
  sigma_path[, 1] <- sigma
  pi_path[, 1] <- pi
  u_path[1] <- u
  
  for (t in seq_len(maxit)) {
    # E-step
    for (k in seq_len(K)) {
      mu_k <- as.numeric(X %*% beta[, k])
      r[, k] <- pi[k] * stats::dnorm(Y, mean = mu_k, sd = sigma[k])
    }
    r <- r / pmax(rowSums(r), .Machine$double.eps)
    
    # A/M-step per component
    for (k in seq_len(K)) {
      pi[k] <- sum(r[, k]) / N
      
      sigma_prev <- sigma[k]
      resid <- Y - as.numeric(X %*% beta[, k])
      
      sigma[k] <- sqrt((sum(r[, k] * (resid^2)) + 2 * nu_square) / (sum(r[, k]) + 2))
      
      g <- as.numeric(crossprod(X, (r[, k] * resid) / (sigma_prev^2)))
      tmp <- beta[, k] + (1 / (u * N)) * g
      
      ord <- order(abs(tmp), decreasing = TRUE)
      if (S[k] < p) tmp[ord[(S[k] + 1):p]] <- 0
      beta[, k] <- tmp
    }
    
    # likelihood + adapt u
    if (LLHP) {
      temp <- rep(0, N)
      psi <- 0
      for (k in seq_len(K)) {
        mu_k <- as.numeric(X %*% beta[, k])
        temp <- temp + pi[k] * stats::dnorm(Y, mean = mu_k, sd = sigma[k])
        psi <- psi + (1 / N) * (nu_square / (sigma[k]^2) + log((sigma[k]^2) / (nu_square^2)) - 1)
      }
      LLLH[t + 1] <- sum(log(pmax(temp, .Machine$double.eps))) - N * psi
      
      if (LLLH[t + 1] <= LLLH[t]) {
        u <- u * 1.2
      }
    }
    
    u_path[t + 1] <- u
    beta_path[, (K * t + 1):(K * t + K)] <- beta
    sigma_path[, t + 1] <- sigma
    pi_path[, t + 1] <- pi
  }
  
  out <- list(
    beta = beta,
    pi = pi,
    time = proc.time()[["elapsed"]] - t0,
    beta_path = beta_path,
    sigma_path = sigma_path,
    pi_path = pi_path,
    u_path = u_path
  )
  
  if (LLHP) {
    out$LLH <- matrix(LLLH / N, ncol = 1)
    out$LLLH <- 5 * matrix(LLLH, ncol = 1)
  } else {
    out$LLH <- NULL
    out$LLLH <- NULL
  }
  
  return(out)
}

simulate_fmr_scenario2 <- function(n = 300,
                                   p = 2000,
                                   balanced = TRUE,
                                   corr = c("indep", "ar1", "exchangeable"),
                                   rho = 0.3,
                                   include_intercept_col1 = FALSE,
                                   seed = 1) {
  corr <- match.arg(corr)
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install MASS: install.packages('MASS')")
  }
  
  set.seed(seed)
  
  # --- mixture settings from screenshot ---
  if (balanced) {
    pi <- c(0.5, 0.5)
    sigma <- c(3, 5)
  } else {
    pi <- c(0.3, 0.7)
    sigma <- c(3, 3)
  }
  K <- 2
  
  # --- covariance Sigma for X ---
  # Screenshot explicitly provides Correlation 1: independent.
  # I also provide ar1/exchangeable as common options (not shown in your screenshot).
  Sigma <- switch(
    corr,
    indep = diag(p),
    ar1 = {
      idx <- seq_len(p)
      outer(idx, idx, function(i, j) rho^abs(i - j))
    },
    exchangeable = {
      # corr(xj, xl)=rho for j!=l, var=1
      m <- matrix(rho, nrow = p, ncol = p)
      diag(m) <- 1
      m
    }
  )
  
  # --- generate X ~ N(0, Sigma) ---
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # optional: treat first column as intercept (x1 = 1)
  if (include_intercept_col1) {
    X[, 1] <- 1
  }
  
  # --- true beta matrices for Scenario 2 ---
  Beta_true <- matrix(0, nrow = p, ncol = K)
  
  M1 <- c(1, 3, 5, 10)
  b1 <- c(3.7, -2.9, 3.2, -2.2)
  
  M2 <- c(1, 2, 6, 7, 8, 25)
  b2 <- c(3.2, 3.3, 4.5, 2.3, -4.0, 2.7)
  
  Beta_true[M1, 1] <- b1
  Beta_true[M2, 2] <- b2
  
  # --- latent class z ~ Multinomial(pi) ---
  z <- sample.int(K, size = n, replace = TRUE, prob = pi)  # values in {1,2}
  
  # --- generate Y: y = x^T beta_z + eps, eps ~ N(0, sigma_z^2) ---
  mu <- numeric(n)
  for (i in seq_len(n)) {
    mu[i] <- sum(X[i, ] * Beta_true[, z[i]])
  }
  Y <- numeric(n)
  idx1 <- which(z == 1)
  idx2 <- which(z == 2)
  
  Y[idx1] <- mu[idx1] + rnorm(length(idx1), 0, sigma[1])
  Y[idx2] <- mu[idx2] + rnorm(length(idx2), 0, sigma[2])
  
  return(list(
    X = X,
    Y = Y,
    z = z,                 # true component labels (hidden in practice)
    Beta_true = Beta_true, # p x 2
    sigma = sigma
  ))
}




# -----------------------------
# Generate AR(1) covariance Sigma of size (p+1)x(p+1)
# Sigma_ij = rho^{|i-j|}
# -----------------------------
ar1_cov <- function(dim, rho) {
  toeplitz(rho^(0:(dim - 1)))
}


# -----------------------------
# Example p=1000
# -----------------------------
set.seed(1)

n <- 200; p <- 1000; rho <- 0.8; tau <- 0.5
B <- 2   # number of replications

true_idx <- c(1, 2, 3, 4)

# Precompute covariance
Sigma <- ar1_cov(p + 1, rho)
R <- chol(Sigma)

# Store results
S_vec <- numeric(B)
min_model_vec <- numeric(B)
S_randomU_vec <- numeric(B)
Joint_vec <- numeric(B)

for (b in 1:B) {
  
  # 1) Generate (U', X)
  Z <- matrix(rnorm(n * (p + 1)), n, p + 1)
  data <- Z %*% R
  
  Uprime <- data[, 1]
  X <- data[, -1, drop = FALSE]
  
  # 2) U
  U <- pnorm(Uprime)
  
  # 3) Generate Y
  eps <- rnorm(n)
  
  sigma_x <- exp(0.5 * (X[, 5] + X[, 6] + X[, 7]))
  
  Y <- 2 * as.numeric(U > 0.4) * (X[, 1])^2 +
    (1 + U) * (X[, 2])^2 +
    (1 + U) * X[, 3] +
    (2 - 3 * U)^2 * X[, 4] +
    sigma_x * eps
  
  #### ===== conditional =====
  est <- estimate_qcor_tau_p(Y, X, U, tau = tau)
  score <- est$omega_star
  ord <- order(score, decreasing = TRUE)
  rank_pos <- integer(length(score))
  rank_pos[ord] <- seq_along(ord)
  
  S_vec[b] <- max(rank_pos[true_idx])
  
  #### ===== our =====
  min_model_vec[b] <- qc_screen_with_min_model_size(
    X, Y, tau, active_index = true_idx
  )$min_model_size
  
  #### ===== Joint =====
  S_list = seq(4,100,4)
  for (s in S_list) {
    OPTS <- list(
      K = 1,
      S = s,        # You can set how many features to keep for each component
      maxit = 50,
      maxitLasso = 500,
      LLHP = TRUE,
      K_T = 2
    )
    
    res <- sEAM_main(X, Y, OPTS)
    
    selected_idx_by_component <- function(beta) {
      K <- ncol(beta)
      lapply(seq_len(K), function(k) which(beta[, k] != 0))
    }
    
    # Usage:
    # res <- sEAM_main(X, Y, OPTS)
    idx_list <- selected_idx_by_component(res$beta)
    if(all(true_idx %in% unlist(idx_list))){
      break
    }
  }
  Joint_vec[b] <- s
  
  print(b)
}

result_table <- data.frame(
  Method = c("CQD", 
             "Conditional copula function", 
             "Joint"),
  Mean = c(mean(min_model_vec),
           mean(S_vec),
           mean(Joint_vec)),
  Median = c(median(min_model_vec),
             median(S_vec),
             median(Joint_vec)),
  SD = c(sd(min_model_vec),
        sd(S_vec),
        sd(Joint_vec))
)

result_table
write.csv(result_table, "Table4_p1000.csv", row.names = FALSE)

# -----------------------------
# Example p=2000
# -----------------------------
set.seed(1)

n <- 200; p <- 2000; rho <- 0.8; tau <- 0.5
B <- 2   # number of replications

true_idx <- c(1, 2, 3, 4)

# Precompute covariance
Sigma <- ar1_cov(p + 1, rho)
R <- chol(Sigma)

# Store results
S_vec <- numeric(B)
min_model_vec <- numeric(B)
S_randomU_vec <- numeric(B)
Joint_vec <- numeric(B)

for (b in 1:B) {
  
  # 1) Generate (U', X)
  Z <- matrix(rnorm(n * (p + 1)), n, p + 1)
  data <- Z %*% R
  
  Uprime <- data[, 1]
  X <- data[, -1, drop = FALSE]
  
  # 2) U
  U <- pnorm(Uprime)
  
  # 3) Generate Y
  eps <- rnorm(n)
  
  sigma_x <- exp(0.5 * (X[, 5] + X[, 6] + X[, 7]))
  
  Y <- 2 * as.numeric(U > 0.4) * (X[, 1])^2 +
    (1 + U) * (X[, 2])^2 +
    (1 + U) * X[, 3] +
    (2 - 3 * U)^2 * X[, 4] +
    sigma_x * eps
  
  #### ===== conditional =====
  est <- estimate_qcor_tau_p(Y, X, U, tau = tau)
  score <- est$omega_star
  ord <- order(score, decreasing = TRUE)
  rank_pos <- integer(length(score))
  rank_pos[ord] <- seq_along(ord)
  
  S_vec[b] <- max(rank_pos[true_idx])
  
  #### ===== our =====
  min_model_vec[b] <- qc_screen_with_min_model_size(
    X, Y, tau, active_index = true_idx
  )$min_model_size
  
  #### ===== Joint =====
  S_list = seq(4,100,4)
  for (s in S_list) {
    OPTS <- list(
      K = 1,
      S = s,        # You can set how many features to keep for each component
      maxit = 50,
      maxitLasso = 500,
      LLHP = TRUE,
      K_T = 2
    )
    
    res <- sEAM_main(X, Y, OPTS)
    
    selected_idx_by_component <- function(beta) {
      K <- ncol(beta)
      lapply(seq_len(K), function(k) which(beta[, k] != 0))
    }
    
    # Usage:
    # res <- sEAM_main(X, Y, OPTS)
    idx_list <- selected_idx_by_component(res$beta)
    if(all(true_idx %in% unlist(idx_list))){
      break
    }
  }
  Joint_vec[b] <- s
  
  print(b)
}

result_table <- data.frame(
  Method = c("CQD", 
             "Conditional copula function", 
             "Joint"),
  Mean = c(mean(min_model_vec),
           mean(S_vec),
           mean(Joint_vec)),
  Median = c(median(min_model_vec),
             median(S_vec),
             median(Joint_vec)),
  SD = c(sd(min_model_vec),
        sd(S_vec),
        sd(Joint_vec))
)

result_table