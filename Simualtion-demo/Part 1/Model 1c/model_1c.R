# Required packages
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

# Check and install missing packages
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

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

# (2.8) Copula Quantile Dependence  —— vectorized version
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
  t_vec <- qnorm((n / (n + 1)) * F_Y(Y))  # transformed Y
  
  # Parallel setup
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


############## compared methods ################
QaSIS <- function(Y, X, tau = 0.5, active_index = NULL, n_workers = max(1, parallel::detectCores() - 5)) {
  p <- ncol(X)
  n <- length(Y)
  
  # Center Y at quantile tau (consistent with the original alpha logic)
  Y <- Y - as.numeric(quantile(Y, tau))
  
  # Standardize X (column mean 0 and variance 1)
  X <- scale(X)
  
  # Parallel setup
  cl <- parallel::makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  
  # Parallel computation of the score for each j
  s_vec <- foreach::foreach(
    j = 1:p,
    .combine = "c",
    .packages = c("splines", "quantreg")
  ) %dopar% {
    x0 <- X[, j]
    knots <- quantile(x0, c(1 / 3, 2 / 3))
    a <- splines::bs(x0, knots = knots, degree = 1)
    fit <- quantreg::rq(Y ~ a, tau = tau)
    sum(fitted(fit)^2) / n
  }
  
  stopCluster(cl)
  
  ranked_index <- order(s_vec, decreasing = TRUE)
  
  if (!is.null(active_index)) {
    rank_positions <- match(active_index, ranked_index)
    min_model_size <- max(rank_positions, na.rm = TRUE)
  } else {
    min_model_size <- NA_integer_
  }
  
  return(min_model_size)
}

QcSIS <- function(Y, X, alpha = 0.5, active_index = NULL) {
  n <- length(Y)
  s_vec=QCSIS(X, Y, alpha, n-1)$w
  
  ranked_index <- order(s_vec, decreasing = TRUE)
  
  if (!is.null(active_index)) {
    rank_positions <- match(active_index, ranked_index)
    min_model_size <- max(rank_positions, na.rm = TRUE)
  } else {
    min_model_size <- NA_integer_
  }
  
  return(min_model_size)
}

wdGscreening <- function(Y, X, ppp = 1,
                         n_workers = max(1, parallel::detectCores() - 5)) {
  Y <- as.matrix(Y)
  n <- nrow(X)
  
  Gaussian_temp <- qnorm((1:n) / (n + 1))
  gt <- function(sq) Gaussian_temp[sq]
  
  ry <- Rfast::Rank(Y)
  rx <- colRanks(as.matrix(X))  # removed parallel=TRUE
  gry <- gt(ry)
  grx <- apply(rx, 2, gt)
  
  X <- as.matrix(grx)
  Y <- gry
  
  G <- Rfast::rmvnorm(n, mu = rep(0, 2), sigma = diag(2))
  
  # Parallel setup
  cl <- parallel::makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  
  # Parallel computation
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

PCSIS <- function(Y, X, active_index = NULL,
                  tau = 0.75,
                  n_workers = max(1, parallel::detectCores() - 5)) {
  stopifnot(is.numeric(Y), is.matrix(X), length(Y) == nrow(X))
  n <- length(Y); p <- ncol(X)
  
  # 1) Preprocessing: center Y and standardize X (with simple protection for zero-variance columns)
  Y <- Y - as.numeric(stats::quantile(Y, tau, names = FALSE))
  mu  <- colMeans(X)
  sdv <- apply(X, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
  X <- sweep(sweep(X, 2, mu, `-`), 2, sdv, `/`)
  
  # 2) Parallel: create the cluster only once
  cl <- parallel::makeCluster(max(1, n_workers))
  doParallel::registerDoParallel(cl)
  on.exit({
    foreach::registerDoSEQ()
    parallel::stopCluster(cl)
  }, add = TRUE)
  
  # 3) First scoring: spline + quantile regression (if it fails, assign score 0 to avoid aborting everything)
  s_vec <- foreach::foreach(
    j = 1:p, .combine = "c",
    .packages = c("splines", "quantreg")
  ) %dopar% {
    x0 <- X[, j]
    # Knots chosen at 1/3 and 2/3 quantiles; if degenerate, use fewer or none
    qs <- unique(stats::quantile(x0, c(1/3, 2/3), names = FALSE, type = 7))
    a  <- tryCatch(splines::bs(x0, knots = qs, degree = 1), error = function(e) NULL)
    if (is.null(a)) return(0)
    fit <- tryCatch(quantreg::rq(Y ~ a, tau = tau), error = function(e) NULL)
    if (is.null(fit)) return(0)
    sum(fitted(fit)^2) / n
  }
  
  # 4) Threshold and candidate set
  if (is.null(active_index) || length(active_index) == 0) {
    threshold <- -Inf
  } else {
    ai <- intersect(unique(as.integer(active_index)), seq_len(p))
    threshold <- suppressWarnings(min(s_vec[ai], na.rm = TRUE))
    if (!is.finite(threshold)) threshold <- -Inf
  }
  idx <- which(s_vec >= threshold)
  if (length(idx) == 0L) return(0L)
  
  # 5) Second scoring: pcov (record 0 if it fails)
  res <- foreach::foreach(
    i = idx, .combine = "c",
    .export = "pcov"
  ) %dopar% {
    val <- tryCatch(pcov(X[, i], Y), error = function(e) 0)
    if (!is.finite(val)) 0 else val
  }
  
  ranked_index   <- order(res, decreasing = TRUE)
  original_index <- idx[ranked_index]
  
  # 6) Minimum model size
  if (!is.null(active_index) && length(active_index)) {
    pos <- match(active_index, original_index)
    if (all(is.na(pos))) return(NA_integer_)
    return(max(pos, na.rm = TRUE))
  } else {
    return(NA_integer_)
  }
}


############## compared methods -- KCQD ################
u4_stat_groups_parallel <- function(X, Y, tau, beta, groups = NULL,
                                    n_workers = max(1, parallel::detectCores() - 5)) {
  n <- length(Y)
  
  ## --------- Internal helper functions ----------
  # d(Y,beta): 1D case -> sign matrix, diag = 0
  d_matrix_1d <- function(y, beta, eps = 1e-12) {
    u <- beta - y
    u_norm <- u / pmax(abs(u), eps)
    D <- outer(u_norm, u_norm, "*")
    diag(D) <- 0
    D
  }
  # C based on multivariate Gaussian kernel: rows are samples i=1..n, columns are features in the group
  c_matrix_vec <- function(Xg) {
    # pairwise squared distances between row-vectors
    G <- tcrossprod(Xg)                       # n x n
    r2 <- rowSums(Xg^2)
    sq <- outer(r2, r2, "+") - 2 * G         # ||xi-xj||^2
    # median heuristic for sigma (take square root)
    sigma <- sqrt(median(sq[upper.tri(sq)]))
    K <- exp(-sq / (2 * sigma^2))
    diag(K) <- 0
    K
  }
  perms4 <- matrix(c(
    1,2,3,4, 1,2,4,3, 1,3,2,4, 1,3,4,2, 1,4,2,3, 1,4,3,2,
    2,1,3,4, 2,1,4,3, 2,3,1,4, 2,3,4,1, 2,4,1,3, 2,4,3,1,
    3,1,2,4, 3,1,4,2, 3,2,1,4, 3,2,4,1, 3,4,1,2, 3,4,2,1,
    4,1,2,3, 4,1,3,2, 4,2,1,3, 4,2,3,1, 4,3,1,2, 4,3,2,1
  ), ncol = 4, byrow = TRUE)
  
  u4_from_C_D_mean <- function(C, D, comb4) {
    total <- 0
    for (k in seq_len(ncol(comb4))) {
      ids <- comb4[, k]
      s_perm <- 0
      for (r in 1:nrow(perms4)) {
        m <- ids[perms4[r, ]]
        u <- m[1]; v <- m[2]; q <- m[3]; r4 <- m[4]
        s_perm <- s_perm + C[u, v] * ( D[u, v] - D[u, q] - D[v, r4] + D[q, r4] )
      }
      total <- total + s_perm / 24
    }
    total / ncol(comb4)
  }
  
  ## --------- Pre-computation (shared by all groups) ----------
  D <- d_matrix_1d(Y, beta)
  comb4 <- combn(n, 4)
  
  ## --------- Group definition ----------
  p <- ncol(X)
  if (is.null(groups)) {
    # By default, split columns into 50 groups in order (20 columns per group)
    stopifnot(p %% 50 == 0)
    gsize <- p / 50
    groups <- split(seq_len(p), rep(1:50, each = gsize))
  } else {
    # User-defined: a list where each element is a vector of column indices for one group
    stopifnot(is.list(groups))
  }
  
  ## --------- Parallel computation for each group ----------
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  
  vals <- foreach::foreach(g = seq_along(groups), .combine = c,
                           .packages = character()) %dopar% {
                             Xg <- X[, groups[[g]], drop = FALSE]   # n x d (here d = 20)
                             C <- c_matrix_vec(Xg)
                             u4_from_C_D_mean(C, D, comb4)
                           }
  as.numeric(vals)  # length = number of groups (default 50)
}


u4_stat_rowgroups_mean <- function(X, Y, tau,
                                   n_row_groups = 20,
                                   n_col_groups = 50,
                                   n_workers_col = max(1, parallel::detectCores() - 5)) {
  n <- nrow(X); p <- ncol(X)
  beta <- as.numeric(quantile(Y, (tau + 1) / 2, type = 7))
  
  # Generate row groups and column groups
  rows_per_group <- n / n_row_groups
  row_groups <- split(seq_len(n), rep(seq_len(n_row_groups), each = rows_per_group))
  
  cols_per_group <- p / n_col_groups
  col_groups <- split(seq_len(p), rep(seq_len(n_col_groups), each = cols_per_group))
  
  # For each row group, run once and extract u4 results for its 50 column groups
  res_mat <- matrix(NA_real_, nrow = n_row_groups, ncol = n_col_groups)
  for (g in seq_along(row_groups)) {
    idx <- row_groups[[g]]
    Xg  <- X[idx, , drop = FALSE]
    Yg  <- Y[idx]
    
    # Inner parallelization is over column groups (avoid nested over-subscription)
    out <- u4_stat_groups_parallel(
      Xg, Yg, tau, beta,
      groups = col_groups,          # explicitly pass the 50 column groups
      n_workers = n_workers_col     # number of workers for column-group parallelism
    )
    # out is a numeric vector of length equal to number of column groups
    res_mat[g, ] <- out
  }
  
  # Take column-wise average over 20 row groups to obtain a 50×1 vector
  colMeans(res_mat)
}

select_groups_by_t <- function(T_vals, u, alpha_G, include_equal = TRUE) {
  T_vals <- as.numeric(T_vals)
  
  # Compute t_ne
  t_candidates <- sort(unique(abs(T_vals)))
  t_ne <- NA_real_
  for (t in t_candidates) {
    num <- u + sum(T_vals <= -t)
    den <- sum(T_vals >=  t)
    ratio <- if (den > 0) num / den else Inf
    if (ratio <= alpha_G) { t_ne <- t; break }
  }
  
  # Select groups according to t_ne
  if (is.na(t_ne)) {
    idx <- integer(0)
  } else {
    if (include_equal) {
      idx <- which(T_vals >= t_ne)
    } else {
      idx <- which(T_vals >  t_ne)
    }
  }
  
  idx
}


KCQD=function(Y, X, alpha, active_index){
  X1=X[1:200,]
  Y1=Y[1:200]
  X2=X[201:400,]
  Y2=Y[201:400]
  T_vals=u4_stat_rowgroups_mean(X1, Y1, alpha)
  cols_per_group = ncol(X) / 50
  
  group_idx=as.vector(select_groups_by_t(T_vals, 0, 0.2))
  element_position=c()
  for(i in group_idx){
    element_position=c(element_position,((i-1)*cols_per_group+1):(i*cols_per_group))
  }
  
  s_vec=u4_stat_rowgroups_mean(X2[,element_position], Y2, alpha,n_col_groups=length(element_position))
  
  ranked_index <- order(s_vec, decreasing = TRUE)
  
  if (!is.null(active_index)) {
    rank_positions <- match(active_index, ranked_index)
    min_model_size <- max(rank_positions, na.rm = TRUE)
  } else {
    min_model_size <- NA_integer_
  }
  
  if(min_model_size==-Inf){
    return(ncol(X))
  }
  
  return(min_model_size)
}


############## simulation 1c p=1000 alpha=0.5 ##############
generate_model_1c <- function(n = 400, p = 1000, alpha = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- outer(1:p, 1:p, function(i, j) 0.8^abs(i - j))
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install 'MASS' package: install.packages('MASS')")
  }
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  X1 <- X[, 1]; X2 <- X[, 2]; X3 <- X[, 3]
  X4 <- X[, 4]; X5 <- X[, 5]; X6 <- X[, 6]; X7 <- X[, 7]
  
  q_alpha <- qnorm(alpha)
  epsilon <- rnorm(n, mean = -q_alpha, sd = 1)
  
  Y <- 2 * (X1^2) +
    2 * (X2^2) +
    exp(5 * X3 + 5 * X4 + 5 * X5 + 5 * X6 + 5 * X7) * epsilon
  list(X = X, Y = Y)
}

set.seed(1234)
n=400
p=1000
rep=2
alpha = 0.5
active.set <- c(1, 2)                      # user-specified active set

# Store results
fs_result = c()
sis_result = c()
dc_sis_result = c()
qa_result = c()
qc_result = c()
wd_result = c()
pc_result = c()
KCQD_result = c()

# Store running times
fs_time = c()
sis_time = c()
dc_sis_time = c()
qa_time = c()
qc_time = c()
wd_time = c()
pc_time = c()
KCQD_time = c()

for (i in 1:rep) {
  data <- generate_model_1c(n, p)
  X <- data$X
  Y <- data$Y
  
  # FS
  t0 <- proc.time()
  result <- qc_screen_with_min_model_size(X, Y, alpha, active_index = active.set)
  fs_time <- c(fs_time, (proc.time() - t0)[["elapsed"]])
  fs_result <- c(fs_result, result$min_model_size)
  
  # SIS
  t0 <- proc.time()
  sis_result <- c(sis_result, max(Rfast::Rank(abs(cor(X, Y, method = 'pearson')), descending = TRUE)[active.set]))
  sis_time <- c(sis_time, (proc.time() - t0)[["elapsed"]])
  
  # DC-SIS
  t0 <- proc.time()
  dc_sis_result <- c(dc_sis_result, max(Rfast::Rank(mdcor(Y, X, 'V'), descending = TRUE)[active.set]))
  dc_sis_time <- c(dc_sis_time, (proc.time() - t0)[["elapsed"]])
  
  # QaSIS
  t0 <- proc.time()
  qa_result <- c(qa_result, QaSIS(Y, X, tau = alpha, active.set))
  qa_time <- c(qa_time, (proc.time() - t0)[["elapsed"]])
  
  # QcSIS
  t0 <- proc.time()
  qc_result <- c(qc_result, QcSIS(Y, X, alpha, active.set))
  qc_time <- c(qc_time, (proc.time() - t0)[["elapsed"]])
  
  # WD
  t0 <- proc.time()
  wd_result <- c(wd_result, max(Rfast::Rank(wdGscreening(Y, X, ppp = 1), descending = TRUE)[active.set]))
  wd_time <- c(wd_time, (proc.time() - t0)[["elapsed"]])
  
  # PC
  t0 <- proc.time()
  pc_result <- c(pc_result, PCSIS(Y, X, active.set))
  pc_time <- c(pc_time, (proc.time() - t0)[["elapsed"]])
  
  # KCQD
  t0 <- proc.time()
  KCQD_result <- c(KCQD_result, KCQD(Y, X, alpha, active.set))
  KCQD_time <- c(KCQD_time, (proc.time() - t0)[["elapsed"]])
  
  print(i)
}

# results
vars <- list(
  fs_result = fs_result,
  sis_result = sis_result,
  dc_sis_result = dc_sis_result,
  qa_result = qa_result,
  qc_result = qc_result,
  KCQD_result = KCQD_result,
  pc_result = pc_result,
  wd_result = wd_result
)

# time
times <- list(
  fs_time = fs_time,
  sis_time = sis_time,
  dc_sis_time = dc_sis_time,
  qa_time = qa_time,
  qc_time = qc_time,
  KCQD_time = KCQD_time,
  pc_time = pc_time,
  wd_time = wd_time
)

# Output in "median (IQR)" format
result_table <- sapply(vars, function(x) {
  sprintf("%.2f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
})

time_table <- sapply(times, function(x) {
  sprintf("%.2f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
})

result_table
time_table

############## simulation 1c p=1000 alpha=0.75 ##############
generate_model_1c <- function(n = 400, p = 1000, alpha = 0.75, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- outer(1:p, 1:p, function(i, j) 0.8^abs(i - j))
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install 'MASS' package: install.packages('MASS')")
  }
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  X1 <- X[, 1]; X2 <- X[, 2]; X3 <- X[, 3]
  X4 <- X[, 4]; X5 <- X[, 5]; X6 <- X[, 6]; X7 <- X[, 7]
  
  q_alpha <- qnorm(alpha)
  epsilon <- rnorm(n, mean = -q_alpha, sd = 1)
  
  Y <- 2 * (X1^2) +
    2 * (X2^2) +
    exp(5 * X3 + 5 * X4 + 5 * X5 + 5 * X6 + 5 * X7) * epsilon
  list(X = X, Y = Y)
}

set.seed(1234)
n=400
p=1000
rep=2
alpha = 0.75
active.set <- c(1, 2)                      # user-specified active set

# Store results
fs_result = c()
sis_result = c()
dc_sis_result = c()
qa_result = c()
qc_result = c()
wd_result = c()
pc_result = c()
KCQD_result = c()

# Store running times
fs_time = c()
sis_time = c()
dc_sis_time = c()
qa_time = c()
qc_time = c()
wd_time = c()
pc_time = c()
KCQD_time = c()

for (i in 1:rep) {
  data <- generate_model_1c(n, p)
  X <- data$X
  Y <- data$Y
  
  # FS
  t0 <- proc.time()
  result <- qc_screen_with_min_model_size(X, Y, alpha, active_index = active.set)
  fs_time <- c(fs_time, (proc.time() - t0)[["elapsed"]])
  fs_result <- c(fs_result, result$min_model_size)
  
  # SIS
  t0 <- proc.time()
  sis_result <- c(sis_result, max(Rfast::Rank(abs(cor(X, Y, method = 'pearson')), descending = TRUE)[active.set]))
  sis_time <- c(sis_time, (proc.time() - t0)[["elapsed"]])
  
  # DC-SIS
  t0 <- proc.time()
  dc_sis_result <- c(dc_sis_result, max(Rfast::Rank(mdcor(Y, X, 'V'), descending = TRUE)[active.set]))
  dc_sis_time <- c(dc_sis_time, (proc.time() - t0)[["elapsed"]])
  
  # QaSIS
  t0 <- proc.time()
  qa_result <- c(qa_result, QaSIS(Y, X, tau = alpha, active.set))
  qa_time <- c(qa_time, (proc.time() - t0)[["elapsed"]])
  
  # QcSIS
  t0 <- proc.time()
  qc_result <- c(qc_result, QcSIS(Y, X, alpha, active.set))
  qc_time <- c(qc_time, (proc.time() - t0)[["elapsed"]])
  
  # WD
  t0 <- proc.time()
  wd_result <- c(wd_result, max(Rfast::Rank(wdGscreening(Y, X, ppp = 1), descending = TRUE)[active.set]))
  wd_time <- c(wd_time, (proc.time() - t0)[["elapsed"]])
  
  # PC
  t0 <- proc.time()
  pc_result <- c(pc_result, PCSIS(Y, X, active.set))
  pc_time <- c(pc_time, (proc.time() - t0)[["elapsed"]])
  
  # KCQD
  t0 <- proc.time()
  KCQD_result <- c(KCQD_result, KCQD(Y, X, alpha, active.set))
  KCQD_time <- c(KCQD_time, (proc.time() - t0)[["elapsed"]])
  
  print(i)
}

# results
vars <- list(
  fs_result = fs_result,
  sis_result = sis_result,
  dc_sis_result = dc_sis_result,
  qa_result = qa_result,
  qc_result = qc_result,
  KCQD_result = KCQD_result,
  pc_result = pc_result,
  wd_result = wd_result
)

# time
times <- list(
  fs_time = fs_time,
  sis_time = sis_time,
  dc_sis_time = dc_sis_time,
  qa_time = qa_time,
  qc_time = qc_time,
  KCQD_time = KCQD_time,
  pc_time = pc_time,
  wd_time = wd_time
)

# Output in "median (IQR)" format
result_table <- sapply(vars, function(x) {
  sprintf("%.2f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
})

time_table <- sapply(times, function(x) {
  sprintf("%.2f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
})

result_table
time_table


############## simulation 1c p=5000 alpha=0.5 ##############
generate_model_1c <- function(n = 400, p = 5000, alpha = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- outer(1:p, 1:p, function(i, j) 0.8^abs(i - j))
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install 'MASS' package: install.packages('MASS')")
  }
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  X1 <- X[, 1]; X2 <- X[, 2]; X3 <- X[, 3]
  X4 <- X[, 4]; X5 <- X[, 5]; X6 <- X[, 6]; X7 <- X[, 7]
  
  q_alpha <- qnorm(alpha)
  epsilon <- rnorm(n, mean = -q_alpha, sd = 1)
  
  Y <- 2 * (X1^2) +
    2 * (X2^2) +
    exp(5 * X3 + 5 * X4 + 5 * X5 + 5 * X6 + 5 * X7) * epsilon
  list(X = X, Y = Y)
}

set.seed(1234)
n=400
p=5000
rep=2
alpha = 0.5
active.set <- c(1, 2)                      # user-specified active set

# Store results
fs_result = c()
sis_result = c()
dc_sis_result = c()
qa_result = c()
qc_result = c()
wd_result = c()
pc_result = c()
KCQD_result = c()

# Store running times
fs_time = c()
sis_time = c()
dc_sis_time = c()
qa_time = c()
qc_time = c()
wd_time = c()
pc_time = c()
KCQD_time = c()

for (i in 1:rep) {
  data <- generate_model_1c(n, p)
  X <- data$X
  Y <- data$Y
  
  # FS
  t0 <- proc.time()
  result <- qc_screen_with_min_model_size(X, Y, alpha, active_index = active.set)
  fs_time <- c(fs_time, (proc.time() - t0)[["elapsed"]])
  fs_result <- c(fs_result, result$min_model_size)
  
  # SIS
  t0 <- proc.time()
  sis_result <- c(sis_result, max(Rfast::Rank(abs(cor(X, Y, method = 'pearson')), descending = TRUE)[active.set]))
  sis_time <- c(sis_time, (proc.time() - t0)[["elapsed"]])
  
  # DC-SIS
  t0 <- proc.time()
  dc_sis_result <- c(dc_sis_result, max(Rfast::Rank(mdcor(Y, X, 'V'), descending = TRUE)[active.set]))
  dc_sis_time <- c(dc_sis_time, (proc.time() - t0)[["elapsed"]])
  
  # QaSIS
  t0 <- proc.time()
  qa_result <- c(qa_result, QaSIS(Y, X, tau = alpha, active.set))
  qa_time <- c(qa_time, (proc.time() - t0)[["elapsed"]])
  
  # QcSIS
  t0 <- proc.time()
  qc_result <- c(qc_result, QcSIS(Y, X, alpha, active.set))
  qc_time <- c(qc_time, (proc.time() - t0)[["elapsed"]])
  
  # WD
  t0 <- proc.time()
  wd_result <- c(wd_result, max(Rfast::Rank(wdGscreening(Y, X, ppp = 1), descending = TRUE)[active.set]))
  wd_time <- c(wd_time, (proc.time() - t0)[["elapsed"]])
  
  # PC
  t0 <- proc.time()
  pc_result <- c(pc_result, PCSIS(Y, X, active.set))
  pc_time <- c(pc_time, (proc.time() - t0)[["elapsed"]])
  
  # KCQD
  t0 <- proc.time()
  KCQD_result <- c(KCQD_result, KCQD(Y, X, alpha, active.set))
  KCQD_time <- c(KCQD_time, (proc.time() - t0)[["elapsed"]])
  
  print(i)
}

# results
vars <- list(
  fs_result = fs_result,
  sis_result = sis_result,
  dc_sis_result = dc_sis_result,
  qa_result = qa_result,
  qc_result = qc_result,
  KCQD_result = KCQD_result,
  pc_result = pc_result,
  wd_result = wd_result
)

# time
times <- list(
  fs_time = fs_time,
  sis_time = sis_time,
  dc_sis_time = dc_sis_time,
  qa_time = qa_time,
  qc_time = qc_time,
  KCQD_time = KCQD_time,
  pc_time = pc_time,
  wd_time = wd_time
)

# Output in "median (IQR)" format
result_table <- sapply(vars, function(x) {
  sprintf("%.2f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
})

time_table <- sapply(times, function(x) {
  sprintf("%.2f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
})

result_table
time_table


############## simulation 1c p=5000 alpha=0.75 ##############
generate_model_1c <- function(n = 400, p = 5000, alpha = 0.75, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- outer(1:p, 1:p, function(i, j) 0.8^abs(i - j))
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install 'MASS' package: install.packages('MASS')")
  }
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  X1 <- X[, 1]; X2 <- X[, 2]; X3 <- X[, 3]
  X4 <- X[, 4]; X5 <- X[, 5]; X6 <- X[, 6]; X7 <- X[, 7]
  
  q_alpha <- qnorm(alpha)
  epsilon <- rnorm(n, mean = -q_alpha, sd = 1)
  
  Y <- 2 * (X1^2) +
    2 * (X2^2) +
    exp(5 * X3 + 5 * X4 + 5 * X5 + 5 * X6 + 5 * X7) * epsilon
  list(X = X, Y = Y)
}

set.seed(1234)
n=400
p=5000
rep=2
alpha = 0.75
active.set <- c(1, 2)                      # user-specified active set

# Store results
fs_result = c()
sis_result = c()
dc_sis_result = c()
qa_result = c()
qc_result = c()
wd_result = c()
pc_result = c()
KCQD_result = c()

# Store running times
fs_time = c()
sis_time = c()
dc_sis_time = c()
qa_time = c()
qc_time = c()
wd_time = c()
pc_time = c()
KCQD_time = c()

for (i in 1:rep) {
  data <- generate_model_1c(n, p)
  X <- data$X
  Y <- data$Y
  
  # FS
  t0 <- proc.time()
  result <- qc_screen_with_min_model_size(X, Y, alpha, active_index = active.set)
  fs_time <- c(fs_time, (proc.time() - t0)[["elapsed"]])
  fs_result <- c(fs_result, result$min_model_size)
  
  # SIS
  t0 <- proc.time()
  sis_result <- c(sis_result, max(Rfast::Rank(abs(cor(X, Y, method = 'pearson')), descending = TRUE)[active.set]))
  sis_time <- c(sis_time, (proc.time() - t0)[["elapsed"]])
  
  # DC-SIS
  t0 <- proc.time()
  dc_sis_result <- c(dc_sis_result, max(Rfast::Rank(mdcor(Y, X, 'V'), descending = TRUE)[active.set]))
  dc_sis_time <- c(dc_sis_time, (proc.time() - t0)[["elapsed"]])
  
  # QaSIS
  t0 <- proc.time()
  qa_result <- c(qa_result, QaSIS(Y, X, tau = alpha, active.set))
  qa_time <- c(qa_time, (proc.time() - t0)[["elapsed"]])
  
  # QcSIS
  t0 <- proc.time()
  qc_result <- c(qc_result, QcSIS(Y, X, alpha, active.set))
  qc_time <- c(qc_time, (proc.time() - t0)[["elapsed"]])
  
  # WD
  t0 <- proc.time()
  wd_result <- c(wd_result, max(Rfast::Rank(wdGscreening(Y, X, ppp = 1), descending = TRUE)[active.set]))
  wd_time <- c(wd_time, (proc.time() - t0)[["elapsed"]])
  
  # PC
  t0 <- proc.time()
  pc_result <- c(pc_result, PCSIS(Y, X, active.set))
  pc_time <- c(pc_time, (proc.time() - t0)[["elapsed"]])
  
  # KCQD
  t0 <- proc.time()
  KCQD_result <- c(KCQD_result, KCQD(Y, X, alpha, active.set))
  KCQD_time <- c(KCQD_time, (proc.time() - t0)[["elapsed"]])
  
  print(i)
}

# results
vars <- list(
  fs_result = fs_result,
  sis_result = sis_result,
  dc_sis_result = dc_sis_result,
  qa_result = qa_result,
  qc_result = qc_result,
  KCQD_result = KCQD_result,
  pc_result = pc_result,
  wd_result = wd_result
)

# time
times <- list(
  fs_time = fs_time,
  sis_time = sis_time,
  dc_sis_time = dc_sis_time,
  qa_time = qa_time,
  qc_time = qc_time,
  KCQD_time = KCQD_time,
  pc_time = pc_time,
  wd_time = wd_time
)

# Output in "median (IQR)" format
result_table <- sapply(vars, function(x) {
  sprintf("%.2f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
})

time_table <- sapply(times, function(x) {
  sprintf("%.2f (%.2f)", median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))
})

result_table
time_table
