# Required package list
pkgs <- c(
  "doParallel",
  "foreach",
  "dcov"
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

# φ'(Φ^{-1}(α)) and four related constants/variances 
phi_prime_at_alpha <- function(alpha) {
  z <- qnorm(alpha)
  -z * dnorm(z)                 # φ'(z) = -z φ(z)
}

# Return M_{⊥α,2}, M_{⊥α,4},  σ^2_{⊥α,1}, σ^2_{⊥α,3}
alpha_constants <- function(alpha) {
  stopifnot(alpha > 0, alpha < 1)
  phip <- phi_prime_at_alpha(alpha)
  
  # M_{⊥α,2}(ω) = (α - α^2)/(4π)
  M2 <- (alpha - alpha^2) / (4 * pi)
  
  # σ^2_{⊥α,3}(ω) = [α - α^2]^2 / (16 π^2)
  sig3_sq <- (alpha - alpha^2)^2 / (16 * pi^2)
  
  # M_{⊥α,4}(ω) = (2α^2 - 4α φ'(Φ^{-1}(α)) + 3[φ'(Φ^{-1}(α))]^2) / (24√3 π)
  M4 <- (2*alpha^2 - 4*alpha*phip + 3*(phip^2)) / (24 * sqrt(3) * pi)
  
  # σ^2_{⊥α,1}(ω)
  # = (α - α^2) * [ (18α^2 - 40α φ' + 25(φ')^2)/(400 π√5)  -  (3φ' - 2α)^2/(432 π^2) ]
  term1 <- (18*alpha^2 - 40*alpha*phip + 25*(phip^2)) / (400 * pi^2 * sqrt(5))
  term2 <- (3*phip - 2*alpha)^2 / (432 * pi^2)
  sig1_sq <- (alpha - alpha^2) * (term1 - term2)
  
  list(M2 = M2, M4 = M4, sig1_sq = sig1_sq, sig3_sq = sig3_sq)
}

# Use the existing estimate_D_star() directly
compute_D_dagger_from_alpha <- function(X, Y, h, alpha, omega=function(u) dnorm(qnorm(u))^2,
                                        delta_n_alpha = 0, eps = 1e-12, n_workers = max(1, parallel::detectCores() - 5)) {
  stopifnot(is.matrix(X), length(Y) == nrow(X))
  n <- length(Y); p <- ncol(X)
  
  # 1) First compute \hat D_α
  F_Y <- ecdf(Y)
  t_vec <- qnorm((n / (n + 1)) * F_Y(Y))  # Transformed Y
  
  # Set up parallel backend
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
  
  # 2) Constants (depending only on α)
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

# Input:
#   D_dagger: length-p vector of \hat D^{\dagger}_α
#   gamma   : γ
# Output: selection result, threshold, etc.
adaptive_threshold_select <- function(D_dagger, gammas) {
  gammas <- as.numeric(gammas)
  p <- length(D_dagger)
  
  # Preprocessing (independent of gamma)
  absD     <- abs(D_dagger)
  t_sorted <- sort(absD, decreasing = FALSE)
  t_all    <- c(0, t_sorted, Inf)   # t_0=0, t_{p+1}=Inf
  
  count_leq_neg_t <- function(t) sum(D_dagger <= -t)
  count_geq_pos_t <- function(t) sum(D_dagger >=  t)
  
  # For a given gamma, solve for tau/m/selected
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
  
  # For each gamma, compute results separately
  details_list <- lapply(gammas, solve_for_gamma)
  names(details_list) <- paste0("gamma=", format(gammas, trim = TRUE))
  
  # Summary table
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
# Classical Gram–Schmidt: orthonormalization on column vectors
orthonormalize <- function(X, eps = 1e-12) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  Y <- matrix(0, n, p)
  # First column
  v <- X[, 1, drop = FALSE]
  nv <- sqrt(sum(v^2))
  if (nv < eps) stop("Column 1 has near-zero norm.")
  Y[, 1] <- (v / nv)[, 1]
  
  # Remaining columns
  for (j in 2:p) {
    Yj <- Y[, 1:(j - 1), drop = FALSE]          # existing orthonormal basis
    xj <- X[, j, drop = FALSE]
    # w = t(Yj) %*% xj; projection = Yj %*% w
    w <- crossprod(Yj, xj)                       # (j-1) x 1
    proj <- Yj %*% w                             
    yj <- xj - proj
    ny <- sqrt(sum(yj^2))
    if (ny < eps) stop(paste("Column", j, "becomes near-zero after projection."))
    Y[, j] <- (yj / ny)[, 1]
  }
  Y
}

get_equi_features <- function(X, eps = 1e-8) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  
  # Column normalization
  scale <- sqrt(colSums(X^2))
  if (any(scale < eps)) stop("Some columns have near-zero norm.")
  Xstd <- sweep(X, 2, scale, "/")
  
  # Covariance (Gram) matrix and its inverse
  sigma <- crossprod(Xstd, Xstd)                 # p x p
  # Take the smallest eigenvalue (symmetric matrix)
  ev <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(ev)
  
  sigma_inv <- tryCatch(solve(sigma),
                        error = function(e) stop("sigma is singular or ill-conditioned."))
  
  # Set sj
  sj <- min(1.0, 2.0 * lambda_min)
  sj <- sj - 1e-5
  
  # A = 2*sj*I - sj^2 * sigma^{-1}
  A <- 2 * sj * diag(p) - (sj^2) * sigma_inv
  
  # C = chol(A) in R gives U s.t. A = t(U) %*% U; aligned with L^T in Python
  C <- tryCatch(chol(A),
                error = function(e) stop("Matrix A is not positive definite; chol failed."))
  
  # Construct U: orthonormalize [Xstd, N(0,1)] and take the last p columns
  Xn <- matrix(rnorm(n * p), n, p)
  XX <- cbind(Xstd, Xn)
  XXo <- orthonormalize(XX)
  U <- XXo[, (p + 1):(2 * p), drop = FALSE]
  
  # Xnew = Xstd %*% (I - sj * sigma_inv) + U %*% C
  Xnew <- Xstd %*% (diag(p) - sj * sigma_inv) + U %*% C
  Xnew
}

PCSIS <- function(Y, X, gammas, active_set,
                  n_workers = max(1, parallel::detectCores() - 5)) {
  # Allow passing scalar gamma
  if (missing(gammas)) stop("Please provide `gammas` (numeric vector or scalar).")
  gammas <- as.numeric(gammas)
  
  # 1) First screening step: select top 50 by dcscore
  dcscore <- mdcor(Y, X, "V")
  p_all   <- ncol(X)
  threshold<- min(dcscore[active_set])
  top_idx <- which(dcscore>=threshold)
  if (length(top_idx) == 0L) return(ncol(X))
  if (length(top_idx)<50) {
    # Sort and take the top 50 indices
    top_idx <- order(dcscore, decreasing = TRUE)[1:50]
    
    # If fewer than 50 indices, use the actual length
    top_idx <- top_idx[!is.na(top_idx)]
  }
  
  Xold <- X[, top_idx, drop = FALSE]
  Xnew <- get_equi_features(Xold)
  
  # 2) Parallel computing
  pold=apply(Xold, 2, function(x) pcov(x, Y))
  pnew=apply(Xnew, 2, function(x) pcov(x, Y))
  
  # 3) Compute res = pcov(Xold[,i],Y)^2 - pcov(Xnew[,i],Y)^2
  res <- pold^2-pnew^2
  
  # 4) Precompute threshold sequence (independent of gamma, can be reused)
  p        <- length(res)
  absD     <- abs(res)
  t_sorted <- sort(absD, decreasing = FALSE)
  t_all    <- c(0, t_sorted, Inf)  # t_0=0, t_{p+1}=Inf
  
  count_leq_neg_t <- function(t) sum(res <= -t)
  count_geq_pos_t <- function(t) sum(res >=  t)
  
  # For a given gamma, solve for tau, m, selected
  solve_for_gamma <- function(gamma) {
    if (!is.finite(gamma) || gamma <= 0) {
      return(list(gamma = gamma, tau = NA_real_, m = NA_integer_, selected_idx = integer(0)))
    }
    m <- 0L
    while (m <= p) {
      t_m <- t_all[m + 1L]
      num <- 1L + count_leq_neg_t(t_m)
      den <- count_geq_pos_t(t_m)
      ratio <- if (den == 0L) Inf else num / den
      if (ratio > gamma && m <= p) {
        m <- m + 1L
      } else {
        break
      }
    }
    tau <- t_all[m]
    selected_in_old <- which(res >= tau)
    selected <- top_idx[selected_in_old]
    list(gamma = gamma, tau = tau, m = m, selected_idx = selected)
  }
  
  # 5) For each gamma, output selection results
  results_list <- lapply(gammas, solve_for_gamma)
  
  
  list(
    details = results_list,       # selected_idx for each gamma
    res = res,                    # difference scores, convenient for checking/plotting
    top_idx = top_idx             # original indices of the top 50 features
  )
}


############## compared methods -- KCQD ################
u4_stat_groups_parallel <- function(X, Y, tau, beta, groups = NULL,
                                   n_workers = max(1, parallel::detectCores() - 5)) {
  n <- length(Y)
  
  ## --------- Internal helper functions ----------
  # d(Y,beta): 1D case -> sign matrix with diag=0
  d_matrix_1d <- function(y, beta, eps = 1e-12) {
    u <- beta - y
    u_norm <- u / pmax(abs(u), eps)
    D <- outer(u_norm, u_norm, "*")
    diag(D) <- 0
    D
  }
  # Gaussian kernel C for vectors: rows are samples i=1..n, columns are features in this group
  c_matrix_vec <- function(Xg) {
    # Pairwise squared distances between row-vectors
    G <- tcrossprod(Xg)                       # n x n
    r2 <- rowSums(Xg^2)
    sq <- outer(r2, r2, "+") - 2 * G         # ||xi-xj||^2
    # Median heuristic for sigma (take the square root)
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
  
  ## --------- Precomputation (shared by all groups) ----------
  D <- d_matrix_1d(Y, beta)
  comb4 <- combn(n, 4)
  
  ## --------- Group definition ----------
  p <- ncol(X)
  if (is.null(groups)) {
    # By default, split columns sequentially into 50 groups (20 columns per group)
    stopifnot(p %% 50 == 0)
    gsize <- p / 50
    groups <- split(seq_len(p), rep(1:50, each = gsize))
  } else {
    # User-defined: list where each element is the column indices of one group
    stopifnot(is.list(groups))
  }
  
  ## --------- Parallel computation for each group ----------
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  
  vals <- foreach::foreach(g = seq_along(groups), .combine = c,
                           .packages = character()) %dopar% {
                             Xg <- X[, groups[[g]], drop = FALSE]   # n x d (here d=20)
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
  
  # For each row group, run once and obtain u4 results for its 50 column groups
  res_mat <- matrix(NA_real_, nrow = n_row_groups, ncol = n_col_groups)
  for (g in seq_along(row_groups)) {
    idx <- row_groups[[g]]
    Xg  <- X[idx, , drop = FALSE]
    Yg  <- Y[idx]
    
    # Inner parallelization is across column groups (avoid excessive nested parallelism)
    out <- u4_stat_groups_parallel(
      Xg, Yg, tau, beta,
      groups = col_groups,          # explicitly pass 50 column groups
      n_workers = n_workers_col     # number of workers for column-group parallelism
    )
    # out is numeric values for each group
    res_mat[g, ] <- out
  }
  
  # Take the column means over 20 row groups to obtain a 50×1 vector
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


KCQD <- function(Y, X, alpha, active_index, gammas) {
  # 1) Split the data
  X1 <- X[1:200, , drop = FALSE]
  Y1 <- Y[1:200]
  X2 <- X[201:400, , drop = FALSE]
  Y2 <- Y[201:400]
  gsize = ncol(X)/50
  
  # 2) First stage: compute group statistics and store (T_vals length = number of groups)
  T_vals <- u4_stat_rowgroups_mean(X1, Y1, alpha)
  
  # 3) Solver for a single gamma
  solve_for_gamma <- function(g) {
    # 3.1 Select groups (original group indices)
    group_idx <- as.integer(select_groups_by_t(T_vals, 0, g))
    group_idx <- unique(group_idx)
    group_idx <- group_idx[!is.na(group_idx)]
    
    if (length(group_idx) == 0L) {
      group_idx = 1
    }
    
    # 3.2 Expand selected groups to their original column indices
    element_position <- unlist(lapply(group_idx, function(i) ((i - 1L) * gsize + 1L):(i * gsize)))
    element_position <- as.integer(element_position)
    # Truncate to the column range of X (for safety)
    element_position <- element_position[element_position >= 1L & element_position <= ncol(X2)]
    
    # 3.3 Second stage: compute statistics on candidate groups (note n_col_groups is number of groups, not columns)
    n_groups2 <- length(element_position)
    s_vec <- u4_stat_rowgroups_mean(X2[, element_position, drop = FALSE], Y2, alpha,
                                    n_col_groups = n_groups2)
    
    # The length of s_vec should match the number of selected groups
    if (length(s_vec) != n_groups2) {
      stop("Length mismatch: s_vec should have the same length as the number of selected groups; please check the return value of u4_stat_rowgroups_mean.")
    }
    
    # 3.4 Adaptive thresholding
    p <- length(s_vec)
    absD     <- abs(s_vec)
    t_sorted <- sort(absD, decreasing = FALSE)
    t_all    <- c(0, t_sorted, Inf)   # t_0=0, t_{p+1}=Inf
    
    count_leq_neg_t <- function(t) sum(s_vec <= -t)
    count_geq_pos_t <- function(t) sum(s_vec >=  t)
    
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
    
    # 3.5 Selection and mapping
    selected_in_svec <- which(s_vec >= tau)         # indices relative to s_vec
    selected_features <- element_position[selected_in_svec] 
    
    list(
      gamma = g, tau = tau, m = m,
      selected_idx = selected_in_svec,     # indices within s_vec
      selected_features = selected_features # original column indices in X
    )
  }
  
  # 4) Compute for each gamma
  details_list <- lapply(gammas, solve_for_gamma)
  names(details_list) <- paste0("gamma=", format(gammas, trim = TRUE))
  
  
  # 6) Also return selected groups/features for each gamma as lists for downstream use
  selected_features_list <- lapply(details_list, function(z) z$selected_features)
  names(selected_features_list) <- names(details_list)
  
  list(
    details = details_list,
    selected_features = selected_features_list  # original column indices for each gamma
  )
}


############## simulation 2b ##############
generate_model_2c <- function(n = 400, p = 2000, n_active) {
  # Set covariance matrix Σ_{ij} = 0.8^{|i-j|}
  Sigma <- outer(1:p, 1:p, function(i, j) 0.8^abs(i - j))
  
  # Generate X ~ N(0, Σ)
  library(MASS)
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # Generate error term ε ~ N(0, 1)
  epsilon <- rnorm(n, mean = 0, sd = 1)
  
  # Construct Y according to Model 2.c
  Y <- 5 * rowSums(log(0.5 * (abs(X[, 1:2]) / abs(1 - X[, 1:2])))) +
    5 * rowSums(sin((pi / 2) * X[, 3:5])) +
    5 * rowSums(exp(0.5 * X[, 6:8])) +
    5 * rowSums(X[, 9:n_active ]) +
    epsilon
  
  return(list(X = X, Y = Y))
}


n <- 400
p <- 2000
alpha <- 0.75
n_active <- 12
nrep_total <- 100
gammas_vec <- c(0.1, 0.15)

hh <- alpha_constants(alpha)
h  <- (hh$M2/(4*hh$M4))^(1/5)*n^(-1/5)
delta_n_alpha <- if (alpha == 0.5) 2.858*n^(-0.1) else 2.455*n^(-0.1)

G <- length(gammas_vec)

## -------- Split 100 repetitions into 5 blocks of 20 --------
n_blocks <- 5
nrep_block <- nrep_total / n_blocks   # 20
stopifnot(nrep_block == floor(nrep_block))


## -------- Function to run one block: return all results for that block --------
run_one_block <- function(block_id, seed) {
  set.seed(seed)
  
  # Initialize containers for this block
  sel_sizes_pc    <- matrix(0L,     nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  fdr_vals_pc     <- matrix(0.0,    nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  all_selected_pc <- matrix(FALSE,  nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  sel_matrix_pc   <- array(FALSE, dim = c(nrep_block, p, G), dimnames = list(NULL, NULL, paste0("g=", gammas_vec)))
  
  sel_sizes_fs    <- matrix(0L,     nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  fdr_vals_fs     <- matrix(0.0,    nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  all_selected_fs <- matrix(FALSE,  nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  sel_matrix_fs   <- array(FALSE, dim = c(nrep_block, p, G), dimnames = list(NULL, NULL, paste0("g=", gammas_vec)))
  
  sel_sizes_KCQD    <- matrix(0L,     nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  fdr_vals_KCQD     <- matrix(0.0,    nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  all_selected_KCQD <- matrix(FALSE,  nrow = nrep_block, ncol = G, dimnames = list(NULL, paste0("g=", gammas_vec)))
  sel_matrix_KCQD   <- array(FALSE, dim = c(nrep_block, p, G), dimnames = list(NULL, NULL, paste0("g=", gammas_vec)))
  
  time_fs   <- numeric(nrep_block)
  time_pc   <- numeric(nrep_block)
  time_kcqd <- numeric(nrep_block)
  
  for (r in 1:nrep_block) {
    cat(sprintf("Block %d | Rep %d / %d\n", block_id, r, nrep_block))
    
    ## Step 1: Generate data
    data <- generate_model_2c(n = n, p = p, n_active )
    X <- data$X
    Y <- data$Y
    
    ## Step 2: FS
    t0 <- proc.time()
    comp <- compute_D_dagger_from_alpha(X, Y, h = h, alpha = alpha, delta_n_alpha = delta_n_alpha)
    res  <- adaptive_threshold_select(comp$D_dagger, gammas = gammas_vec)
    time_fs[r] <- (proc.time() - t0)[["elapsed"]]
    
    for (g in seq_len(G)) {
      S_hat <- res$details[[g]]$selected_idx
      sel_sizes_fs[r, g] <- length(S_hat)
      false_discoveries  <- sum(S_hat > n_active)
      fdr_vals_fs[r, g]  <- false_discoveries / max(1, length(S_hat))
      all_selected_fs[r, g] <- all(1:n_active %in% S_hat)
      if (length(S_hat) > 0) sel_matrix_fs[r, S_hat, g] <- TRUE
    }
    
    ## Step 3: PC
    t0 <- proc.time()
    res <- PCSIS(Y, X, gammas = gammas_vec, active_set = c(1:n_active))
    time_pc[r] <- (proc.time() - t0)[["elapsed"]]
    
    for (g in seq_len(G)) {
      S_hat <- res$details[[g]]$selected_idx
      sel_sizes_pc[r, g] <- length(S_hat)
      false_discoveries  <- sum(S_hat > n_active)
      fdr_vals_pc[r, g]  <- false_discoveries / max(1, length(S_hat))
      all_selected_pc[r, g] <- all(1:n_active %in% S_hat)
      if (length(S_hat) > 0) sel_matrix_pc[r, S_hat, g] <- TRUE
    }
    
    ## Step 4: KCQD
    t0 <- proc.time()
    res <- KCQD(Y, X, alpha, 1:n_active , gammas_vec)
    time_kcqd[r] <- (proc.time() - t0)[["elapsed"]]
    
    for (g in seq_len(G)) {
      S_hat <- res$details[[g]]$selected_idx
      sel_sizes_KCQD[r, g] <- length(S_hat)
      false_discoveries    <- sum(S_hat > n_active)
      fdr_vals_KCQD[r, g]  <- false_discoveries / max(1, length(S_hat))
      all_selected_KCQD[r, g] <- all(1:n_active %in% S_hat)
      if (length(S_hat) > 0) sel_matrix_KCQD[r, S_hat, g] <- TRUE
    }
  }
  
  list(
    block_id = block_id, seed = seed,
    sel_sizes_pc = sel_sizes_pc, fdr_vals_pc = fdr_vals_pc, all_selected_pc = all_selected_pc, sel_matrix_pc = sel_matrix_pc,
    sel_sizes_fs = sel_sizes_fs, fdr_vals_fs = fdr_vals_fs, all_selected_fs = all_selected_fs, sel_matrix_fs = sel_matrix_fs,
    sel_sizes_KCQD = sel_sizes_KCQD, fdr_vals_KCQD = fdr_vals_KCQD, all_selected_KCQD = all_selected_KCQD, sel_matrix_KCQD = sel_matrix_KCQD,
    time_fs = time_fs, time_pc = time_pc, time_kcqd = time_kcqd
  )
}

## -------- Run each block and save as .rds --------
blk <- run_one_block(block_id = 1, seed = 1234)
saveRDS(blk, file = sprintf("sim_block_%02d.rds", 1234))
blk <- run_one_block(block_id = 2, seed = 1235)
saveRDS(blk, file = sprintf("sim_block_%02d.rds", 1235))
blk <- run_one_block(block_id = 3, seed = 1236)
saveRDS(blk, file = sprintf("sim_block_%02d.rds", 1236))
blk <- run_one_block(block_id = 4, seed = 1237)
saveRDS(blk, file = sprintf("sim_block_%02d.rds", 1237))
blk <- run_one_block(block_id = 5, seed = 1238)
saveRDS(blk, file = sprintf("sim_block_%02d.rds", 1238))

## -------- Read and combine all blocks --------
files <- sprintf("sim_block_%02d.rds", c(1234,1235,1236,1237,1238))
blocks <- lapply(files, readRDS)

## Helper: row-bind matrix-like objects
bind_mtx <- function(field) do.call(rbind, lapply(blocks, `[[`, field))

## Helper: concatenate 3D arrays along the first dimension (rep dimension)
bind_arr <- function(field) {
  arrs <- lapply(blocks, `[[`, field)
  stopifnot(length(unique(lapply(arrs, function(a) dim(a)[2]))) == 1)  # consistent p
  stopifnot(length(unique(lapply(arrs, function(a) dim(a)[3]))) == 1)  # consistent G
  p_local <- dim(arrs[[1]])[2]
  G_local <- dim(arrs[[1]])[3]
  n1s <- sapply(arrs, function(a) dim(a)[1])
  out <- array(FALSE, dim = c(sum(n1s), p_local, G_local),
               dimnames = list(NULL, NULL, dimnames(arrs[[1]])[[3]]))
  idx <- 1
  for (a in arrs) {
    nblk <- dim(a)[1]
    out[idx:(idx + nblk - 1), , ] <- a
    idx <- idx + nblk
  }
  out
}

## Combine to obtain the final results over 100 repetitions
sel_sizes_pc_all     <- bind_mtx("sel_sizes_pc")
fdr_vals_pc_all      <- bind_mtx("fdr_vals_pc")
all_selected_pc_all  <- bind_mtx("all_selected_pc")
sel_matrix_pc_all    <- bind_arr("sel_matrix_pc")

sel_sizes_fs_all     <- bind_mtx("sel_sizes_fs")
fdr_vals_fs_all      <- bind_mtx("fdr_vals_fs")
all_selected_fs_all  <- bind_mtx("all_selected_fs")
sel_matrix_fs_all    <- bind_arr("sel_matrix_fs")

sel_sizes_KCQD_all    <- bind_mtx("sel_sizes_KCQD")
fdr_vals_KCQD_all     <- bind_mtx("fdr_vals_KCQD")
all_selected_KCQD_all <- bind_mtx("all_selected_KCQD")
sel_matrix_KCQD_all   <- bind_arr("sel_matrix_KCQD")

time_fs_all   <- unlist(lapply(blocks, `[[`, "time_fs"))
time_pc_all   <- unlist(lapply(blocks, `[[`, "time_pc"))
time_kcqd_all <- unlist(lapply(blocks, `[[`, "time_kcqd"))

## ----- Step 3: Summarize -----
# FS
median_size_fs <- apply(sel_sizes_fs_all, 2, median)
Pj_fs <- sapply(1:G, function(g) {
  colMeans(sel_matrix_fs_all[, 1:n_active, g, drop = FALSE])
})
prob_all_fs <- colMeans(all_selected_fs_all)
emp_fdr_fs  <- colMeans(fdr_vals_fs_all)

names(median_size_fs) <- paste0("g=", gammas_vec)
colnames(Pj_fs)       <- paste0("g=", gammas_vec)
names(prob_all_fs)    <- paste0("g=", gammas_vec)
names(emp_fdr_fs)     <- paste0("g=", gammas_vec)

# PC
median_size_pc <- apply(sel_sizes_pc_all, 2, median)
Pj_pc <- sapply(1:G, function(g) {
  colMeans(sel_matrix_pc_all[, 1:n_active, g, drop = FALSE])
})
prob_all_pc <- colMeans(all_selected_pc_all)
emp_fdr_pc  <- colMeans(fdr_vals_pc_all)

names(median_size_pc) <- paste0("g=", gammas_vec)
colnames(Pj_pc)       <- paste0("g=", gammas_vec)
names(prob_all_pc)    <- paste0("g=", gammas_vec)
names(emp_fdr_pc)     <- paste0("g=", gammas_vec)

# KCQD
median_size_kcqd <- apply(sel_sizes_KCQD_all, 2, median)
Pj_kcqd <- sapply(1:G, function(g) {
  colMeans(sel_matrix_KCQD_all[, 1:n_active, g, drop = FALSE])
})
prob_all_kcqd <- colMeans(all_selected_KCQD_all)
emp_fdr_kcqd  <- colMeans(fdr_vals_KCQD_all)

names(median_size_kcqd) <- paste0("g=", gammas_vec)
colnames(Pj_kcqd)       <- paste0("g=", gammas_vec)
names(prob_all_kcqd)    <- paste0("g=", gammas_vec)
names(emp_fdr_kcqd)     <- paste0("g=", gammas_vec)

## ----- Save results -----
summary_results <- list(
  FS   = list(median_size = median_size_fs, Pj = Pj_fs,
              prob_all = prob_all_fs, emp_fdr = emp_fdr_fs),
  PC   = list(median_size = median_size_pc, Pj = Pj_pc,
              prob_all = prob_all_pc, emp_fdr = emp_fdr_pc),
  KCQD = list(median_size = median_size_kcqd, Pj = Pj_kcqd,
              prob_all = prob_all_kcqd, emp_fdr = emp_fdr_kcqd)
)

saveRDS(summary_results, file = "model_3c_alpha_075.rds")

# Extract summary_results and tidy into a table
library(dplyr)
library(tidyr)
library(purrr)

df <- map_dfr(names(summary_results), function(method) {
  res <- summary_results[[method]]
  tibble(
    Method      = method,
    Gamma       = names(res$median_size),
    Median_Size = as.numeric(res$median_size),
    Prob_All    = as.numeric(res$prob_all),
    Emp_FDR     = as.numeric(res$emp_fdr)
  )
})

print(df)

method_names <- names(summary_results)

df_pj_wide <- purrr::map_dfr(method_names, function(method) {
  Pj <- summary_results[[method]]$Pj
  Pj <- as.matrix(Pj)                 # ensure it is a matrix
  gammas <- colnames(Pj)
  K <- nrow(Pj)
  feature_cols <- paste0("F", seq_len(K))
  
  purrr::map_dfr(seq_along(gammas), function(j) {
    vals <- as.numeric(Pj[, j])       
    tibble(
      Method = method,
      Gamma  = gammas[j],
      !!!setNames(as.list(round(vals, 2)), feature_cols)
    )
  })
})

df_pj_wide




