# 需要的包列表
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

# 检查并安装缺失的包
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# 加载所有包
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

# (2.8) Copula Quantile Dependence  —— 向量化版本
estimate_D_star <- function(s_vec, t_vec, x_vec, h = 0.1, alpha = 0.5, omega = function(u) 1) {
  n <- length(s_vec)
  F_Xj <- ecdf(x_vec)
  u_hat <- (n / (n + 1)) * F_Xj(x_vec)
  
  # 对所有 u_hat 一次性计算 C(u_i, alpha | s_vec, t_vec)
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

# 主函数：外层并行
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
  
  # 并行设置
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
  
  # centered Y at quantile tau (和原来的 alpha 逻辑一致)
  Y <- Y - as.numeric(quantile(Y, tau))
  
  # 标准化 X（列均值0方差1）
  X <- scale(X)
  
  # 并行设置
  cl <- parallel::makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  
  # 并行计算每个 j 的 score
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
  rx <- colRanks(as.matrix(X))  # 去掉了 parallel=TRUE
  gry <- gt(ry)
  grx <- apply(rx, 2, gt)
  
  X <- as.matrix(grx)
  Y <- gry
  
  G <- Rfast::rmvnorm(n, mu = rep(0, 2), sigma = diag(2))
  
  # 并行设置
  cl <- parallel::makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  
  # 并行计算
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
  
  # 1) 预处理：居中Y、标准化X（对0方差列做简易保护）
  Y <- Y - as.numeric(stats::quantile(Y, tau, names = FALSE))
  mu  <- colMeans(X)
  sdv <- apply(X, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
  X <- sweep(sweep(X, 2, mu, `-`), 2, sdv, `/`)
  
  # 2) 并行：只创建一次集群
  cl <- parallel::makeCluster(max(1, n_workers))
  doParallel::registerDoParallel(cl)
  on.exit({
    foreach::registerDoSEQ()
    parallel::stopCluster(cl)
  }, add = TRUE)
  
  # 3) 第一次打分：样条 + 分位回归（失败就记0分，别让整体报错）
  s_vec <- foreach::foreach(
    j = 1:p, .combine = "c",
    .packages = c("splines", "quantreg")
  ) %dopar% {
    x0 <- X[, j]
    # 结点按 1/3、2/3 分位，退化就少用/不用
    qs <- unique(stats::quantile(x0, c(1/3, 2/3), names = FALSE, type = 7))
    a  <- tryCatch(splines::bs(x0, knots = qs, degree = 1), error = function(e) NULL)
    if (is.null(a)) return(0)
    fit <- tryCatch(quantreg::rq(Y ~ a, tau = tau), error = function(e) NULL)
    if (is.null(fit)) return(0)
    sum(fitted(fit)^2) / n
  }
  
  # 4) 阈值与候选集合
  if (is.null(active_index) || length(active_index) == 0) {
    threshold <- -Inf
  } else {
    ai <- intersect(unique(as.integer(active_index)), seq_len(p))
    threshold <- suppressWarnings(min(s_vec[ai], na.rm = TRUE))
    if (!is.finite(threshold)) threshold <- -Inf
  }
  idx <- which(s_vec >= threshold)
  if (length(idx) == 0L) return(0L)
  
  # 5) 第二次打分：pcov（失败记0）
  res <- foreach::foreach(
    i = idx, .combine = "c",
    .export = "pcov"
  ) %dopar% {
    val <- tryCatch(pcov(X[, i], Y), error = function(e) 0)
    if (!is.finite(val)) 0 else val
  }
  
  ranked_index   <- order(res, decreasing = TRUE)
  original_index <- idx[ranked_index]
  
  # 6) 最小模型大小
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
  
  ## --------- 内部函数 ----------
  # d(Y,beta): 1D 情况 -> sign 矩阵，diag=0
  d_matrix_1d <- function(y, beta, eps = 1e-12) {
    u <- beta - y
    u_norm <- u / pmax(abs(u), eps)
    D <- outer(u_norm, u_norm, "*")
    diag(D) <- 0
    D
  }
  # 向量高斯核的 C：行是样本 i=1..n，列是该组的特征
  c_matrix_vec <- function(Xg) {
    # pairwise squared distances between row-vectors
    G <- tcrossprod(Xg)                       # n x n
    r2 <- rowSums(Xg^2)
    sq <- outer(r2, r2, "+") - 2 * G         # ||xi-xj||^2
    # median heuristic for sigma (开根号！)
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
  
  ## --------- 预计算（所有组共享） ----------
  D <- d_matrix_1d(Y, beta)
  comb4 <- combn(n, 4)
  
  ## --------- 组定义 ----------
  p <- ncol(X)
  if (is.null(groups)) {
    # 默认按列顺序均分为 50 组（每组 20 列）
    stopifnot(p %% 50 == 0)
    gsize <- p / 50
    groups <- split(seq_len(p), rep(1:50, each = gsize))
  } else {
    # 用户自定义：list，每个元素是一组的列下标
    stopifnot(is.list(groups))
  }
  
  ## --------- 并行计算每一组 ----------
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  
  vals <- foreach::foreach(g = seq_along(groups), .combine = c,
                           .packages = character()) %dopar% {
                             Xg <- X[, groups[[g]], drop = FALSE]   # n x d (这里 d=20)
                             C <- c_matrix_vec(Xg)
                             u4_from_C_D_mean(C, D, comb4)
                           }
  as.numeric(vals)  # 长度 = 组数（默认 50）
}


u4_stat_rowgroups_mean <- function(X, Y, tau,
                                   n_row_groups = 20,
                                   n_col_groups = 50,
                                   n_workers_col = max(1, parallel::detectCores() - 5)) {
  n <- nrow(X); p <- ncol(X)
  beta <- as.numeric(quantile(Y, (tau + 1) / 2, type = 7))
  
  # 生成行分组与列分组
  rows_per_group <- n / n_row_groups
  row_groups <- split(seq_len(n), rep(seq_len(n_row_groups), each = rows_per_group))
  
  cols_per_group <- p / n_col_groups
  col_groups <- split(seq_len(p), rep(seq_len(n_col_groups), each = cols_per_group))
  
  # 对每个“行组”跑一遍，取出该行组下 50 个列组的 u4 结果
  res_mat <- matrix(NA_real_, nrow = n_row_groups, ncol = n_col_groups)
  for (g in seq_along(row_groups)) {
    idx <- row_groups[[g]]
    Xg  <- X[idx, , drop = FALSE]
    Yg  <- Y[idx]
    
    # 内层并行还是列组并行（避免双层并行过度订阅）
    out <- u4_stat_groups_parallel(
      Xg, Yg, tau, beta,
      groups = col_groups,          # 显式传入 50 个列组
      n_workers = n_workers_col     # 列组并行的 worker 数
    )
    # out 是 data.frame: group, u4, sigma
    res_mat[g, ] <- out
  }
  
  # 对 20 个行组在列方向取平均，得到 50×1 的向量
  colMeans(res_mat)
}

select_groups_by_t <- function(T_vals, u, alpha_G, include_equal = TRUE) {
  T_vals <- as.numeric(T_vals)
  
  # 计算 t_ne
  t_candidates <- sort(unique(abs(T_vals)))
  t_ne <- NA_real_
  for (t in t_candidates) {
    num <- u + sum(T_vals <= -t)
    den <- sum(T_vals >=  t)
    ratio <- if (den > 0) num / den else Inf
    if (ratio <= alpha_G) { t_ne <- t; break }
  }
  
  # 根据 t_ne 选组
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

## ==== config ====
n <- 400; p <- 5000
rep_total <- 100
chunks <- 5
rep_per_chunk <- rep_total / chunks  # 20
seeds <- c(1234, 1235, 1236, 1237, 1238)  # 每个 chunk 一个 seed
alpha <- 0.5                               # 如需 0.5 就改成 0.5
active.set <- c(1, 2)                       # 你给定的 active.set

## ==== 生成器（保持你的写法）====
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

## ==== 容器 ====
init_containers_1c <- function() {
  list(
    # 结果
    fs_result   = c(),
    sis_result  = c(),
    dc_sis_result = c(),
    qa_result   = c(),
    qc_result   = c(),
    wd_result   = c(),
    pc_result   = c(),
    KCQD_result = c(),
    # 时间
    fs_time   = c(),
    sis_time  = c(),
    dc_sis_time = c(),
    qa_time   = c(),
    qc_time   = c(),
    wd_time   = c(),
    pc_time   = c(),
    KCQD_time = c()
  )
}

## ==== 跑一个 chunk（20 次）并保存为 .RData ====
run_chunk_1c <- function(chunk_id, seed, n, p, rep_per_chunk, alpha, active.set) {
  set.seed(seed)
  C <- init_containers_1c()
  
  for (i in seq_len(rep_per_chunk)) {
    dat <- generate_model_1c(n, p, alpha = alpha)
    X <- dat$X; Y <- dat$Y
    
    # FS
    t0 <- proc.time()
    res_fs <- qc_screen_with_min_model_size(X, Y, alpha, active_index = active.set)
    C$fs_time   <- c(C$fs_time,   (proc.time() - t0)[["elapsed"]])
    C$fs_result <- c(C$fs_result, res_fs$min_model_size)
    
    # SIS
    t0 <- proc.time()
    C$sis_result <- c(C$sis_result, max(Rfast::Rank(abs(cor(X, Y, method = "pearson")), descending = TRUE)[active.set]))
    C$sis_time   <- c(C$sis_time,   (proc.time() - t0)[["elapsed"]])
    
    # DC-SIS (ensure mdcor is available; otherwise replace with column-wise dcor)
    t0 <- proc.time()
    C$dc_sis_result <- c(C$dc_sis_result, max(Rfast::Rank(mdcor(Y, X, "V"), descending = TRUE)[active.set]))
    C$dc_sis_time   <- c(C$dc_sis_time,   (proc.time() - t0)[["elapsed"]])
    
    # QaSIS
    t0 <- proc.time()
    C$qa_result <- c(C$qa_result, QaSIS(Y, X, tau = alpha, active.set))
    C$qa_time   <- c(C$qa_time,   (proc.time() - t0)[["elapsed"]])
    
    # QcSIS
    t0 <- proc.time()
    C$qc_result <- c(C$qc_result, QcSIS(Y, X, alpha, active.set))
    C$qc_time   <- c(C$qc_time,   (proc.time() - t0)[["elapsed"]])
    
    # WD
    t0 <- proc.time()
    C$wd_result <- c(C$wd_result, max(Rfast::Rank(wdGscreening(Y, X, ppp = 1), descending = TRUE)[active.set]))
    C$wd_time   <- c(C$wd_time,   (proc.time() - t0)[["elapsed"]])
    
    # PC
    #t0 <- proc.time()
    #C$pc_result <- c(C$pc_result, PCSIS(Y, X, active.set))
    #C$pc_time   <- c(C$pc_time,   (proc.time() - t0)[["elapsed"]])
    
    # KCQD (fix original typo: avoid KCQD_result <- c(KCQD_result, ...))
    t0 <- proc.time()
    C$KCQD_result <- c(C$KCQD_result, KCQD(Y, X, alpha, active.set))
    C$KCQD_time   <- c(C$KCQD_time,   (proc.time() - t0)[["elapsed"]])
    
    cat(sprintf("model1c | chunk %d | iter %d/%d done\n", chunk_id, i, rep_per_chunk))
  }
  
  save_path <- file.path(sprintf("model1c_chunk_%d_seed_%s_alpha_%s.RData", chunk_id, seed, gsub("\\.", "_", as.character(alpha))))
  save(C, n, p, rep_per_chunk, alpha, active.set, seed, file = save_path)
  invisible(save_path)
}

## ==== 开跑（可分 5 次不同会话分别跑）====
run_chunk_1c(chunk_id = 1, seed = 1234, n = n, p = p,
             rep_per_chunk = rep_per_chunk, alpha = alpha, active.set = active.set)
run_chunk_1c(chunk_id = 2, seed = 1235, n = n, p = p,
             rep_per_chunk = rep_per_chunk, alpha = alpha, active.set = active.set)
run_chunk_1c(chunk_id = 3, seed = 1236, n = n, p = p,
             rep_per_chunk = rep_per_chunk, alpha = alpha, active.set = active.set)
run_chunk_1c(chunk_id = 4, seed = 1237, n = n, p = p,
             rep_per_chunk = rep_per_chunk, alpha = alpha, active.set = active.set)
run_chunk_1c(chunk_id = 5, seed = 1238, n = n, p = p,
             rep_per_chunk = rep_per_chunk, alpha = alpha, active.set = active.set)


## ==== 工具 ====
fmt_median_iqr <- function(x) sprintf("%.2f (%.2f)", stats::median(x, na.rm = TRUE), IQR(x, na.rm = TRUE))

summarize_named_vecs <- function(lst_named_vecs) {
  sapply(lst_named_vecs, function(x) fmt_median_iqr(x))
}

merge_flat_lists <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  keys <- union(names(a), names(b))
  out <- vector("list", length(keys)); names(out) <- keys
  for (k in keys) {
    xa <- a[[k]]; xb <- b[[k]]
    if (is.null(xa)) out[[k]] <- xb
    else if (is.null(xb)) out[[k]] <- xa
    else out[[k]] <- c(xa, xb)
  }
  out
}

## ==== Read and merge all chunk results ====
files <- list.files(pattern = "^model1c_chunk_\\d+_seed_.*_alpha_.*\\.RData$", full.names = TRUE)
stopifnot(length(files) > 0)

C_all <- NULL
for (f in files) {
  env <- new.env()
  load(f, envir = env)   # contains C, alpha, active.set
  C_all <- merge_flat_lists(C_all, env$C)
}
# If you later run multiple alpha values and want separate summaries, you can use alpha as a key

## ==== Assemble and summarize ====
vars <- list(
  fs_result   = C_all$fs_result,
  sis_result  = C_all$sis_result,
  dc_sis_result = C_all$dc_sis_result,
  qa_result   = C_all$qa_result,
  qc_result   = C_all$qc_result,
  KCQD_result = C_all$KCQD_result,
  pc_result   = C_all$pc_result,
  wd_result   = C_all$wd_result
)

times <- list(
  fs_time   = C_all$fs_time,
  sis_time  = C_all$sis_time,
  dc_sis_time = C_all$dc_sis_time,
  qa_time   = C_all$qa_time,
  qc_time   = C_all$qc_time,
  KCQD_time = C_all$KCQD_time,
  pc_time   = C_all$pc_time,
  wd_time   = C_all$wd_time
)

result_table <- summarize_named_vecs(vars)
time_table   <- summarize_named_vecs(times)

result_table <- vapply(result_table,
               function(x) if (length(x) == 0) NA_character_ else as.character(x),
               character(1))

time_table <- vapply(time_table,
                       function(x) if (length(x) == 0) NA_character_ else as.character(x),
                       character(1))

result_df <- data.frame(Metric = names(result_table), Value = unname(result_table), row.names = NULL)
time_df   <- data.frame(Metric = names(time_table),   Value = unname(time_table),   row.names = NULL)

## ==== Save results ====
alpha_str <- gsub("\\.", "_", as.character(alpha))
save_path_rdata <- sprintf("model1c_merged_n%d_p%d_alpha%s.RData", n, p, alpha_str)
save(vars, times, result_table, time_table, file = save_path_rdata)

write.csv(result_df, sprintf("model1c_results_n%d_p%d_alpha%s.csv", n, p, alpha_str), row.names = FALSE)
write.csv(time_df,   sprintf("model1c_time_n%d_p%d_alpha%s.csv",    n, p, alpha_str), row.names = FALSE)

## ==== Quick view + length check (expected 100) ====
print(result_table)
print(time_table)
cat("lens:", length(vars$fs_result), length(vars$sis_result), length(vars$dc_sis_result), "\n")
