# Required packages
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

# φ'(Φ^{-1}(α)) and four constants/variances
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

# The existing estimate_D_star() will be directly called
compute_D_dagger_from_alpha <- function(X, Y, h, alpha, omega=function(u) dnorm(qnorm(u))^2,
                                        delta_n_alpha = 0, eps = 1e-12, n_workers = max(1, parallel::detectCores() - 5)) {
  stopifnot(is.matrix(X), length(Y) == nrow(X))
  n <- length(Y)
  
  # 1) First compute \hat D_α
  F_Y <- ecdf(Y)
  t_vec <- qnorm((n / (n + 1)) * F_Y(Y))  # Transformed Y
  
  F_Xj <- ecdf(X[, 1])
  s_vec <- qnorm((n / (n + 1)) * F_Xj(X[, 1]))
  D_alpha_vals <- estimate_D_star(s_vec, t_vec, X[, 1], h = h, alpha = alpha, omega = omega)

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


############## simulation ##############
# Define sequences of alpha and n
alphas <- seq(0.1, 0.9, by = 0.05)
ns <- c(100, 200, 400, 800, 1600, 3200)

# Store results
results <- data.frame(alpha = numeric(),
                      n = integer(),
                      mean_res1 = numeric())

# Main loop
for (alpha in alphas) {
  hh <- alpha_constants(alpha)
  for (n in ns) {
    h <- (hh$M2 / (4 * hh$M4))^(1/5) * n^(-1/5)
    
    res1 <- numeric(200)
    for (rep in 1:200) {
      X <- as.matrix(rnorm(n), nrow = n)
      Y <- as.matrix(rnorm(n), nrow = n)
      
      comp <- compute_D_dagger_from_alpha(X, Y, h = h, alpha = alpha, delta_n_alpha = 0)
      res1[rep] <- -comp$D_dagger
    }
    
    # Save the mean
    results <- rbind(results,
                     data.frame(alpha = alpha,
                                n = n,
                                mean_res1 = mean(res1)))
    cat("Done: alpha =", alpha, ", n =", n, "\n")
  }
}

# Inspect results
print(results)

# Optional: save to CSV
write.csv(results, "mean_res1_results.csv", row.names = FALSE)

results <- read.csv("mean_res1_results.csv")


library(dplyr)
library(ggplot2)

# 1) Preprocessing: take logarithm (avoid log(0) with a very small protection)
plot_df <- results %>%
  mutate(
    alpha_f = factor(sprintf("%.2f", alpha)),             
    logn    = log(n),
    logm    = log(pmax(mean_res1, .Machine$double.eps))  
  )

# 2) Plot: scatter + linear regression for each α + confidence band
p <- ggplot(plot_df, aes(x = logn, y = logm, color = alpha_f, fill = alpha_f)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  labs(
    x = expression(log(n)),
    y = expression(log(delta[n]^"*"~"("*alpha*")")),
    color = expression(alpha),
    fill  = expression(alpha)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(),
    legend.position = "right"
  )


print(p)


# Save figure
ggsave("FigureB1.png", p, width = 9, height = 6, dpi = 300)

