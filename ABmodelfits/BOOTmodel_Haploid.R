BOOTmodel_Haploid <- function(model.fit, Nboot, out.dir, out.name)
{
  allow.neg.intercept <- "no"
  
  # ---------------------------------------------------------------------------
  # 1) Read model output
  # ---------------------------------------------------------------------------
  if (is.character(model.fit) && length(model.fit) == 1) {
    fit <- dget(model.fit)
  } else if (is.list(model.fit)) {
    fit <- model.fit
  } else {
    stop("model.fit must be either a filepath (character) or a model output list.")
  }
  
  required_names <- c("estimates", "pedigree", "settings", "model")
  if (!all(required_names %in% names(fit))) {
    stop("Input object does not look like ABneutralHaploid output.")
  }
  
  if (!identical(fit$model, "ABneutralHaploid.R")) {
    stop("This bootstrap function is only for model == 'ABneutralHaploid.R'.")
  }
  
  # ---------------------------------------------------------------------------
  # 2) Extract settings and baseline fit
  # ---------------------------------------------------------------------------
  settings <- fit$settings
  est <- fit$estimates
  ped <- as.data.frame(fit$pedigree)
  
  need_cols <- c("div.obs", "div.pred", "residual")
  if (!all(need_cols %in% colnames(ped))) {
    stop("fit$pedigree must contain columns: div.obs, div.pred, residual.")
  }
  
  eqp <- as.numeric(as.character(settings[settings[, 1] == "eqp", 2]))
  eqp.weight <- as.numeric(as.character(settings[settings[, 1] == "eqp.weight", 2]))
  optim.method <- as.character(settings[settings[, 1] == "optim.method", 2])
  p0uu <- as.numeric(as.character(settings[settings[, 1] == "p0uu", 2]))
  
  if (length(eqp) == 0 || length(eqp.weight) == 0 || length(optim.method) == 0 || length(p0uu) == 0) {
    stop("Could not parse eqp/eqp.weight/optim.method/p0uu from settings.")
  }
  
  base <- est[1, , drop = FALSE]
  param_names <- c("alpha", "beta", "intercept")
  if (!all(param_names %in% colnames(base))) {
    stop("fit$estimates does not contain alpha, beta, intercept.")
  }
  
  # ---------------------------------------------------------------------------
  # 3) Internal model helpers (2-state haploid)
  # ---------------------------------------------------------------------------
  make_genmatrix <- function(alpha, beta) {
    # state order: (UU, MM)
    matrix(c(
      1 - alpha, alpha,
      beta,      1 - beta
    ), nrow = 2, byrow = TRUE)
  }
  
  pair_divergence <- function(p1, p2) {
    # probability two lineages differ
    p1[1] * p2[2] + p1[2] * p2[1]
  }
  
  calc_Dt1t2 <- function(pedigree_mat, p0uu, alpha, beta) {
    p0mm <- 1 - p0uu
    svGzero <- c(p0uu, p0mm)  # (UU, MM)
    
    G <- make_genmatrix(alpha, beta)
    
    t0  <- as.numeric(pedigree_mat[, 1])
    dt1 <- as.numeric(pedigree_mat[, 2] - pedigree_mat[, 1])
    dt2 <- as.numeric(pedigree_mat[, 3] - pedigree_mat[, 1])
    
    exps <- sort(unique(c(t0, dt1, dt2)))
    pow  <- setNames(lapply(exps, function(k) expm::`%^%`(G, k)), as.character(exps))
    
    basis <- rbind(
      UU = c(1, 0),
      MM = c(0, 1)
    )
    
    Dt1t2 <- numeric(NROW(pedigree_mat))
    
    for (i in seq_len(NROW(pedigree_mat))) {
      M0 <- pow[[as.character(t0[i])]]
      M1 <- pow[[as.character(dt1[i])]]
      M2 <- pow[[as.character(dt2[i])]]
      
      svt0 <- drop(svGzero %*% M0)
      P1 <- basis %*% M1
      P2 <- basis %*% M2
      
      dUU <- pair_divergence(P1["UU", ], P2["UU", ])
      dMM <- pair_divergence(P1["MM", ], P2["MM", ])
      
      Dt1t2[i] <- svt0[1] * dUU + svt0[2] * dMM
    }
    
    puu_inf <- beta / (alpha + beta)
    list(puuinf = puu_inf, Dt1t2 = Dt1t2)
  }
  
  LSE_intercept <- function(param_int, pedigree_mat) {
    alpha     <- param_int[1]
    beta      <- param_int[2]
    intercept <- param_int[3]
    
    d <- calc_Dt1t2(pedigree_mat, p0uu = p0uu, alpha = alpha, beta = beta)
    
    sum((pedigree_mat[, 4] - intercept - d$Dt1t2)^2) +
      eqp.weight * NROW(pedigree_mat) * ((d$puuinf - eqp)^2)
  }
  
  # ---------------------------------------------------------------------------
  # 4) Bootstrap loop
  # ---------------------------------------------------------------------------
  set_start <- c(
    alpha = as.numeric(base[1, "alpha"]),
    beta = as.numeric(base[1, "beta"]),
    intercept = as.numeric(base[1, "intercept"])
  )
  
  boot_mat <- matrix(NA_real_, nrow = Nboot, ncol = 6)
  colnames(boot_mat) <- c("alpha", "beta", "intercept", "PrMMinf", "PrUUinf", "value")
  
  good <- 0L
  
  for (b in seq_len(Nboot)) {
    message("Bootstrap iteration: ", b, "/", Nboot)
    
    ped_b <- as.matrix(ped)
    ped_b[, 4] <- ped_b[, "div.pred"] + sample(ped_b[, "residual"], nrow(ped_b), replace = TRUE)
    
    opt <- try(
      suppressWarnings(
        optimx::optimx(
          par = set_start,
          fn = function(par) LSE_intercept(par, ped_b),
          method = optim.method
        )
      ),
      silent = TRUE
    )
    
    if (inherits(opt, "try-error") || NROW(opt) < 1) next
    opt <- as.data.frame(opt)[1, , drop = FALSE]
    
    # robust extraction
    if (all(c("alpha", "beta", "intercept") %in% colnames(opt))) {
      a  <- as.numeric(opt$alpha[1])
      bt <- as.numeric(opt$beta[1])
      ic <- as.numeric(opt$intercept[1])
    } else {
      num_cols <- which(vapply(opt, is.numeric, logical(1)))
      if (length(num_cols) < 3) next
      vals <- as.numeric(opt[1, num_cols[1:3]])
      a <- vals[1]; bt <- vals[2]; ic <- vals[3]
    }
    
    if ("convcode" %in% colnames(opt)) {
      if (!isTRUE(opt$convcode[1] == 0)) next
    }
    
    if (!is.finite(a) || !is.finite(bt) || !is.finite(ic)) next
    if (a <= 0 || bt <= 0) next
    if (allow.neg.intercept == "no" && ic <= 0) next
    
    PrMMinf <- a  / (a + bt)
    PrUUinf <- bt / (a + bt)
    val <- if ("value" %in% colnames(opt)) as.numeric(opt$value[1]) else NA_real_
    
    good <- good + 1L
    boot_mat[good, ] <- c(a, bt, ic, PrMMinf, PrUUinf, val)
  }
  
  boot_mat <- boot_mat[seq_len(good), , drop = FALSE]
  boot_df <- as.data.frame(boot_mat)
  
  if (nrow(boot_df) == 0) stop("No successful bootstrap fits (N.good.boots = 0).")
  
  # ---------------------------------------------------------------------------
  # 5) Summary
  # ---------------------------------------------------------------------------
  stat_names <- c("alpha", "beta", "beta/alpha", "intercept", "PrMMinf", "PrUUinf")
  
  se_vec <- c(
    sd(boot_df$alpha, na.rm = TRUE),
    sd(boot_df$beta, na.rm = TRUE),
    sd(boot_df$beta / boot_df$alpha, na.rm = TRUE),
    sd(boot_df$intercept, na.rm = TRUE),
    sd(boot_df$PrMMinf, na.rm = TRUE),
    sd(boot_df$PrUUinf, na.rm = TRUE)
  )
  
  ci_mat <- rbind(
    quantile(boot_df$alpha, probs = c(0.025, 0.975), na.rm = TRUE),
    quantile(boot_df$beta, probs = c(0.025, 0.975), na.rm = TRUE),
    quantile(boot_df$beta / boot_df$alpha, probs = c(0.025, 0.975), na.rm = TRUE),
    quantile(boot_df$intercept, probs = c(0.025, 0.975), na.rm = TRUE),
    quantile(boot_df$PrMMinf, probs = c(0.025, 0.975), na.rm = TRUE),
    quantile(boot_df$PrUUinf, probs = c(0.025, 0.975), na.rm = TRUE)
  )
  
  se_out <- cbind(SE = se_vec, CI2.5 = ci_mat[, 1], CI97.5 = ci_mat[, 2])
  rownames(se_out) <- stat_names
  
  out <- list(
    standard.errors = se_out,
    boot.base = base[1, , drop = FALSE],
    settings = settings,
    N.boots = Nboot,
    N.good.boots = good,
    boot.results = boot_df,
    model = fit$model
  )
  
  dput(out, file = paste0(out.dir, "/", out.name, ".Rdata"))
  return(out)
}
