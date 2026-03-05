
#rm(list=ls())



ABneutralHaploid <- function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{
  allow.neg.intercept <- "no"
  
  # ---- Basic checks ---------------------------------------------------------
  if (is.null(p0uu) || is.null(eqp)) {
    stop("Both eqp value AND p0uu have to be supplied")
  }
  
  p0mm <- 1 - p0uu
  p0um <- 0  # kept for settings output compatibility
  
  if (!isTRUE(all.equal(sum(c(p0mm, p0um, p0uu), na.rm = TRUE), 1))) {
    stop("The initial state probabilities don't sum to 1")
  }
  
  pedigree <- as.matrix(pedigree.data)
  
  # ---- Helpers (2-state) ----------------------------------------------------
  make_genmatrix <- function(alpha, beta) {
    # State order: (UU, MM)
    matrix(c(
      1 - alpha, alpha,
      beta,      1 - beta
    ), nrow = 2, byrow = TRUE)
  }
  
  pair_divergence <- function(p1, p2) {
    # Probability that two tips differ (U vs M)
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
    
    # Basis rows = start state UU / MM (2-state)
    basis <- rbind(
      UU = c(1, 0),
      MM = c(0, 1)
    )
    
    Dt1t2 <- numeric(NROW(pedigree_mat))
    
    for (i in seq_len(NROW(pedigree_mat))) {
      M0 <- pow[[as.character(t0[i])]]
      M1 <- pow[[as.character(dt1[i])]]
      M2 <- pow[[as.character(dt2[i])]]
      
      svt0 <- drop(svGzero %*% M0)  # length-2
      P1   <- basis %*% M1          # 2x2
      P2   <- basis %*% M2
      
      dUU <- pair_divergence(P1["UU", ], P2["UU", ])
      dMM <- pair_divergence(P1["MM", ], P2["MM", ])
      
      Dt1t2[i] <- svt0[1] * dUU + svt0[2] * dMM
    }
    
    # Equilibrium Pr(UU) for 2-state chain:
    puu_inf <- beta / (alpha + beta)
    
    list(puuinf = puu_inf, Dt1t2 = Dt1t2)
  }
  
  # ---- Objective function ---------------------------------------------------
  LSE_intercept <- function(param_int) {
    alpha     <- param_int[1]
    beta      <- param_int[2]
    intercept <- param_int[3]
    
    d <- calc_Dt1t2(pedigree, p0uu = p0uu, alpha = alpha, beta = beta)
    
    sum((pedigree[, 4] - intercept - d$Dt1t2)^2) +
      eqp.weight * NROW(pedigree) * ((d$puuinf - eqp)^2)
  }
  
  # Robust extraction of parameter values from optimx output
  get_pars_from_optimx <- function(opt_df, par_names) {
    if (all(par_names %in% names(opt_df))) {
      return(as.numeric(opt_df[1, par_names]))
    }
    num_cols <- which(vapply(opt_df, is.numeric, logical(1)))
    if (length(num_cols) < length(par_names)) {
      stop("optimx output does not contain enough numeric columns to recover parameters.")
    }
    as.numeric(opt_df[1, num_cols[seq_along(par_names)]])
  }
  
  # ---- Optimization loop ----------------------------------------------------
  optim.method <- "Nelder-Mead"
  final_list <- vector("list", Nstarts)
  
  par_names <- c("alpha", "beta", "intercept")
  
  for (s in seq_len(Nstarts)) {
    alpha.start     <- 10^(stats::runif(1, log10(10^-9), log10(10^-2)))
    beta.start      <- 10^(stats::runif(1, log10(10^-9), log10(10^-2)))
    intercept.start <- stats::runif(1, 0, max(pedigree[, 4]))
    
    # Name the parameter vector so optimx returns named columns
    param_int0 <- c(
      alpha = alpha.start,
      beta = beta.start,
      intercept = intercept.start
    )
    
    message("Progress: ", s / Nstarts, "\n")
    
    opt_out <- suppressWarnings(
      optimx::optimx(par = param_int0, fn = LSE_intercept, method = optim.method)
    )
    
    opt <- as.data.frame(opt_out)
    if (NROW(opt) < 1) stop("optimx returned no rows (optimization failed).")
    
    opt <- opt[1, , drop = FALSE]
    
    pars <- get_pars_from_optimx(opt, par_names)
    names(pars) <- par_names
    
    alphafinal <- pars["alpha"]
    betfinal   <- pars["beta"]
    
    PrMMinf <- alphafinal / (alphafinal + betfinal)
    PrUUinf <- betfinal   / (alphafinal + betfinal)
    
    opt$PrMMinf <- PrMMinf
    opt$PrUUinf <- PrUUinf
    
    opt$alpha.start     <- alpha.start
    opt$beta.start      <- beta.start
    opt$intercept.start <- intercept.start
    
    for (nm in par_names) {
      if (!(nm %in% names(opt))) opt[[nm]] <- pars[[nm]]
    }
    
    final_list[[s]] <- opt
  }
  
  final <- do.call(rbind, final_list)
  
  # ---- Post-fit: LSQ part only (no equilibrium term) ------------------------
  lsqpart <- numeric(NROW(final))
  for (l in seq_len(NROW(final))) {
    d <- calc_Dt1t2(
      pedigree, p0uu = p0uu,
      alpha = as.numeric(final[l, "alpha"]),
      beta  = as.numeric(final[l, "beta"])
    )
    lsqpart[l] <- sum((pedigree[, 4] - as.numeric(final[l, "intercept"]) - d$Dt1t2)^2)
  }
  final$value.part <- lsqpart
  
  if ("value" %in% names(final)) {
    final <- final[order(final[, "value"]), , drop = FALSE]
  }
  
  # ---- Flag/filter results --------------------------------------------------
  alpha_ok <- final[, "alpha"] > 0
  beta_ok  <- final[, "beta"] > 0
  conv_ok  <- if ("convcode" %in% names(final)) final[, "convcode"] == 0 else rep(TRUE, NROW(final))
  
  if (allow.neg.intercept == "yes") {
    keep <- which(alpha_ok & beta_ok & conv_ok)
  } else {
    keep <- which(alpha_ok & beta_ok & (final[, "intercept"] > 0) & conv_ok)
  }
  
  final.1 <- final[keep, , drop = FALSE]
  final.2 <- final[setdiff(seq_len(NROW(final)), keep), , drop = FALSE]
  
  # ---- Predicted values using best model -----------------------------------
  best <- final.1[1, , drop = FALSE]
  
  d_best <- calc_Dt1t2(
    pedigree, p0uu = p0uu,
    alpha = as.numeric(best[1, "alpha"]),
    beta  = as.numeric(best[1, "beta"])
  )
  
  Dt1t2 <- d_best$Dt1t2
  intercept <- as.numeric(best[1, "intercept"])
  Residual <- pedigree[, 4] - intercept - Dt1t2
  
  delta.t <- pedigree[, 2] + pedigree[, 3] - 2 * pedigree[, 1]
  pedigree <- cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
  colnames(pedigree)[c(4, 5, 6, 7)] <- c("div.obs", "delta.t", "div.pred", "residual")
  
  # ---- Info about settings --------------------------------------------------
  info <- c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
  info2 <- c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
  info.out <- data.frame(Para = info, Setting = info2)
  
  # ---- Generating theoretical fit grid -------------------------------------
  alpha <- as.numeric(best[1, "alpha"])
  beta  <- as.numeric(best[1, "beta"])
  
  time1 <- seq(1, max(c(pedigree[, 2], pedigree[, 3])))
  time2 <- seq(1, max(c(pedigree[, 2], pedigree[, 3])))
  
  time.out <- expand.grid(time1, time2)
  time0 <- rep(0, nrow(time.out))
  
  pedigree.new <- as.matrix(cbind(time0, time.out))
  pedigree.new <- cbind(pedigree.new, pedigree.new[, 2] + pedigree.new[, 3] - 2 * pedigree.new[, 1])
  pedigree.new <- pedigree.new[!duplicated(pedigree.new[, 4]), , drop = FALSE]
  pedigree.new <- pedigree.new[, 1:3, drop = FALSE]
  
  d_sim <- calc_Dt1t2(
    pedigree.new, p0uu = p0uu,
    alpha = alpha, beta = beta
  )
  
  pedigree.new <- cbind(
    pedigree.new,
    d_sim$Dt1t2 + intercept,
    pedigree.new[, 2] + pedigree.new[, 3] - 2 * pedigree.new[, 1]
  )
  colnames(pedigree.new) <- c("time0", "time1", "time2", "div.sim", "delta.t")
  pedigree.new <- pedigree.new[order(pedigree.new[, "delta.t"]), , drop = FALSE]
  
  # ---- Output ---------------------------------------------------------------
  model <- "ABneutralHaploid.R"
  
  abfree.out <- list(final.1, final.2, pedigree, info.out, model, pedigree.new)
  names(abfree.out) <- c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")
  
  dput(abfree.out, paste0(out.dir, "/", out.name, ".Rdata"))
  
  return(abfree.out)
}



ABneutralSOMAsimple <- function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{
  allow.neg.intercept <- "no"
  
  # ---- Basic checks ---------------------------------------------------------
  if (is.null(p0uu) || is.null(eqp)) {
    stop("Both eqp value AND p0uu have to be supplied")
  }
  
  p0mm <- 1 - p0uu
  p0um <- 0
  
  if (!isTRUE(all.equal(sum(c(p0mm, p0um, p0uu), na.rm = TRUE), 1))) {
    stop("The initial state probabilities don't sum to 1")
  }
  
  pedigree <- as.matrix(pedigree.data)
  
  # ---- Helpers (compact + reusable) ----------------------------------------
  make_genmatrix <- function(alpha, beta) {
    matrix(c(
      (1 - alpha)^2,            2 * (1 - alpha) * alpha,  alpha^2,
      beta * (1 - alpha),       (1 - alpha) * (1 - beta) + alpha * beta, alpha * (1 - beta),
      beta^2,                   2 * (1 - beta) * beta,   (1 - beta)^2
    ), nrow = 3, byrow = TRUE)
  }
  
  pair_divergence <- function(p1, p2) {
    0.5 * (p1[1] * p2[2] + p1[2] * p2[1] + p1[2] * p2[3] + p1[3] * p2[2]) +
      1.0 * (p1[1] * p2[3] + p1[3] * p2[1])
  }
  
  calc_Dt1t2 <- function(pedigree_mat, p0uu, alpha, beta, weight) {
    p0mm <- 1 - p0uu
    svGzero <- c(p0uu, weight * p0mm, (1 - weight) * p0mm)  # (UU, UM-ish, MM-ish)
    
    G <- make_genmatrix(alpha, beta)
    
    t0  <- as.numeric(pedigree_mat[, 1])
    dt1 <- as.numeric(pedigree_mat[, 2] - pedigree_mat[, 1])
    dt2 <- as.numeric(pedigree_mat[, 3] - pedigree_mat[, 1])
    
    exps <- sort(unique(c(t0, dt1, dt2)))
    pow  <- setNames(lapply(exps, function(k) expm::`%^%`(G, k)), as.character(exps))
    
    basis <- rbind(
      UU = c(1, 0, 0),
      UM = c(0, 1, 0),
      MM = c(0, 0, 1)
    )
    
    Dt1t2 <- numeric(NROW(pedigree_mat))
    
    for (i in seq_len(NROW(pedigree_mat))) {
      M0 <- pow[[as.character(t0[i])]]
      M1 <- pow[[as.character(dt1[i])]]
      M2 <- pow[[as.character(dt2[i])]]
      
      svt0 <- drop(svGzero %*% M0)     # length-3
      P1   <- basis %*% M1             # 3x3, rows = start-state
      P2   <- basis %*% M2
      
      dUU <- pair_divergence(P1["UU", ], P2["UU", ])
      dUM <- pair_divergence(P1["UM", ], P2["UM", ])
      dMM <- pair_divergence(P1["MM", ], P2["MM", ])
      
      Dt1t2[i] <- svt0[1] * dUU + svt0[2] * dUM + svt0[3] * dMM
    }
    
    puu_inf <- (beta^2) / ((alpha + beta)^2)
    list(puuinf = puu_inf, Dt1t2 = Dt1t2)
  }
  
  # ---- Objective function ---------------------------------------------------
  LSE_intercept <- function(param_int) {
    alpha     <- param_int[1]
    beta      <- param_int[2]
    weight    <- param_int[3]
    intercept <- param_int[4]
    
    d <- calc_Dt1t2(pedigree, p0uu = p0uu, alpha = alpha, beta = beta, weight = weight)
    
    sum((pedigree[, 4] - intercept - d$Dt1t2)^2) +
      eqp.weight * NROW(pedigree) * ((d$puuinf - eqp)^2)
  }
  
  # Robust extraction of parameter values from optimx output
  get_pars_from_optimx <- function(opt_df, par_names) {
    # If optimx preserved parameter names, use them
    if (all(par_names %in% names(opt_df))) {
      return(as.numeric(opt_df[1, par_names]))
    }
    # Fallback: take the first length(par_names) numeric columns
    num_cols <- which(vapply(opt_df, is.numeric, logical(1)))
    if (length(num_cols) < length(par_names)) {
      stop("optimx output does not contain enough numeric columns to recover parameters.")
    }
    as.numeric(opt_df[1, num_cols[seq_along(par_names)]])
  }
  
  # ---- Optimization loop ----------------------------------------------------
  optim.method <- "Nelder-Mead"
  final_list <- vector("list", Nstarts)
  
  par_names <- c("alpha", "beta", "weight", "intercept")
  
  for (s in seq_len(Nstarts)) {
    alpha.start     <- 10^(stats::runif(1, log10(10^-9), log10(10^-2)))
    beta.start      <- 10^(stats::runif(1, log10(10^-9), log10(10^-2)))
    weight.start    <- stats::runif(1, 0, 0.5)
    intercept.start <- stats::runif(1, 0, max(pedigree[, 4]))
    
    # IMPORTANT: name the parameter vector so optimx returns named columns
    param_int0 <- c(
      alpha = alpha.start,
      beta = beta.start,
      weight = weight.start,
      intercept = intercept.start
    )
    
    message("Progress: ", s / Nstarts, "\n")
    
    opt_out <- suppressWarnings(
      optimx::optimx(par = param_int0, fn = LSE_intercept, method = optim.method)
    )
    
    opt <- as.data.frame(opt_out)
    if (NROW(opt) < 1) stop("optimx returned no rows (optimization failed).")
    
    # Keep first row (as in your original)
    opt <- opt[1, , drop = FALSE]
    
    pars <- get_pars_from_optimx(opt, par_names)
    names(pars) <- par_names
    
    alphafinal <- pars["alpha"]
    betfinal   <- pars["beta"]
    
    PrMMinf <- (alphafinal^2) / ((alphafinal + betfinal)^2)
    PrUMinf <- (2 * alphafinal * betfinal) / ((alphafinal + betfinal)^2)
    PrUUinf <- (betfinal^2) / ((alphafinal + betfinal)^2)
    
    opt$PrMMinf <- PrMMinf
    opt$PrUMinf <- PrUMinf
    opt$PrUUinf <- PrUUinf
    
    opt$alpha.start     <- alpha.start
    opt$beta.start      <- beta.start
    opt$weight.start    <- weight.start
    opt$intercept.start <- intercept.start
    
    # If parameter columns weren't present, add them explicitly (for downstream code)
    for (nm in par_names) {
      if (!(nm %in% names(opt))) opt[[nm]] <- pars[[nm]]
    }
    
    final_list[[s]] <- opt
  }
  
  final <- do.call(rbind, final_list)
  
  # ---- Post-fit: compute "value.part" (LSQ without equilibrium term) --------
  lsqpart <- numeric(NROW(final))
  for (l in seq_len(NROW(final))) {
    d <- calc_Dt1t2(
      pedigree, p0uu = p0uu,
      alpha  = as.numeric(final[l, "alpha"]),
      beta   = as.numeric(final[l, "beta"]),
      weight = as.numeric(final[l, "weight"])
    )
    lsqpart[l] <- sum((pedigree[, 4] - as.numeric(final[l, "intercept"]) - d$Dt1t2)^2)
  }
  final$value.part <- lsqpart
  
  if ("value" %in% names(final)) {
    final <- final[order(final[,"value"]), , drop = FALSE]
  }
  
  # ---- Flag/filter results --------------------------------------------------
  alpha_ok <- final[, "alpha"] > 0
  beta_ok  <- final[, "beta"] > 0
  conv_ok  <- if ("convcode" %in% names(final)) final[, "convcode"] == 0 else rep(TRUE, NROW(final))
  
  if (allow.neg.intercept == "yes") {
    keep <- which(alpha_ok & beta_ok & conv_ok)
  } else {
    keep <- which(alpha_ok & beta_ok & (final[, "intercept"] > 0) & conv_ok)
  }
  
  final.1 <- final[keep, , drop = FALSE]
  final.2 <- final[setdiff(seq_len(NROW(final)), keep), , drop = FALSE]
  
  # ---- Predicted values using best model -----------------------------------
  best <- final.1[1, , drop = FALSE]
  
  d_best <- calc_Dt1t2(
    pedigree, p0uu = p0uu,
    alpha  = as.numeric(best[1, "alpha"]),
    beta   = as.numeric(best[1, "beta"]),
    weight = as.numeric(best[1, "weight"])
  )
  
  Dt1t2 <- d_best$Dt1t2
  intercept <- as.numeric(best[1, "intercept"])
  Residual <- pedigree[, 4] - intercept - Dt1t2
  
  delta.t <- pedigree[, 2] + pedigree[, 3] - 2 * pedigree[, 1]
  pedigree <- cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
  colnames(pedigree)[c(4, 5, 6, 7)] <- c("div.obs", "delta.t", "div.pred", "residual")
  
  # ---- Info about settings --------------------------------------------------
  info <- c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
  info2 <- c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
  info.out <- data.frame(Para = info, Setting = info2)
  
  # ---- Generating theoretical fit grid -------------------------------------
  alpha <- as.numeric(best[1, "alpha"])
  beta  <- as.numeric(best[1, "beta"])
  weight <- as.numeric(best[1, "weight"])
  
  time1 <- seq(1, max(c(pedigree[, 2], pedigree[, 3])))
  time2 <- seq(1, max(c(pedigree[, 2], pedigree[, 3])))
  
  time.out <- expand.grid(time1, time2)
  time0 <- rep(0, nrow(time.out))
  
  pedigree.new <- as.matrix(cbind(time0, time.out))
  pedigree.new <- cbind(pedigree.new, pedigree.new[, 2] + pedigree.new[, 3] - 2 * pedigree.new[, 1])
  pedigree.new <- pedigree.new[!duplicated(pedigree.new[, 4]), , drop = FALSE]
  pedigree.new <- pedigree.new[, 1:3, drop = FALSE]
  
  d_sim <- calc_Dt1t2(
    pedigree.new, p0uu = p0uu,
    alpha = alpha, beta = beta, weight = weight
  )
  
  pedigree.new <- cbind(
    pedigree.new,
    d_sim$Dt1t2 + intercept,
    pedigree.new[, 2] + pedigree.new[, 3] - 2 * pedigree.new[, 1]
  )
  colnames(pedigree.new) <- c("time0", "time1", "time2", "div.sim", "delta.t")
  pedigree.new <- pedigree.new[order(pedigree.new[, "delta.t"]), , drop = FALSE]
  
  # ---- Output ---------------------------------------------------------------
  model <- "ABneutralSOMAsimple.R"
  
  abfree.out <- list(final.1, final.2, pedigree, info.out, model, pedigree.new)
  names(abfree.out) <- c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")
  
  dput(abfree.out, paste0(out.dir, "/", out.name, ".Rdata"))
  
  return(abfree.out)
}


### Model comparison function
#This needs to be modified with respect to the original AlpaBeta function
#As in haploid model is there is one parameter less
#Note: does not currently work for non-neutral models
FtestRSSHaploid <- function(pedigree.select, pedigree.null) 
{
    est <- dget(pedigree.select)
    estN <- dget(pedigree.null)
    if (estN$model == "ABnull.R") {
        RSSf <- est$estimates[1, "value"]
        RSSr <- sum((estN$pedigree[, "residual"])^2)
        Npara_r <- 1
        Npara_f <- 4
        dfF <- length(est$pedigree[, "residual"]) - 4
        dfR <- length(estN$pedigree[, "residual"]) - 1
        dfN <- dfR - dfF
        Fvalue <- ((RSSr - RSSf)/(Npara_f - Npara_r))/(RSSf/dfF)
        pvalue <- pf(Fvalue, dfN, dfF, lower.tail = FALSE)
        output <- c(RSSf, RSSr, dfF, dfR, Fvalue, pvalue)
        names(output) <- c("RSS_F", "RSS_R", "df_F", "df_R", 
            "Fvalue", "pvalue")
        outfinal <- list(output, est$estimates, estN$estimates)
        names(outfinal) <- c("Ftest", "est.selection", "est.neutral")
    }
    if (estN$model != "ABnull.R") {
        RSSf <- est$estimates[1, "value"]
        RSSr <- estN$estimates[1, "value"]
        Npara_r <- 3
        Npara_f <- 4
        dfF <- length(est$pedigree[, "residual"]) - 4
        dfR <- length(estN$pedigree[, "residual"]) - 3
        dfN <- dfR - dfF
        Fvalue <- ((RSSr - RSSf)/(Npara_f - Npara_r))/(RSSf/dfF)
        pvalue <- pf(Fvalue, dfN, dfF, lower.tail = FALSE)
        output <- c(RSSf, RSSr, dfF, dfR, Fvalue, pvalue)
        names(output) <- c("RSS_F", "RSS_R", "df_F", "df_R", 
            "Fvalue", "pvalue")
        outfinal <- list(output, est$estimates, estN$estimates)
        names(outfinal) <- c("Ftest", "est.selection", "est.neutral")
    }
    outfinal
}


