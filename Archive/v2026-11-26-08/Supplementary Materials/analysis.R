# Script containing code written in R (v4.4.3 [2025-02-28 ucrt]; R Core Team, 2025) to reproduce the analyses used by the paper ’Bayesian Factor Score Reliability: Conceptualizing a Bayesian Measure of Factor Score Indeterminacy’ and its supplementary materials.

# Author: Emil Niclas Meyer-Hansen (email: emil098meyerhansen@gmail.com)
# Date: 2026/11/26
# Copyright: 2026 The author. License: Except where otherwise indicated, all contents of this document and associated files are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). All data and software, including but not necessarily limited to, source code, executable code, code snippets, code chunks, algorithms, and/or scripts, attributable to this document and/or its associated files are expressly excluded from the foregoing license, and unless otherwise indicated, are instead licensed under the GNU General Public License, version 3 (GPL-3.0) but permission is granted to re-license these materials under later versions. The processed data sampled from the Open-Source Psychometrics Project (2014; provided in Petras & Dantlgraber, 2022) is excluded from the foregoing licenses and is instead licensed under the GNU General Public License, version 2 (GPL-2.0). By engaging with this document and/or any associated files, which include, but are not necessarily limited to, downloading, using, viewing, and/or distributing any of them, in parts of whole, you agree to comply with the applicable license terms for the respective content types.


# NOTE: While this file contains the code to reproduce all analyses, for a user to reproduce the analyses themselves, they may need to install the software listed below (e.g., R and the required R packages) and change file paths to the specified data and models.
# For more details on software specifications, see the supplementary materials.
# Due to the complexity of Bayesian factor analysis models, users should expect that the step involving model fitting will take several hours.
# Comments have been added throughout to help link the output of the code to the values/figures/tables reported in the main paper and supplementary materials.





# Setup -------------------------------------------------------------------

# Load (or install) the required R packages
library(pacman)
pacman::p_load(
  tidyverse,
  brms,
  parallel,
  doParallel,
  foreach,
  random,
  Polychrome,
  posterior,
  einsum,
  tikzDevice,
  knitr,
  IPV,
  nleqslv,
  ggdist,
  bayesplot,
  lavaan,
  blavaan,
  scrutiny,
  weights,
  MASS,
  datawizard,
  rstudioapi,
  ggbreak,
  effectsize,
  ggplot2,
  readr,
  latex2exp,
  bayestestR,
  kableExtra,
  openxlsx,
  trialr
)


# Set working directory
# setwd("")


# Specify tikz options
options(
  tikzDefaultEngine = "pdftex",
  tikzLatex = "pdflatex",
  tikzLatexPackages = c(
    "\\usepackage[utf8]{inputenc}",
    "\\usepackage[T1]{fontenc}",
    "\\usepackage{lmodern}",
    "\\usepackage{amsmath,amssymb}",
    "\\usepackage{amsfonts}",
    "\\usepackage{graphicx}",
    "\\usepackage{xcolor}",
    "\\usepackage{newunicodechar}",
    "\\newunicodechar{−}{-}",
    "\\newunicodechar{;}{;}",
    "\\usepackage{pgfplots}"
  )
)


# R and Rstudio version
r_version <- sessionInfo()$R.version$version.string
r_version <- stringr::str_replace(r_version, "R version ", "v")
r_version <- stringr::str_replace(r_version, "\\(", "\\[")
r_version <- stringr::str_replace(r_version, "\\)", "\\]")
rstudio_version <- paste0('v', rstudioapi::versionInfo()$version)


# Document Options
knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(fig.pos = "H", out.extra = "")
digits <- list(
  text = 3,
  figure = 3,
  table = 3
)
options(digits = digits$figure)
big_mark <- ","


# Specify cores
cores <- chains <- max(c(1, parallel::detectCores()-1)) # NOTE: 15 cores were used for the analysis
doParallel::registerDoParallel(cores = cores)
options(mc.cores = cores)
rstan_options(auto_write = TRUE)
posterior_draws <- 40000
warmup_draws <- 4000
total_draws <- ceiling( ( posterior_draws/chains ) + warmup_draws )
posterior_draws <- (total_draws - warmup_draws) * chains


# Specify sample size
sample_size <- 2000


# Specify NHST
alpha <- .05
beta <- .2


# Specify RNG 
SEED <- 322047077 # random::randomNumbers(1, min = 1, max = 1000000000, col = 1) %>% as.vector() # NOTE: Commented-out function used to generate the specified seed. Promotes reproducibility but does not guarantee it.
set.seed(SEED)


# Specify colors used for traceplots
traceplot_color_chart <- Polychrome::palette36.colors(cores) %>% as.vector()
golden_ratio <- (1 + sqrt(5)) / 2
fig_dpi <- 320


# Define QoL Formatting function
qformat <- function(x = NULL, scientific = FALSE, trim = TRUE, big.mark = ","){
  # Exception handling
  stopifnot(length(x) > 0)
  stopifnot(is.numeric(x))
  stopifnot(length(scientific) == 1)
  stopifnot(is.logical(scientific))
  stopifnot(length(trim) == 1)
  stopifnot(is.logical(trim))
  stopifnot(length(big.mark) == 1)
  stopifnot(is.character(big.mark))
  
  return(format(x, scientific = scientific, trim = trim, big.mark = big_mark))
}


# Custom functions for conservative rounding
round_down2 <- function(x = NULL, digits = 3, add = FALSE, threshold = 0){
  require(scrutiny)
  require(weights)
  stopifnot(length(threshold) == 1)
  stopifnot(is.numeric(threshold))
  stopifnot(is.finite(threshold))
  x_rounded <- ifelse(
    is.na(x), NA, ifelse(
      x < threshold, weights::rd(scrutiny::round_up(x = x, digits = digits), digits = digits, add = add, max = digits), weights::rd(scrutiny::round_down(x = x, digits = digits), digits = digits, add = add, max = digits)
    )
  )
  return(x_rounded)
}
round_up2 <- function(x = NULL, digits = 3, add = FALSE, threshold = 0){
  require(scrutiny)
  require(weights)
  stopifnot(length(threshold) == 1)
  stopifnot(is.numeric(threshold))
  stopifnot(is.finite(threshold))
  x_rounded <- ifelse(
    is.na(x), NA, ifelse(
      x > threshold, weights::rd(scrutiny::round_up(x = x, digits = digits), digits = digits, add = add, max = digits), weights::rd(scrutiny::round_down(x = x, digits = digits), digits = digits, add = add, max = digits)
    )
  )
  return(x_rounded)
}


# Custom function for getting the critical value from a normal (or Student t) distribution
critical.value <- function(alpha = NULL, df = Inf){
  stopifnot(length(alpha) == 1)
  stopifnot(is.numeric(alpha))
  stopifnot(alpha > 0)
  stopifnot(alpha < 1)
  stopifnot(length(df) == 1)
  stopifnot(is.numeric(df) | is.integer(df))
  stopifnot(df >= 1)
  
  if(is.finite(df)){
    return(qt(1-alpha/2, df = df))
  }
  else(
    return(qnorm(1-alpha/2))
  )
}


# Define standard Maximum-likelihood functions
var2 <- function(x = NULL, x_mean = NULL, na.rm = FALSE){
  stopifnot(is.numeric(x))
  stopifnot(is.logical(na.rm))
  if(!is.null(x_mean)){
    stopifnot(is.numeric(x_mean))
  }else{
    x_mean <- mean(x, na.rm = na.rm)
  }
  x_var <- mean((x - x_mean)^2, na.rm = na.rm)
  return(x_var)
}
sd2 <- function(x = NULL, x_mean = NULL, na.rm = FALSE){
  stopifnot(is.numeric(x))
  stopifnot(is.logical(na.rm))
  if(!is.null(x_mean)){
    stopifnot(is.numeric(x_mean))
  }else{
    x_mean <- mean(x, na.rm = na.rm)
  }
  x_sd <- sqrt(var2(x = x, x_mean = x_mean, na.rm = FALSE))
  return(x_sd)
}
cov2 <- function(x = NULL, y = NULL, x_mean = NULL, y_mean = NULL, na.rm = FALSE){
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(is.logical(na.rm))
  if(!is.null(x_mean)){
    stopifnot(is.numeric(x_mean))
  }else{
    x_mean <- mean(x, na.rm = na.rm)
  }
  if(!is.null(y_mean)){
    stopifnot(is.numeric(y_mean))
  }else{
    y_mean <- mean(y, na.rm = na.rm)
  }
  xy_cov <- mean((x - x_mean)*(y - y_mean), na.rm = na.rm)
  return(xy_cov)
}
cor2 <- function(x = NULL, y = NULL, x_mean = NULL, y_mean = NULL, x_sd = NULL, y_sd = NULL, na.rm = FALSE){
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(is.logical(na.rm))
  if(!is.null(x_mean)){
    stopifnot(is.numeric(x_mean))
  }else{
    x_mean <- mean(x, na.rm = na.rm)
  }
  if(!is.null(y_mean)){
    stopifnot(is.numeric(y_mean))
  }else{
    y_mean <- mean(y, na.rm = na.rm)
  }
  if(!is.null(x_sd)){
    stopifnot(is.numeric(x_sd))
    stopifnot(x_sd > 0)
  }else{
    x_sd <- sd2(x = x, x_mean = x_mean, na.rm = na.rm)
  }
  if(!is.null(y_sd)){
    stopifnot(is.numeric(y_sd))
    stopifnot(y_sd > 0)
  }else{
    y_sd <- sd2(x = y, x_mean = y_mean, na.rm = na.rm)
  }
  xy_cov <- cov2(x, y, x_mean = x_mean, y_mean = y_mean, na.rm = na.rm)
  xy_cor <- xy_cov/(x_sd*y_sd)
  return(xy_cor)
}
mse2 <- function(y = NULL, y_pred = NULL, na.rm = FALSE){
  stopifnot(is.numeric(y))
  stopifnot(is.numeric(y_pred))
  stopifnot(is.logical(na.rm))
  return(mean((y - y_pred)^2, na.rm = na.rm))
}
rmse2 <- function(y = NULL, y_pred = NULL, na.rm = FALSE){
  stopifnot(is.numeric(y))
  stopifnot(is.numeric(y_pred))
  stopifnot(is.logical(na.rm))
  return(sqrt(mse2(y = y, y_pred = y_pred, na.rm = na.rm)))
}
standardize2 <- function(x = NULL, x_mean = NULL, x_sd = NULL, na.rm = FALSE){
  stopifnot(is.numeric(x))
  stopifnot(is.logical(na.rm))
  if(!is.null(x_mean)){
    stopifnot(is.numeric(x_mean))
  }else{
    x_mean <- mean(x, na.rm = na.rm)
  }
  if(!is.null(x_sd)){
    stopifnot(is.numeric(x_sd))
    stopifnot(x_sd > 0)
  }else{
    x_sd <- sd2(x = x, x_mean = x_mean, na.rm = na.rm)
  }
  x_standardized <- (x - x_mean) / x_sd
  return(x_standardized)
}


# Define converter functions for the FSR
fsr_to_fsr_adj <- function(fsr = NULL, n = NULL, p = NULL){
  if(is.null(fsr)) stop("Argument 'fsr' cannot be NULL.")
  if(!length(fsr) > 0) stop("Argument 'fsr' cannot must have positive length.")
  if(!is.numeric(fsr)) stop("Argument 'fsr' must be numeric.")
  if(!all(fsr >= 0)) stop("Argument 'fsr' cannot be negative.")
  if(!all(fsr <= 1)) stop("Argument 'fsr' cannot exceed unity.")
  
  if(is.null(n)) stop("Argument 'n' cannot be NULL.")
  if(!length(n) == 1 & length(n) != length(fsr)) stop("Argument 'n' must have length one or length equal to argument 'fsr'.")
  if(!is.numeric(n)) stop("Argument 'n' must be numeric.")
  n <- round(n)
  if(!all(n > 0)) stop("Argument 'n' must be positive.")
  
  if(is.null(p)) stop("Argument 'p' cannot be NULL.")
  if(!length(p) == 1 & length(p) != length(fsr)) stop("Argument 'p' must have length one or length equal to argument 'fsr'.")
  if(!is.numeric(p)) stop("Argument 'p' must be numeric.")
  if(!all(p >= 1)) stop("Argument 'p' must be positive.")
  if(!all(n > p)) stop("Argument 'n' must exceed argument 'p'.")
  
  adj_factor <- (n - 1) / (n - p - 1)
  fsr_adj <- 1 - (1 - fsr) * adj_factor
  return(fsr_adj)
}
fsr_to_mcc <- function(fsr = NULL){
  if(is.null(fsr)) stop("Argument 'fsr' cannot be NULL.")
  if(!length(fsr) > 0) stop("Argument 'fsr' cannot must have positive length.")
  if(!is.numeric(fsr)) stop("Argument 'fsr' must be numeric.")
  if(!all(fsr >= 0)) stop("Argument 'fsr' cannot be negative.")
  if(!all(fsr <= 1)) stop("Argument 'fsr' cannot exceed unity.")
  mcc <- 2*fsr - 1
  return(mcc)
}
mcc_to_fsr <- function(mcc = NULL){
  if(is.null(mcc)) stop("Argument 'mcc' cannot be NULL.")
  if(!length(mcc) > 0) stop("Argument 'mcc' cannot must have positive length.")
  if(!is.numeric(mcc)) stop("Argument 'mcc' must be numeric.")
  if(!all(abs(mcc) <= 1)) stop("Argument 'mcc' cannot exceed unity in absolute value.")
  fsr <- (mcc + 1)/2
  return(fsr)
}
fsr_to_vif <- function(fsr = NULL){
  if(is.null(fsr)) stop("Argument 'fsr' cannot be NULL.")
  if(!length(fsr) > 0) stop("Argument 'fsr' cannot must have positive length.")
  if(!is.numeric(fsr)) stop("Argument 'fsr' must be numeric.")
  if(!all(fsr >= 0)) stop("Argument 'fsr' cannot be negative.")
  if(!all(fsr <= 1)) stop("Argument 'fsr' cannot exceed unity.")
  vif <- 1/(1-fsr)
  return(vif)
}


# Define function to create results for summary table
get_summary_table <- function(x = NULL, ci = 0.95){
  require(tidyverse)
  require(bayestestR)
  x_hdi <- vapply(seq_len(ncol(x)), function(j) {
    hj <- bayestestR::hdi(x[, j], ci = ci)
    c(hj$CI_low, hj$CI_high)
  }, numeric(2))  
  output <- paste0(
    apply(x, 2, mean) %>% unname() %>% round_down2(digits = digits$table),
    " (",
    apply(x, 2, sd2) %>% unname() %>% round_up2(digits = digits$table),
    ") [", 
    x_hdi[1,] %>% round_down2(digits = digits$table),
    "; ",
    x_hdi[2,] %>% round_up2(digits = digits$table),
    "]"
  )
  return(output)
}


# Define function for estimating the shape and rate parameter for the gamma distribution based on mode (i.e., 1) and variance (i.e., 0.304^2) constraint
get_gamma_params_from_mode_sd <- function(mode = NULL, sd = NULL, a_init = 1, b_init = 1, verbose = FALSE){
  require(nleqslv)
  f <- function(x){
    a <- x[1]
    b <- x[2]
    c((a - 1)/b - mode,
      a/b^2 - sd^2)
  }
  f_solved <- nleqslv::nleqslv(x = c(a_init, b_init), fn = f)
  if(verbose){
    return(f_solved)
  }
  else{
    c(a = f_solved$x[1], b = f_solved$x[2])
  }
}


# Define d-wise operators for addition, subtraction, and the Hadamard product
dwise <- function(f) function(x, y) {
  stopifnot(is.array(x), length(x) %% length(y) == 0L)
  return(f(x, as.vector(y)))
}
`%+d%` <- dwise(`+`)
`%-d%` <- dwise(`-`)
`%od%` <- dwise(`*`)


# Define the d-wise transpose operator
td <- function(x) return(aperm(x, c(2, 1, 3)))


# Define the d-wise operator for matrix multiplication 
`%*d%` <- local({
  f <- einsum::einsum_generator("ijs, jls->ils")
  function(x, y) {
    stopifnot(dim(x)[2] == dim(y)[1], dim(x)[3] == dim(y)[3])
    structure(f(x, y), dimnames = list(dimnames(x)[[1]], dimnames(y)[[2]], NULL))
  }
})


# Define the d-wise (generalized) inverse matrix operator
invd <- function(x, symmetric = FALSE, tol = .Machine$double.eps^0.5) {
  stopifnot(length(dim(x)) == 3L, dim(x)[1] == dim(x)[2])
  d <- dim(x)[3]
  out <- array(0, dim(x), dimnames = dimnames(x))
  bad <- integer(0)
  for(s in seq_len(d)){
    xs <- x[,, s]
    out[,, s] <- if(rcond(xs) < tol){
      bad <- c(bad, s); MASS::ginv(xs)
    } else if(symmetric) chol2inv(chol(xs)) else solve(xs)
  }
  if(length(bad)){
    warning(sprintf("Singular matrices at %d/%d slices (first: %s). Pseudo-inverse used. Show caution when interpreting the results.", length(bad), d, paste(utils::head(bad, 5), collapse = ", ")))
  }
  return(structure(out, bad = bad))
}
invd <- function(x, symmetric = FALSE, tol = .Machine$double.eps^0.5) {
  require(MASS)
  stopifnot(length(dim(x)) == 3L, dim(x)[1] == dim(x)[2])
  k <- dim(x)[1]
  d <- dim(x)[3]
  out <- array(0, dim(x), dimnames = dimnames(x))
  bad <- integer(0)
  for(s in seq_len(d)){
    xs <- matrix(x[,, s], k, k)
    out[,, s] <- if(rcond(xs) < tol){
      bad <- c(bad, s); MASS::ginv(xs)
    } else if(symmetric) chol2inv(chol(xs)) else solve(xs)
  }
  if(length(bad)){
    warning(sprintf("Singular matrices at %d/%d slices (first: %s). Pseudo-inverse used. Show caution when interpreting the results.", length(bad), d, paste(utils::head(bad, 5), collapse = ", ")))
  }
  return(structure(out, bad = bad))
}

# Define the d-wise square-root operator
sqrtd <- function(x) {
  stopifnot(length(dim(x)) == 3L, dim(x)[1] == dim(x)[2])
  return(array(apply(x, 3, sqrt), dim(x), dimnames(x)))
}


# Define utility function for getting parameter matrices from draws
get_parameter_matrices <- function(brms_fit = NULL, draws = NULL, latent = NULL, overt = NULL, prior = NULL, ndraws = NULL){
  require(brms)
  
  # Input validation
  if(!is.null(draws)){
    if(!is.matrix(draws) & !is.data.frame(draws)) stop("Argument 'draws' must be a matrix or data.frame object.")
    if(!dim(draws)[1] > 0) stop("Argument 'draws' must have a positive number of rows.")
    if(!dim(draws)[2] > 0) stop("Argument 'draws' must have a positive number of columns.")
  }
  else if(!is.null(brms_fit)){
    if(!brms::is.brmsfit(brms_fit)) stop("Argument 'brms_fit' must be a brmsfit object.")
    draws <- brms::as_draws_df(brms_fit)
  }
  else{
    stop("Arguments 'brms_fit' and 'draws' cannot both be NULL.")
  }
  
  if(is.null(latent)) stop("Argument 'latent' cannot be NULL.")
  if(!length(latent) > 0) stop("Argument 'latent' must have positive length.")
  if(!is.character(latent)) stop("Argument 'latent' must be a character vector.")
  
  if(is.null(overt)) stop("Argument 'overt' cannot be NULL.")
  if(!length(overt) > 0) stop("Argument 'overt' must have positive length.")
  if(!is.character(overt)) stop("Argument 'overt' must be a character vector.")
  
  prior <- match.arg(arg = prior, choices = c("no", "yes", "only"))
  
  # Organize input
  overt <- sort(overt)
  latent <- sort(latent)
  overt_util <- overt %>% stringr::str_remove_all("\\_")
  latent_util <- latent %>% stringr::str_remove_all("\\_")
  draws_colnames <- sort(colnames(draws))
  draws <- draws[,draws_colnames]
  
  
  # Input validation of optional arguments
  if(!is.null(ndraws)){
    if(length(ndraws) != 1) stop("Argument 'ndraws' must have a length of one.")
    if(!is.numeric(ndraws)) stop("Argument 'ndraws' must be numeric")
    if(ndraws > nrow(draws)) stop("Argument 'ndraws' cannot be greater than the number of draws specified by the objects in the 'brms_fit' or 'draws' argument.")
    draws <- draws[1:ndraws,]
  }
  
  # Specify number of draws
  o <- length(overt)
  l <- length(latent)
  d <- nrow(draws)
  
  
  # Define utility functions
  col_of <- function(x) as.numeric(draws[[x]])
  fill_coef <- function(nm, cn, d) {
    P   <- prod(dim(nm))
    A   <- array(0, c(dim(nm), d), dimnames = c(dimnames(nm), list(NULL)))
    hit <- which(nm %in% cn)
    if (length(hit)) {
      got <- vapply(nm[hit], col_of, numeric(d))          # d x k
      stopifnot(nrow(got) == d)
      A[outer(hit, (seq_len(d) - 1L) * P, "+")] <- t(got)
    }
    A
  }
  nmat <- function(r, c, f) {
    m <- outer(r, c, f)
    dimnames(m) <- list(r, c)
    m
  }
  fill_cov <- function(nm_sd, nm_cor, cn, nm_sd_col = nm_sd, nm_cor_alt = nm_cor, d) {
    p <- nrow(nm_cor); q <- ncol(nm_cor)
    stopifnot(length(nm_sd) == p, length(nm_sd_col) == q,
              nrow(nm_cor_alt) == q, ncol(nm_cor_alt) == p)
    
    for (nm in list(nm_sd, nm_sd_col)) {
      have <- nm %in% cn
      if (any(!have)) warning("no sigma for: ", paste(nm[!have], collapse = ", "))
    }
    
    sd_r <- matrix(fill_coef(nm_sd, cn, d = d), p, d)             # p x 1 x d  ->  p x d
    sd_c <- matrix(fill_coef(nm_sd_col, cn, d = d), q, d)             # q x 1 x d  ->  q x d
    
    rho <- fill_coef(nm_cor, cn, d = d) + aperm(fill_coef(nm_cor_alt, cn, d = d), c(2, 1, 3))
    if (identical(dimnames(nm_cor)[[1]], dimnames(nm_cor)[[2]]))
      rho <- rho + array(diag(1, p, p), c(p, p, d))
    
    Si <- aperm(array(sd_r, c(p, d, q)), c(1, 3, 2))           # Si[i,j,t] = sd_r[i,t]
    Sj <- aperm(array(sd_c, c(q, d, p)), c(3, 1, 2))           # Sj[i,j,t] = sd_c[j,t]
    S  <- rho * Si * Sj
    dimnames(S) <- dimnames(nm_cor)
    S
  }
  
  if(prior != "only"){
    Psi_Y_X <- fill_cov(
      nmat(latent_util, "1", function(r, c) paste0("sigma_", r)),
      nmat(latent_util, overt_util, function(r, c) paste0("rescor__", r, "__", c)),
      draws_colnames,
      nm_sd_col  = nmat(overt_util, "1", function(r, c) paste0("sigma_", r)),
      nm_cor_alt = nmat(overt_util, latent_util, function(r, c) paste0("rescor__", r, "__", c)),
      d = d
    )
    output_posterior <- list(
      Lambda = fill_coef(nmat(overt_util, latent, function(r, c) paste0("bsp_", r, "_mi", c)), draws_colnames, d = d),
      Beta = fill_coef(nmat(latent_util, latent, function(r, c) paste0("bsp_", r, "_mi", c)), draws_colnames, d = d),
      Omega = fill_coef(nmat(overt_util, overt_util, function(r, c) paste0("b_", r, "_", c)), draws_colnames, d = d),
      Delta = fill_coef(nmat(latent, overt_util, function(r, c) paste0("b_", r, "_", c)), draws_colnames, d = d),
      Psi_Y = fill_cov(nmat(latent_util, "1", function(r, c) paste0("sigma_", r)), nmat(latent_util, latent_util, function(r, c) paste0("rescor__", r, "__", c)), draws_colnames, d = d),
      Psi_X = fill_cov(nmat(overt_util, "1", function(r, c) paste0("sigma_", r)), nmat(overt_util, overt_util, function(r, c) paste0("rescor__", r, "__", c)), draws_colnames, d = d),
      Psi_Y_X = Psi_Y_X,
      Psi_X_Y <- aperm(Psi_Y_X, c(2, 1, 3)), # NOTE: Currently doesn't output anything for some reason
      alpha_X = fill_coef(nmat(overt_util, "1", function(r, c) paste0("b_", r, "_Intercept")), draws_colnames, d = d) |> aperm(c(2, 1, 3)),
      alpha_Y = fill_coef(nmat(latent_util, "1", function(r, c) paste0("b_", r, "_Intercept")), draws_colnames, d = d) |> aperm(c(2, 1, 3))
    )
  }
  if(prior != "no"){
    prior_Psi_Y_X <- fill_cov(
      nmat(latent_util, "1", function(r, c) paste0("prior_sigma_", r)),
      nmat(latent_util, overt_util, function(r, c) paste0("prior_rescor__", r, "__", c)),
      draws_colnames,
      nm_sd_col  = nmat(overt_util, "1", function(r, c) paste0("prior_sigma_", r)),
      nm_cor_alt = nmat(overt_util, latent_util, function(r, c) paste0("prior_rescor__", r, "__", c)),
      d = d
    )
    output_prior <- list(
      prior_Lambda = fill_coef(nmat(overt_util, latent, function(r, c) paste0("prior_bsp_", r, "_mi", c)), draws_colnames, d = d),
      prior_Beta = fill_coef(nmat(latent_util, latent, function(r, c) paste0("prior_bsp_", r, "_mi", c)), draws_colnames, d = d),
      prior_Omega = fill_coef(nmat(overt_util, overt_util, function(r, c) paste0("prior_b_", r, "_", c)), draws_colnames, d = d),
      prior_Delta = fill_coef(nmat(latent, overt_util, function(r, c) paste0("prior_b_", r, "_", c)), draws_colnames, d = d),
      prior_Psi_Y = fill_cov(nmat(latent_util, "1", function(r, c) paste0("prior_sigma_", r)), nmat(latent_util, latent_util, function(r, c) paste0("prior_rescor__", r, "__", c)), draws_colnames, d = d),
      prior_Psi_X = fill_cov(nmat(overt_util, "1", function(r, c) paste0("prior_sigma_", r)), nmat(overt_util, overt_util, function(r, c) paste0("prior_rescor__", r, "__", c)), draws_colnames, d = d),
      prior_Psi_Y_X = prior_Psi_Y_X,
      prior_Psi_X_Y <- aperm(prior_Psi_Y_X, c(2, 1, 3)), # NOTE: Currently doesn't output anything for some reason
      prior_alpha_X = fill_coef(nmat(overt_util, "1", function(r, c) paste0("prior_b_", r, "_Intercept")), draws_colnames, d = d) |> aperm(c(2, 1, 3)),
      prior_alpha_Y = fill_coef(nmat(latent_util, "1", function(r, c) paste0("prior_b_", r, "_Intercept")), draws_colnames, d = d) |> aperm(c(2, 1, 3))
    )
  }
  if(prior == "only") return(output_prior)
  else if(prior == "no") return(output_posterior)
  else return(c(output_posterior, output_prior))
}

# Function for computing the Bayesian factor score reliability (BFSR) index
# arg 'brms_fit': A brmsfit object containing the Bayesian measurement model
# arg 'overt': A positive length character vector containing the names of the overt variables
# arg 'latent': A positive length character vector containing the names of the latent variables
# arg 'method': A length one character vector containing the name of the method for generating the sets of plausible values. Defaults to (and currently only supports) the method commonly attributed to Thompson and Thurstone
# arg 'sample_prior': A length one character vector specifying whether to compute the BSFR prior, the BSFR posteior, or both. This presupposes that the brmsfit object was fit using a compatible 'sample_prior' argument
# arg 'ndraws': An optional integer scalar specifying the number of draws used to generate the sets of plausible values. The specified value must be compatible with the actual number of post-warmup draws in the brmsfit object
# Output: A named l by d array of factor score reliability indices, where l is the number of latent variables and d is the number of draws
fsr_bayes <- function(brms_fit = NULL, draws = NULL, overt = NULL, latent = NULL, method = NULL, prior = NULL,  ndraws = NULL){
  require(tidyverse)
  require(einsum)
  require(brms)
  require(posterior)
  
  # Input validation
  if(!is.null(draws)){
    if(!is.matrix(draws) & !is.data.frame(draws)) stop("Argument 'draws' must be a matrix or data.frame object.")
    if(!dim(draws)[1] > 0) stop("Argument 'draws' must have a positive number of rows.")
    if(!dim(draws)[2] > 0) stop("Argument 'draws' must have a positive number of columns.")
  }
  else if(!is.null(brms_fit)){
    if(!brms::is.brmsfit(brms_fit)) stop("Argument 'brms_fit' must be a brmsfit object.")
    draws <- brms::as_draws_df(brms_fit)
  }
  else{
    stop("Arguments 'brms_fit' and 'draws' cannot both be NULL.")
  }
  
  if(is.null(overt)) stop("Argument 'overt' cannot be NULL")
  if(!is.character(overt) | !is.vector(overt)) stop("Argument 'overt' must be a character vector")
  if(!length(overt) > 0) stop("Argument 'overt' must have a positive length")
  
  if(is.null(latent)) stop("Argument 'latent' cannot be NULL")
  if(!is.character(latent) | !is.vector(latent)) stop("Argument 'latent' must be a character vector")
  if(!length(latent) > 0) stop("Argument 'latent' must have a positive length")
  
  method <- match.arg(arg = method, choices = c("regression"))
  
  prior <- match.arg(arg = prior, choices = c("no", "yes", "only"))
  
  
  # Organize input
  overt <- sort(overt)
  latent <- sort(latent)
  draws_colnames <- sort(colnames(draws))
  draws <- draws[,draws_colnames]
  
  # Input validation of the optional arguments
  if(!is.null(ndraws)){
    if(length(ndraws) != 1) stop("Argument 'ndraws' must have a length of one")
    if(!is.numeric(ndraws)) stop("Argument 'ndraws' must be numeric")
    if(ndraws > nrow(draws)) stop("Argument 'ndraws' cannot be greater than the number of draws specified by the objects in the 'brms_fit' or 'draws' argument.")
    draws <- draws[sample(1:nrow(draws), ndraws, replace = FALSE),]
  }
  
  # Get dimensions
  o <- length(overt)
  l <- length(latent)
  d <- nrow(draws)
  
  
  # Get parameter matrices
  parameter_matrices <- get_parameter_matrices(draws = draws, latent = latent, overt = overt, prior = prior, ndraws = ndraws)
  I_ood <- array(diag(1, o, o), dim = c(o, o, d))
  I_lld <- array(diag(1, l, l), dim = c(l, l, d))
  one1ld <- array(matrix(1, 1, l), dim = c(1, l, d))
  
  if(prior != "only"){
    Lambda <- parameter_matrices$Lambda
    Beta <- parameter_matrices$Beta
    Omega <- parameter_matrices$Omega
    Delta <- parameter_matrices$Delta
    Psi_Y <- parameter_matrices$Psi_Y
    Psi_X <- parameter_matrices$Psi_X
    Psi_Y_X <- parameter_matrices$Psi_Y_X
    Psi_X_Y <- aperm(parameter_matrices$Psi_Y_X, c(2, 1, 3))
    
    # Compute additional parameters
    Omega_star <- invd(I_ood %-d% Omega); dimnames(Omega_star)[[1]] <- overt; dimnames(Omega_star)[[2]] <- overt
    Beta_star <- invd(I_lld %-d% Beta); dimnames(Beta_star)[[1]] <- latent; dimnames(Beta_star)[[2]] <- latent
    G_star <- invd(I_lld %-d% Beta %-d% (Delta %*d% Omega_star %*d% Lambda)); dimnames(G_star)[[1]] <- latent; dimnames(G_star)[[2]] <- latent
    H_star <- invd(I_ood %-d% Omega %-d% (Lambda %*d% Beta_star %*d% Delta)); dimnames(H_star)[[1]] <- overt; dimnames(H_star)[[2]] <- overt
    G <- Omega_star %*d% Lambda %*d% G_star; dimnames(G)[[1]] <- overt; dimnames(G)[[2]] <- latent
    H <- Beta_star %*d% Delta %*d% H_star; dimnames(H)[[1]] <- latent; dimnames(H)[[2]] <- overt
    Sigma_Y <- G_star %*d% Psi_Y %*d% td(G_star) + H %*d% Psi_X %*d% td(H) + G_star %*d% Psi_Y_X %*d% td(H) + H %*d% Psi_X_Y %*d% td(G_star); dimnames(Sigma_Y)[[1]] <- latent; dimnames(Sigma_Y)[[2]] <- latent
    Sigma_X <- G %*d% Psi_Y %*d% td(G) + H_star %*d% Psi_X %*d% td(H_star)  + G %*d% Psi_Y_X %*d% td(H_star) + H_star %*d% Psi_X_Y %*d% td(G); dimnames(Sigma_X)[[1]] <- overt; dimnames(Sigma_X)[[2]] <- overt
    Sigma_Y_X <- G_star %*d% Psi_Y %*d% td(G) + H %*d% Psi_X %*d% td(H_star) + G_star %*d% Psi_Y_X %*d% td(H_star) + H %*d% Psi_X_Y %*d% td(G); dimnames(Sigma_Y_X)[[1]] <- latent; dimnames(Sigma_Y_X)[[2]] <- overt
    Sigma_X_Y <- G %*d% Psi_Y %*d% td(G_star) + H_star %*d% Psi_X %*d% td(H) + G %*d% Psi_Y_X %*d% td(H) + H_star %*d% Psi_X_Y %*d% td(G_star); dimnames(Sigma_X_Y)[[1]] <- overt; dimnames(Sigma_X_Y)[[2]] <- latent
    
    W_TT <- Sigma_Y_X %*d% invd(Sigma_X); dimnames(W_TT)[[1]] <- latent; dimnames(W_TT)[[2]] <- overt
    W_TT_star <- W_TT %*d% Sigma_X_Y; dimnames(W_TT_star)[[1]] <- latent; dimnames(W_TT_star)[[2]] <- latent
    bfsr <- one1ld %*d% (W_TT_star %od% I_lld) %*d% invd(Sigma_Y %od% I_lld); dimnames(bfsr)[[2]] <- paste0("fsr_", latent)
    bfsr <- posterior::as_draws_df(aperm(bfsr, c(3, 1, 2))); colnames(bfsr)[1:l] <- paste0("fsr_", latent)
    if(posterior::is_draws(draws)){
      bfsr$.chain <- draws$.chain
      bfsr$.iteration <- draws$.iteration
      bfsr$.draw <- draws$.draw
    }
  }
  if(prior != "no"){
    prior_Lambda <- parameter_matrices$prior_Lambda
    prior_Beta <- parameter_matrices$prior_Beta
    prior_Omega <- parameter_matrices$prior_Omega
    prior_Delta <- parameter_matrices$prior_Delta
    prior_Psi_Y <- parameter_matrices$prior_Psi_Y
    prior_Psi_X <- parameter_matrices$prior_Psi_X
    prior_Psi_Y_X <- parameter_matrices$prior_Psi_Y_X
    prior_Psi_X_Y <- aperm(parameter_matrices$prior_Psi_Y_X, c(2, 1, 3))
    
    # Compute additional parameters
    prior_Omega_star <- invd(I_ood %-d% prior_Omega); dimnames(prior_Omega_star)[[1]] <- overt; dimnames(prior_Omega_star)[[2]] <- overt
    prior_Beta_star <- invd(I_lld %-d% prior_Beta); dimnames(prior_Beta_star)[[1]] <- latent; dimnames(prior_Beta_star)[[2]] <- latent
    prior_G_star <- invd(I_lld %-d% prior_Beta %-d% (prior_Delta %*d% prior_Omega_star %*d% prior_Lambda)); dimnames(prior_G_star)[[1]] <- latent; dimnames(prior_G_star)[[2]] <- latent
    prior_H_star <- invd(I_ood %-d% prior_Omega %-d% (prior_Lambda %*d% prior_Beta_star %*d% prior_Delta)); dimnames(prior_H_star)[[1]] <- overt; dimnames(prior_H_star)[[2]] <- overt
    prior_G <- prior_Omega_star %*d% prior_Lambda %*d% prior_G_star; dimnames(prior_G)[[1]] <- overt; dimnames(prior_G)[[2]] <- latent
    prior_H <- prior_Beta_star %*d% prior_Delta %*d% prior_H_star; dimnames(prior_H)[[1]] <- latent; dimnames(prior_H)[[2]] <- overt
    prior_Sigma_Y <- prior_G_star %*d% prior_Psi_Y %*d% td(prior_G_star) + prior_H %*d% prior_Psi_X %*d% td(prior_H) + prior_G_star %*d% prior_Psi_Y_X %*d% td(prior_H) + prior_H %*d% prior_Psi_X_Y %*d% td(prior_G_star); dimnames(prior_Sigma_Y)[[1]] <- latent; dimnames(prior_Sigma_Y)[[2]] <- latent
    prior_Sigma_X <- prior_G %*d% prior_Psi_Y %*d% td(prior_G) + prior_H_star %*d% prior_Psi_X %*d% td(prior_H_star)  + prior_G %*d% prior_Psi_Y_X %*d% td(prior_H_star) + prior_H_star %*d% prior_Psi_X_Y %*d% td(prior_G); dimnames(prior_Sigma_X)[[1]] <- overt; dimnames(prior_Sigma_X)[[2]] <- overt
    prior_Sigma_Y_X <- prior_G_star %*d% prior_Psi_Y %*d% td(prior_G) + prior_H %*d% prior_Psi_X %*d% td(prior_H_star) + prior_G_star %*d% prior_Psi_Y_X %*d% td(prior_H_star) + prior_H %*d% prior_Psi_X_Y %*d% td(prior_G); dimnames(prior_Sigma_Y_X)[[1]] <- latent; dimnames(prior_Sigma_Y_X)[[2]] <- overt
    prior_Sigma_X_Y <- prior_G %*d% prior_Psi_Y %*d% td(prior_G_star) + prior_H_star %*d% prior_Psi_X %*d% td(prior_H) + prior_G %*d% prior_Psi_Y_X %*d% td(prior_H) + prior_H_star %*d% prior_Psi_X_Y %*d% td(prior_G_star); dimnames(prior_Sigma_X_Y)[[1]] <- overt; dimnames(prior_Sigma_X_Y)[[2]] <- latent
    
    prior_W_TT <- prior_Sigma_Y_X %*d% invd(prior_Sigma_X); dimnames(prior_W_TT)[[1]] <- latent; dimnames(prior_W_TT)[[2]] <- overt
    prior_W_TT_star <- prior_W_TT %*d% prior_Sigma_X_Y; dimnames(prior_W_TT_star)[[1]] <- latent; dimnames(prior_W_TT_star)[[2]] <- latent
    prior_bfsr <- one1ld %*d% (prior_W_TT_star %od% I_lld) %*d% invd(prior_Sigma_Y %od% I_lld); dimnames(prior_bfsr)[[2]] <- paste0("prior_", "fsr_", latent)
    prior_bfsr <- posterior::as_draws_df(aperm(prior_bfsr, c(3, 1, 2))); colnames(prior_bfsr)[1:l] <- paste0("prior_", "fsr_", latent)
    if(posterior::is_draws(draws)){
      prior_bfsr$.chain <- draws$.chain
      prior_bfsr$.iteration <- draws$.iteration
      prior_bfsr$.draw <- draws$.draw
    }
  }
  
  if(prior == "only") return(prior_bfsr)
  else if(prior == "no") return(bfsr)
  else return(posterior::bind_draws(bfsr, prior_bfsr, along = "variable"))
}





# Load empirical data -----------------------------------------------------

# Load HEXACO items (see documentation: https://ipip.ori.org/newhexaco_pi_key.htm)
data(HEXACO, package = "IPV")
df <- HEXACO
rm(HEXACO)
# dim(df)
# |N| = 20,365   |K| = 244

# colSums(is.na(df))
# A few missing values for some items


# Save the data
# write.csv(df, "Demonstration/Data/HEXACO.csv", row.names = FALSE)


# Limit data to variables of interest
df <- df %>%
  dplyr::select(-c(V1, V2, country, elapse)) %>%
  na.omit()
# dim(df)
# |N| = 20,326   |K| = 240





# Sample data -------------------------------------------------------------

# Subset data into training, test, and unused samples
df$org_row_nr <- row.names(df)
df$subset <- sample(
  c(rep("Sampled", sample_size), rep("Unused", nrow(df) - sample_size)),
  nrow(df),
  replace = FALSE
)
df_sample <- df %>% dplyr::filter(subset == "Sampled") %>% dplyr::select(-c(subset))
row.names(df_sample) <- NULL
dim(df_sample)
# |N| = 1000  |K| = 241
write.csv(df_sample, "Demonstration/Data/df_sample.csv", row.names = FALSE)





# Load sampled data -------------------------------------------------------


# Specify list of latent phenomena names
trait <- "honesty_humility"
facets <- c("sincerity", "fairness", "greed_avoidance", "modesty")


# Load sampled data
df_sample <- readr::read_csv("Demonstration/Data/df_sample.csv") %>% dplyr::select(-c(org_row_nr)) %>% as.data.frame()


# Specify list of indicator names
df_colnames <- colnames(df_sample)
honesty_humility_indicators <- df_colnames[grepl("H_", df_colnames)]
sincerity_indicators <- honesty_humility_indicators[grepl("H_Sinc", honesty_humility_indicators)]
fairness_indicators <- honesty_humility_indicators[grepl("H_Fair", honesty_humility_indicators)]
greed_avoidance_indicators <- honesty_humility_indicators[grepl("H_Gree", honesty_humility_indicators)]
modesty_indicators <- honesty_humility_indicators[grepl("H_Mode", honesty_humility_indicators)]


# Limit data to the variables of interest
df_sample <- df_sample %>% dplyr::select(dplyr::any_of(honesty_humility_indicators)) %>% as.data.frame()


# Standardize data using frequentist formula
df_sample_freq <- apply(df_sample, 2, datawizard::standardize) %>% as.data.frame()


# Standardize data using likelihood formula
df_sample_bayes <- apply(df_sample, 2, standardize2) %>% as.data.frame()


# Get specifications
overt <- sort(honesty_humility_indicators)
latent <- sort(c(trait, facets))
o <- length(overt); l <- length(latent); d <- posterior_draws
n <- nrow(df_sample)


# Facet indicator to constrain for identification
# sample(facets, 1, replace = FALSE)
# "fairness"

# Sincerity indicator to constrain for identification
# sample(sincerity_indicators, 1, replace = FALSE)
# "H_Sinc_Sinc5"

# Fairness indicator to constrain for identification
# sample(fairness_indicators, 1, replace = FALSE)
# "H_Fair_Fair7"

# Greed-avoidance indicator to constrain for identification
# sample(greed_avoidance_indicators, 1, replace = FALSE)
# "H_Gree_Gree6"

# Modesty indicator to constrain for identification
# sample(modesty_indicators, 1, replace = FALSE)
# "H_Mode_Mode1"





# Bayesian CFA ------------------------------------------------------------


# Introduce latent variables to impute
df_sample_bayes$honesty_humility <- as.numeric(NA)
df_sample_bayes$sincerity <- as.numeric(NA)
df_sample_bayes$fairness <- as.numeric(NA)
df_sample_bayes$greed_avoidance <- as.numeric(NA)
df_sample_bayes$modesty <- as.numeric(NA)


# Specify brms-based model for Bayesian confirmatory factor analysis
bglm_honesty_humility_model <- brms::bf(honesty_humility | mi() ~ 1) +
  brms::bf(sincerity | mi() ~ 1 + mi(honesty_humility)) +
  brms::bf(fairness | mi() ~ 1 + mi(honesty_humility)) +
  brms::bf(greed_avoidance | mi() ~ 1 + mi(honesty_humility)) +
  brms::bf(modesty | mi() ~ 1 + mi(honesty_humility)) +
  
  brms::bf(H_Sinc_Sinc1 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc2 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc3 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc4 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc5 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc6 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc7 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc8 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc9 ~ 1 + mi(sincerity)) +
  brms::bf(H_Sinc_Sinc10 ~ 1 + mi(sincerity)) +
  
  brms::bf(H_Fair_Fair1 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair2 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair3 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair4 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair5 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair6 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair7 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair8 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair9 ~ 1 + mi(fairness)) +
  brms::bf(H_Fair_Fair10 ~ 1 + mi(fairness)) +
  
  brms::bf(H_Gree_Gree1 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree2 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree3 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree4 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree5 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree6 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree7 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree8 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree9 ~ 1 + mi(greed_avoidance)) +
  brms::bf(H_Gree_Gree10 ~ 1 + mi(greed_avoidance)) +
  
  brms::bf(H_Mode_Mode1 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode2 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode3 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode4 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode5 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode6 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode7 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode8 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode9 ~ 1 + mi(modesty)) +
  brms::bf(H_Mode_Mode10 ~ 1 + mi(modesty)) +
  
  brms::set_rescor(rescor = FALSE)


# Specify priors
# brms::get_prior(bglm_honesty_humility_model, data = df_sample_bayes, family = gaussian()) # NOTE: Provides a helpful starting point
# 1/critical.value(0.001)
# 0.304
# get_gamma_params_from_mode_sd(1, .304) # NOTE: Used for dispersion parameter
#    a    b 
# 12.7 11.7
# get_gamma_params_from_mode_sd(0, .304) # NOTE: Used for dispersion parameter
#    a    b 
# 1.00 3.29
bglm_honesty_humility_prior <- c(
  # Specify intercept priors
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc1"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc2"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc3"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc4"),
  brms::set_prior("constant(0)", class = "Intercept", resp = "HSincSinc5"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc6"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc7"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc8"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc9"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HSincSinc10"),
  
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair1"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair2"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair3"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair4"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair5"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair6"),
  brms::set_prior("constant(0)", class = "Intercept", resp = "HFairFair7"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair8"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair9"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HFairFair10"),
  
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree1"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree2"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree3"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree4"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree5"),
  brms::set_prior("constant(0)", class = "Intercept", resp = "HGreeGree6"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree7"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree8"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree9"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HGreeGree10"),
  
  brms::set_prior("constant(0)", class = "Intercept", resp = "HModeMode1"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode2"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode3"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode4"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode5"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode6"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode7"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode8"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode9"),
  brms::set_prior("normal(0, 1)", class = "Intercept", resp = "HModeMode10"),
  
  brms::set_prior("constant(0)", class = "Intercept", resp = "sincerity"),
  brms::set_prior("constant(0)", class = "Intercept", resp = "fairness"),
  brms::set_prior("constant(0)", class = "Intercept", resp = "greedavoidance"),
  brms::set_prior("constant(0)", class = "Intercept", resp = "modesty"),
  brms::set_prior("constant(0)", class = "Intercept", resp = "honestyhumility"),
  
  # Specify coefficient priors
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc1"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc2"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc3"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc4"),
  brms::set_prior("constant(1)", class = "b", coef = "misincerity", resp = "HSincSinc5"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc6"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc7"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc8"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc9"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "misincerity", resp = "HSincSinc10"),
  
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair1"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair2"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair3"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair4"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair5"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair6"),
  brms::set_prior("constant(1)", class = "b", coef = "mifairness", resp = "HFairFair7"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair8"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair9"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mifairness", resp = "HFairFair10"),
  
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree1"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree2"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree3"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree4"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree5"),
  brms::set_prior("constant(1)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree6"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree7"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree8"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree9"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "migreed_avoidance", resp = "HGreeGree10"),
  
  brms::set_prior("constant(1)", class = "b", coef = "mimodesty", resp = "HModeMode1"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode2"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode3"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode4"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode5"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode6"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode7"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode8"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode9"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mimodesty", resp = "HModeMode10"),
  
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mihonesty_humility", resp = "sincerity"),
  brms::set_prior("constant(1)", class = "b", coef = "mihonesty_humility", resp = "fairness"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mihonesty_humility", resp = "greedavoidance"),
  brms::set_prior("normal(0, 0.304)", class = "b", coef = "mihonesty_humility", resp = "modesty"),
  
  # Specify sigma priors
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc1", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc2", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc3", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc4", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc5", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc6", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc7", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc8", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc9", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HSincSinc10", lb = 0),
  
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair1", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair2", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair3", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair4", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair5", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair6", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair7", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair8", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair9", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HFairFair10", lb = 0),
  
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree1", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree2", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree3", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree4", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree5", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree6", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree7", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree8", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree9", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HGreeGree10", lb = 0),
  
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode1", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode2", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode3", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode4", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode5", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode6", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode7", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode8", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode9", lb = 0),
  brms::set_prior("gamma(12.7, 11.7)", class = "sigma", resp = "HModeMode10", lb = 0),
  
  brms::set_prior("gamma(1.00, 3.29)", class = "sigma", resp = "sincerity", lb = 0),
  brms::set_prior("gamma(1.00, 3.29)", class = "sigma", resp = "fairness", lb = 0),
  brms::set_prior("gamma(1.00, 3.29)", class = "sigma", resp = "greedavoidance", lb = 0),
  brms::set_prior("gamma(1.00, 3.29)", class = "sigma", resp = "modesty", lb = 0),
  brms::set_prior("gamma(1.00, 3.29)", class = "sigma", resp = "honestyhumility", lb = 0)
)

# Fit model
t0 <- Sys.time()
bglm_honesty_humility_fit <- brms::brm(
  bglm_honesty_humility_model,
  prior = bglm_honesty_humility_prior,
  family = gaussian(),
  data = df_sample_bayes,
  sample_prior = "yes",
  seed = SEED,
  cores = cores,
  chains = chains,
  iter = total_draws,
  warmup = warmup_draws,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file = "Demonstration/Models/bglm_honesty_humility_fit",
  file_refit = "always"
)
t1 <- Sys.time()
t1 - t0
# NOTE: Computation took approximately 12.4 hours with 15 cores


# Load model (alternative)
# bglm_honesty_humility_fit <- readRDS("Demonstration/Models/bglm_honesty_humility_fit.rds")


# Check priors
bglm_honesty_humility_fit %>% brms::prior_summary()


# Check diagnostics
bglm_honesty_humility_fit_diagnostics <- posterior::summarize_draws(
  bglm_honesty_humility_fit %>% brms::as_draws(),
  posterior::default_convergence_measures(),
  posterior::default_mcse_measures(),
  .cores = cores
) %>% as.data.frame() %>% na.omit()
max(bglm_honesty_humility_fit_diagnostics$rhat) # The largest rhat is 1.01
min(bglm_honesty_humility_fit_diagnostics$ess_bulk) # The smallest ess_bulk is 960
min(bglm_honesty_humility_fit_diagnostics$ess_tail) # The smallest ess_tail is 1412
max(bglm_honesty_humility_fit_diagnostics$mcse_mean) # The largest mcse is 7.2


# Assess Traceplots of the posterior
fig_width <- 16
fig_height <- fig_width/golden_ratio
bglm_honesty_humility_fit_intercept_traceplot <- bglm_honesty_humility_fit %>%
  bayesplot::mcmc_trace(regex_pars = ("^Intercept_"), np = nuts_params(bglm_honesty_humility_fit)) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 3000, 500), labels = ~ .x + warmup_draws + 1) +
  scale_y_continuous(breaks = seq(-10, 10, 0.05)) +
  scale_color_manual(values = traceplot_color_chart) +
  labs(
    x = "Draw",
    y = "Estimate",
    title = latex2exp::TeX("\\textbf{FIGURE}: Traceplots of estimated intercepts ($\\hat{\\alpha_{Y}}$, $\\hat{\\alpha_{X}}$)"),
    caption = latex2exp::TeX("\\textbf{NOTE}: Warmup draws per chain = 4,000. Post-warmup draws per chain = 2,667. Total post-warmup draws = 40,005."),
    color = "Chain"
  ) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()); bglm_honesty_humility_fit_intercept_traceplot
ggsave("Demonstration/Plots/bglm_honesty_humility_fit_intercept_traceplot.png", plot = bglm_honesty_humility_fit_intercept_traceplot, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi)
ggsave("Demonstration/Plots/bglm_honesty_humility_fit_intercept_traceplot.pdf", plot = bglm_honesty_humility_fit_intercept_traceplot, width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
bglm_honesty_humility_fit_bsp_traceplot <- bglm_honesty_humility_fit %>%
  bayesplot::mcmc_trace(regex_pars = ("^bsp_"), np = nuts_params(bglm_honesty_humility_fit)) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 3000, 500), labels = ~ .x + warmup_draws + 1) +
  scale_y_continuous(breaks = seq(-10, 10, 0.1)) +
  scale_color_manual(values = traceplot_color_chart) +
  labs(
    x = "Draw",
    y = "Estimate",
    title = latex2exp::TeX("\\textbf{FIGURE}: Traceplots of estimated factor loadings ($\\hat{\\lambda}$, $\\hat{\\beta}$)"),
    caption = latex2exp::TeX("\\textbf{NOTE}: Warmup draws per chain = 4,000. Post-warmup draws per chain = 2,667. Total post-warmup draws = 40,005."),
    color = "Chain"
  ) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()); bglm_honesty_humility_fit_bsp_traceplot
ggsave("Demonstration/Plots/bglm_honesty_humility_fit_bsp_traceplot.png", plot = bglm_honesty_humility_fit_bsp_traceplot, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi)
ggsave("Demonstration/Plots/bglm_honesty_humility_fit_bsp_traceplot.pdf", plot = bglm_honesty_humility_fit_bsp_traceplot, width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
bglm_honesty_humility_fit_sigma_traceplot <- bglm_honesty_humility_fit %>%
  bayesplot::mcmc_trace(regex_pars = ("^sigma_"), np = nuts_params(bglm_honesty_humility_fit)) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 3000, 500), labels = ~ .x + warmup_draws + 1) +
  scale_y_continuous(breaks = seq(-10, 10, 0.05)) +
  scale_color_manual(values = traceplot_color_chart) +
  labs(
    x = "Draw",
    y = "Estimate",
    title = latex2exp::TeX("\\textbf{FIGURE}: Traceplots of estimated dispersions ($\\hat{\\psi_{Y}}$, $\\hat{\\psi_{X}}$)"),
    caption = latex2exp::TeX("\\textbf{NOTE}: Warmup draws per chain = 4,000. Post-warmup draws per chain = 2,667. Total post-warmup draws = 40,005."),
    color = "Chain"
  ) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()); bglm_honesty_humility_fit_sigma_traceplot
ggsave("Demonstration/Plots/bglm_honesty_humility_fit_sigma_traceplot.png", plot = bglm_honesty_humility_fit_sigma_traceplot, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi)
ggsave("Demonstration/Plots/bglm_honesty_humility_fit_sigma_traceplot.pdf", plot = bglm_honesty_humility_fit_sigma_traceplot, width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
# NOTE: Traceplots looks good, not reported for brevity.





# Frequentist CFA ---------------------------------------------------------


# Specify lavaan-based model for frequentist confirmatory factor analysis
glm_honesty_humility_model <- paste0('
  honesty_humility =~ NA*sincerity + 1*fairness + NA*greed_avoidance + NA*modesty',
                                     paste0('\n sincerity =~ 1*H_Sinc_Sinc5 + NA*', paste0(sincerity_indicators[sincerity_indicators != 'H_Sinc_Sinc5'], collapse = " + NA*")),
                                     paste0('\n fairness =~ 1*H_Fair_Fair7 + NA*', paste0(fairness_indicators[fairness_indicators != 'H_Fair_Fair7'], collapse = " + NA*")),
                                     paste0('\n greed_avoidance =~ 1*H_Gree_Gree6 + NA*', paste0(greed_avoidance_indicators[greed_avoidance_indicators != 'H_Gree_Gree6'], collapse = " + NA*")),
                                     paste0('\n modesty =~ 1*H_Mode_Mode1 + NA*', paste0(modesty_indicators[modesty_indicators != 'H_Mode_Mode1'], collapse = " + NA*")),
                                     '\n honesty_humility ~ 0*1',
                                     '\n sincerity ~ 0*1',
                                     '\n fairness ~ 0*1',
                                     '\n greed_avoidance ~ 0*1',
                                     '\n modesty ~ 0*1',
                                     '\n H_Sinc_Sinc5 ~ 0*1',
                                     '\n H_Fair_Fair7 ~ 0*1',
                                     '\n H_Gree_Gree6 ~ 0*1',
                                     '\n H_Mode_Mode1 ~ 0*1',
                                     paste0('\n ', sincerity_indicators[sincerity_indicators != 'H_Sinc_Sinc5'], ' ~ NA*1', collapse = ""),
                                     paste0('\n ', fairness_indicators[fairness_indicators != 'H_Fair_Fair7'], ' ~ NA*1', collapse = ""),
                                     paste0('\n ', greed_avoidance_indicators[greed_avoidance_indicators != 'H_Gree_Gree6'], ' ~ NA*1', collapse = ""),
                                     paste0('\n ', modesty_indicators[modesty_indicators != 'H_Mode_Mode1'], ' ~ NA*1', collapse = ""),
                                     '\n honesty_humility ~~ NA*honesty_humility',
                                     '\n sincerity ~~ NA*sincerity',
                                     '\n fairness ~~ NA*fairness',
                                     '\n greed_avoidance ~~ NA*greed_avoidance',
                                     '\n modesty ~~ NA*modesty',
                                     paste0('\n ', honesty_humility_indicators, ' ~~ NA*', honesty_humility_indicators, collapse = "")
)


# Fit model to validate specification, sample size, and number of parameters
# cfa_fit <- lavaan::cfa(
#   glm_honesty_humility_model,
#   data = df_sample_freq,
#   std.lv = FALSE
# )
# summary(cfa_fit, standardized = TRUE)
# lavaan::lavInspect(cfa_fit, what = "ntotal") # NOTE: n = 2,000 (As expected)
# lavaan::lavInspect(cfa_fit, what = "npar") # NOTE: p = 120 (As expected)



# Fit bootstrapped-based model
t0 <- Sys.time()
glm_honesty_humility_fit <- foreach::foreach(i = 1:posterior_draws, .combine = "rbind") %dopar% {
  require(lavaan)
  while(TRUE){
    df_sample_freq_i <- df_sample_freq[sample(1:nrow(df_sample_freq), nrow(df_sample_freq), replace = TRUE),]
    row.names(df_sample_freq_i) <- NULL
    tryCatch(
      expr = {
        cfa_fit_i <- lavaan::cfa(
          glm_honesty_humility_model,
          data = df_sample_freq_i,
          std.lv = FALSE
        )
      }
    )
    if(!is.null(cfa_fit_i) & lavaan::lavInspect(cfa_fit_i, what = "converged") & lavInspect(cfa_fit_i, "post.check")) break
  }
  cfa_fit_i_estimates <- lavaan::lavInspect(cfa_fit_i, what = "est")
  cfa_fit_i_psi_sqrt <- sqrt(cfa_fit_i_estimates$psi)
  cfa_fit_i_theta_sqrt <- sqrt(cfa_fit_i_estimates$theta)
  
  cfa_fit_i_fsr <- lavaan::lavInspect(cfa_fit_i, what = "fs.reliability")
  return(data.frame(
    # Latent Intercepts
    b_honestyhumility_Intercept = cfa_fit_i_estimates$alpha["honesty_humility",],
    b_sincerity_Intercept = cfa_fit_i_estimates$alpha["sincerity",],
    b_fairness_Intercept = cfa_fit_i_estimates$alpha["fairness",],
    b_greedavoidance_Intercept = cfa_fit_i_estimates$alpha["greed_avoidance",],
    b_modesty_Intercept = cfa_fit_i_estimates$alpha["modesty",],
    
    # Indicator Intercepts
    b_HSincSinc1_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc1",],
    b_HSincSinc2_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc2",],
    b_HSincSinc3_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc3",],
    b_HSincSinc4_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc4",],
    b_HSincSinc5_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc5",],
    b_HSincSinc6_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc6",],
    b_HSincSinc7_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc7",],
    b_HSincSinc8_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc8",],
    b_HSincSinc9_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc9",],
    b_HSincSinc10_Intercept = cfa_fit_i_estimates$nu["H_Sinc_Sinc10",],
    
    b_HFairFair1_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair1",],
    b_HFairFair2_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair2",],
    b_HFairFair3_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair3",],
    b_HFairFair4_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair4",],
    b_HFairFair5_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair5",],
    b_HFairFair6_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair6",],
    b_HFairFair7_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair7",],
    b_HFairFair8_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair8",],
    b_HFairFair9_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair9",],
    b_HFairFair10_Intercept = cfa_fit_i_estimates$nu["H_Fair_Fair10",],
    
    b_HGreeGree1_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree1",],
    b_HGreeGree2_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree2",],
    b_HGreeGree3_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree3",],
    b_HGreeGree4_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree4",],
    b_HGreeGree5_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree5",],
    b_HGreeGree6_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree6",],
    b_HGreeGree7_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree7",],
    b_HGreeGree8_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree8",],
    b_HGreeGree9_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree9",],
    b_HGreeGree10_Intercept = cfa_fit_i_estimates$nu["H_Gree_Gree10",],
    
    b_HModeMode1_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode1",],
    b_HModeMode2_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode2",],
    b_HModeMode3_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode3",],
    b_HModeMode4_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode4",],
    b_HModeMode5_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode5",],
    b_HModeMode6_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode6",],
    b_HModeMode7_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode7",],
    b_HModeMode8_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode8",],
    b_HModeMode9_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode9",],
    b_HModeMode10_Intercept = cfa_fit_i_estimates$nu["H_Mode_Mode10",],
    
    # Beta
    bsp_sincerity_mihonesty_humility = cfa_fit_i_estimates$beta["sincerity","honesty_humility"],
    bsp_fairness_mihonesty_humility = cfa_fit_i_estimates$beta["fairness","honesty_humility"],
    bsp_greedavoidance_mihonesty_humility = cfa_fit_i_estimates$beta["greed_avoidance","honesty_humility"],
    bsp_modesty_mihonesty_humility = cfa_fit_i_estimates$beta["modesty","honesty_humility"],
    
    # Lambda
    bsp_HSincSinc1_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc1","sincerity"],
    bsp_HSincSinc2_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc2","sincerity"],
    bsp_HSincSinc3_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc3","sincerity"],
    bsp_HSincSinc4_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc4","sincerity"],
    bsp_HSincSinc5_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc5","sincerity"],
    bsp_HSincSinc6_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc6","sincerity"],
    bsp_HSincSinc7_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc7","sincerity"],
    bsp_HSincSinc8_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc8","sincerity"],
    bsp_HSincSinc9_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc9","sincerity"],
    bsp_HSincSinc10_misincerity = cfa_fit_i_estimates$lambda["H_Sinc_Sinc10","sincerity"],
    
    bsp_HFairFair1_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair1","fairness"],
    bsp_HFairFair2_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair2","fairness"],
    bsp_HFairFair3_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair3","fairness"],
    bsp_HFairFair4_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair4","fairness"],
    bsp_HFairFair5_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair5","fairness"],
    bsp_HFairFair6_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair6","fairness"],
    bsp_HFairFair7_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair7","fairness"],
    bsp_HFairFair8_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair8","fairness"],
    bsp_HFairFair9_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair9","fairness"],
    bsp_HFairFair10_mifairness = cfa_fit_i_estimates$lambda["H_Fair_Fair10","fairness"],
    
    bsp_HGreeGree1_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree1","greed_avoidance"],
    bsp_HGreeGree2_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree2","greed_avoidance"],
    bsp_HGreeGree3_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree3","greed_avoidance"],
    bsp_HGreeGree4_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree4","greed_avoidance"],
    bsp_HGreeGree5_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree5","greed_avoidance"],
    bsp_HGreeGree6_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree6","greed_avoidance"],
    bsp_HGreeGree7_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree7","greed_avoidance"],
    bsp_HGreeGree8_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree8","greed_avoidance"],
    bsp_HGreeGree9_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree9","greed_avoidance"],
    bsp_HGreeGree10_migreed_avoidance = cfa_fit_i_estimates$lambda["H_Gree_Gree10","greed_avoidance"],
    
    bsp_HModeMode1_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode1","modesty"],
    bsp_HModeMode2_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode2","modesty"],
    bsp_HModeMode3_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode3","modesty"],
    bsp_HModeMode4_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode4","modesty"],
    bsp_HModeMode5_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode5","modesty"],
    bsp_HModeMode6_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode6","modesty"],
    bsp_HModeMode7_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode7","modesty"],
    bsp_HModeMode8_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode8","modesty"],
    bsp_HModeMode9_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode9","modesty"],
    bsp_HModeMode10_mimodesty = cfa_fit_i_estimates$lambda["H_Mode_Mode10","modesty"],
    
    
    # Latent dispersions
    sigma_honestyhumility = cfa_fit_i_psi_sqrt["honesty_humility", "honesty_humility"],
    sigma_sincerity = cfa_fit_i_psi_sqrt["sincerity", "sincerity"],
    sigma_greedavoidance = cfa_fit_i_psi_sqrt["fairness", "fairness"],
    sigma_fairness = cfa_fit_i_psi_sqrt["greed_avoidance", "greed_avoidance"],
    sigma_modesty = cfa_fit_i_psi_sqrt["modesty", "modesty"],
    
    
    # Indicator dispersions
    sigma_HSincSinc1 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc1", "H_Sinc_Sinc1"],
    sigma_HSincSinc2 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc2", "H_Sinc_Sinc2"],
    sigma_HSincSinc3 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc3", "H_Sinc_Sinc3"],
    sigma_HSincSinc4 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc4", "H_Sinc_Sinc4"],
    sigma_HSincSinc5 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc5", "H_Sinc_Sinc5"],
    sigma_HSincSinc6 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc6", "H_Sinc_Sinc6"],
    sigma_HSincSinc7 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc7", "H_Sinc_Sinc7"],
    sigma_HSincSinc8 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc8", "H_Sinc_Sinc8"],
    sigma_HSincSinc9 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc9", "H_Sinc_Sinc9"],
    sigma_HSincSinc10 = cfa_fit_i_theta_sqrt["H_Sinc_Sinc10", "H_Sinc_Sinc10"],
    
    sigma_HFairFair1 = cfa_fit_i_theta_sqrt["H_Fair_Fair1", "H_Fair_Fair1"],
    sigma_HFairFair2 = cfa_fit_i_theta_sqrt["H_Fair_Fair2", "H_Fair_Fair2"],
    sigma_HFairFair3 = cfa_fit_i_theta_sqrt["H_Fair_Fair3", "H_Fair_Fair3"],
    sigma_HFairFair4 = cfa_fit_i_theta_sqrt["H_Fair_Fair4", "H_Fair_Fair4"],
    sigma_HFairFair5 = cfa_fit_i_theta_sqrt["H_Fair_Fair5", "H_Fair_Fair5"],
    sigma_HFairFair6 = cfa_fit_i_theta_sqrt["H_Fair_Fair6", "H_Fair_Fair6"],
    sigma_HFairFair7 = cfa_fit_i_theta_sqrt["H_Fair_Fair7", "H_Fair_Fair7"],
    sigma_HFairFair8 = cfa_fit_i_theta_sqrt["H_Fair_Fair8", "H_Fair_Fair8"],
    sigma_HFairFair9 = cfa_fit_i_theta_sqrt["H_Fair_Fair9", "H_Fair_Fair9"],
    sigma_HFairFair10 = cfa_fit_i_theta_sqrt["H_Fair_Fair10", "H_Fair_Fair10"],
    
    sigma_HGreeGree1 = cfa_fit_i_theta_sqrt["H_Gree_Gree1", "H_Gree_Gree1"],
    sigma_HGreeGree2 = cfa_fit_i_theta_sqrt["H_Gree_Gree2", "H_Gree_Gree2"],
    sigma_HGreeGree3 = cfa_fit_i_theta_sqrt["H_Gree_Gree3", "H_Gree_Gree3"],
    sigma_HGreeGree4 = cfa_fit_i_theta_sqrt["H_Gree_Gree4", "H_Gree_Gree4"],
    sigma_HGreeGree5 = cfa_fit_i_theta_sqrt["H_Gree_Gree5", "H_Gree_Gree5"],
    sigma_HGreeGree6 = cfa_fit_i_theta_sqrt["H_Gree_Gree6", "H_Gree_Gree6"],
    sigma_HGreeGree7 = cfa_fit_i_theta_sqrt["H_Gree_Gree7", "H_Gree_Gree7"],
    sigma_HGreeGree8 = cfa_fit_i_theta_sqrt["H_Gree_Gree8", "H_Gree_Gree8"],
    sigma_HGreeGree9 = cfa_fit_i_theta_sqrt["H_Gree_Gree9", "H_Gree_Gree9"],
    sigma_HGreeGree10 = cfa_fit_i_theta_sqrt["H_Gree_Gree10", "H_Gree_Gree10"],
    
    sigma_HModeMode1 = cfa_fit_i_theta_sqrt["H_Mode_Mode1", "H_Mode_Mode1"],
    sigma_HModeMode2 = cfa_fit_i_theta_sqrt["H_Mode_Mode2", "H_Mode_Mode2"],
    sigma_HModeMode3 = cfa_fit_i_theta_sqrt["H_Mode_Mode3", "H_Mode_Mode3"],
    sigma_HModeMode4 = cfa_fit_i_theta_sqrt["H_Mode_Mode4", "H_Mode_Mode4"],
    sigma_HModeMode5 = cfa_fit_i_theta_sqrt["H_Mode_Mode5", "H_Mode_Mode5"],
    sigma_HModeMode6 = cfa_fit_i_theta_sqrt["H_Mode_Mode6", "H_Mode_Mode6"],
    sigma_HModeMode7 = cfa_fit_i_theta_sqrt["H_Mode_Mode7", "H_Mode_Mode7"],
    sigma_HModeMode8 = cfa_fit_i_theta_sqrt["H_Mode_Mode8", "H_Mode_Mode8"],
    sigma_HModeMode9 = cfa_fit_i_theta_sqrt["H_Mode_Mode9", "H_Mode_Mode9"],
    sigma_HModeMode10 = cfa_fit_i_theta_sqrt["H_Mode_Mode10", "H_Mode_Mode10"],
    
    
    # Factor score reliability
    fsr_honesty_humility = cfa_fit_i_fsr["honesty_humility"],
    fsr_sincerity = cfa_fit_i_fsr["sincerity"],
    fsr_fairness = cfa_fit_i_fsr["fairness"],
    fsr_greed_avoidance = cfa_fit_i_fsr["greed_avoidance"],
    fsr_modesty = cfa_fit_i_fsr["modesty"]
  ))
}
t1 <- Sys.time()
t1 - t0
# NOTE: Computation took approximately 37.3 mins with 15 cores


# Save results
write.csv(glm_honesty_humility_fit, "Demonstration/Results/glm_honesty_humility_fit.csv", row.names = FALSE)





# Load models -------------------------------------------------------------

# Load Bayesian model
bglm_honesty_humility_fit <- readRDS("Demonstration/Models/bglm_honesty_humility_fit.rds")

# Extract and sort posterior draws
bglm_honesty_humility_draws <- brms::as_draws_df(bglm_honesty_humility_fit)
bglm_honesty_humility_draws <- bglm_honesty_humility_draws[,sort(colnames(bglm_honesty_humility_draws))]


# Load frequentist results
glm_honesty_humility_fit <- readr::read_csv("Demonstration/Results/glm_honesty_humility_fit.csv") %>% as.data.frame()


# Apply finite-sample adjustment to the Frequentist FSR indices
glm_honesty_humility_fit$fsr_fairness_adj <- fsr_to_fsr_adj(glm_honesty_humility_fit$fsr_fairness, n = 2000, p = 120)
glm_honesty_humility_fit$fsr_greed_avoidance_adj <- fsr_to_fsr_adj(glm_honesty_humility_fit$fsr_greed_avoidance, n = 2000, p = 120)
glm_honesty_humility_fit$fsr_honesty_humility_adj <- fsr_to_fsr_adj(glm_honesty_humility_fit$fsr_honesty_humility, n = 2000, p = 120)
glm_honesty_humility_fit$fsr_modesty_adj <- fsr_to_fsr_adj(glm_honesty_humility_fit$fsr_modesty, n = 2000, p = 120)
glm_honesty_humility_fit$fsr_sincerity_adj <- fsr_to_fsr_adj(glm_honesty_humility_fit$fsr_sincerity, n = 2000, p = 120)





# Compute BFSR index ------------------------------------------------------


# Compute BFSR index
bglm_honesty_humility_fsr <- fsr_bayes(draws = bglm_honesty_humility_draws, overt = overt, latent = latent, prior = "yes")


# Check diagnostics
bglm_honesty_humility_fsr_diagnostics <- posterior::summarize_draws(
  bglm_honesty_humility_fsr %>% brms::as_draws(),
  posterior::default_convergence_measures(),
  posterior::default_mcse_measures(),
  .cores = cores
) %>% as.data.frame() %>% na.omit()
max(bglm_honesty_humility_fsr_diagnostics$rhat) # The largest rhat is 1.01
min(bglm_honesty_humility_fsr_diagnostics$ess_bulk) # The smallest ess_bulk is 1764
min(bglm_honesty_humility_fsr_diagnostics$ess_tail) # The smallest ess_tail is 3217
max(bglm_honesty_humility_fsr_diagnostics$mcse_mean) # The largest mcse is 0.000761


# Assess Traceplots of the posterior
fig_width <- 16
fig_height <- fig_width/golden_ratio
bglm_honesty_humility_fsr_traceplot <- bglm_honesty_humility_fsr %>%
  bayesplot::mcmc_trace(regex_pars = ("^fsr_"), np = nuts_params(bglm_honesty_humility_fit)) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 3000, 500), labels = ~ .x + warmup_draws + 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.02)) +
  scale_color_manual(values = traceplot_color_chart) +
  labs(
    x = "Draw",
    y = "Estimate",
    title = latex2exp::TeX("\\textbf{FIGURE}: Traceplots of estimated factor score reliabilities ($\\hat{\\rho}^2$)"),
    caption = latex2exp::TeX("\\textbf{NOTE}: Warmup draws per chain = 4,000. Post-warmup draws per chain = 2,667. Total post-warmup draws = 40,005."),
    color = "Chain"
  ) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()); bglm_honesty_humility_fsr_traceplot
ggsave("Demonstration/Plots/bglm_honesty_humility_fsr_traceplot.png", plot = bglm_honesty_humility_fsr_traceplot, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi)
ggsave("Demonstration/Plots/bglm_honesty_humility_fsr_traceplot.pdf", plot = bglm_honesty_humility_fsr_traceplot, width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
# NOTE: Traceplots looks good. Not reported for brevity.



# Construct and save BFSR plot
fig_width <- 7
fig_height <- fig_width/golden_ratio
# fig_width  <- 5
# fig_height <- fig_width*golden_ratio
n_distributions <- 3
n_latent <- 5
plot_df <- data.frame(
  fsr = c(
    bglm_honesty_humility_fsr$fsr_fairness, bglm_honesty_humility_fsr$prior_fsr_fairness, glm_honesty_humility_fit$fsr_fairness_adj,
    bglm_honesty_humility_fsr$fsr_greed_avoidance, bglm_honesty_humility_fsr$prior_fsr_greed_avoidance, glm_honesty_humility_fit$fsr_greed_avoidance_adj,
    bglm_honesty_humility_fsr$fsr_modesty, bglm_honesty_humility_fsr$prior_fsr_modesty, glm_honesty_humility_fit$fsr_modesty_adj,
    bglm_honesty_humility_fsr$fsr_sincerity, bglm_honesty_humility_fsr$prior_fsr_sincerity, glm_honesty_humility_fit$fsr_sincerity_adj,
    bglm_honesty_humility_fsr$fsr_honesty_humility, bglm_honesty_humility_fsr$prior_fsr_honesty_humility, glm_honesty_humility_fit$fsr_honesty_humility_adj
  ),
  latent = c(
    rep("Fairness", n_distributions*nrow(bglm_honesty_humility_fsr)),
    rep("Greed-avoidance", n_distributions*nrow(bglm_honesty_humility_fsr)),
    rep("Modesty", n_distributions*nrow(bglm_honesty_humility_fsr)),
    rep("Sincerity", n_distributions*nrow(bglm_honesty_humility_fsr)),
    rep("Honesty-humility", n_distributions*nrow(bglm_honesty_humility_fsr))
  ) %>% factor(levels = c("Fairness", "Greed-avoidance", "Modesty", "Sincerity", "Honesty-humility")),
  distribution = c(
    rep(rep(c("Posterior", "Prior", "Maximum-likelihood"), each = nrow(bglm_honesty_humility_fsr)), n_latent)
  ) %>% factor(levels = c("Prior", "Maximum-likelihood", "Posterior"))
)
honesty_humility_bfsr_plot <- plot_df %>% ggplot(aes(x = fsr, y = latent, fill = distribution, color = distribution, slab_linetype = distribution)) +
  theme_classic() +
  ggdist::stat_slab(
    alpha = .60,
    color = "black",
    normalize = "groups"
  ) +
  ggdist::stat_pointinterval(
    point_interval = "mode_hdci",
    .width = c(0.95),
    outline_bars = FALSE,
    color = "black",
    fill = "black",
    alpha = .70,
    position = ggplot2::position_nudge(y = -0.05)
  ) +
  ggdist::stat_pointinterval(
    point_interval = "mode_hdci",
    .width = c(0.95),
    outline_bars = FALSE,
    alpha = .70,
    position = ggplot2::position_nudge(y = -0.05)
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  # scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.05)) + # NOTE: Used instead for the truncated plot
  labs(
    # x = latex2exp::TeX("$\\hat{\\rho}^2$"),
    x = "$\\hat{\\rho}^2$", # Used instead when rendering plot as .tikz
    y = "Latent",
    fill = "Distribution",
    color = "Distribution",
    slab_linetype = "Distribution"
  ) +
  theme(
    axis.text.x.top  = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank(),
    axis.title.x = element_text(margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 10)),
    panel.grid.minor = element_blank()
  ); honesty_humility_bfsr_plot
ggsave("Demonstration/Plots/honesty_humility_bfsr_plot.png", plot = honesty_humility_bfsr_plot, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi)
ggsave("Demonstration/Plots/honesty_humility_bfsr_plot.pdf", plot = honesty_humility_bfsr_plot, width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
# tikzDevice::tikz("Demonstration/Plots/honesty_humility_bfsr_plot.tikz", standAlone = FALSE, width = fig_width, height = fig_height, sanitize = FALSE)
# print(honesty_humility_bfsr_plot)
# dev.off()
# NOTE: This is the full (i.e., non-truncated) version of figure 2 from the main paper. It is reported as figure 1 in the supplementary materials.


# Truncated plot for ease of viewing
honesty_humility_bfsr_plot <- honesty_humility_bfsr_plot +
  ggbreak::scale_x_break(c(0.099, 0.699), scales = "fixed", expand = FALSE, space = 0.25)
# ggsave("Demonstration/Plots/honesty_humility_bfsr_plot_truncated.png", plot = honesty_humility_bfsr_plot, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi)
# ggsave("Demonstration/Plots/honesty_humility_bfsr_plot_truncated.pdf", plot = honesty_humility_bfsr_plot, width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
tikzDevice::tikz("Demonstration/Plots/honesty_humility_bfsr_plot_truncated.tikz", standAlone = FALSE, width = fig_width, height = fig_height, sanitize = FALSE)
print(honesty_humility_bfsr_plot)
dev.off()
# NOTE: This is the truncated plot reported as figure 2 of the main paper.



# Save results
write.csv(bglm_honesty_humility_fsr, "Demonstration/Results/bglm_honesty_humility_fsr.csv", row.names = FALSE)





# Summarize BFSR index ----------------------------------------------------


# Load BFSR results
bglm_honesty_humility_fsr <- readr::read_csv("Demonstration/Results/bglm_honesty_humility_fsr.csv") %>% as.data.frame()


# Compute BFSR posterior summary
bfsr_posterior_summary <- data.frame(
  honesty_humility = paste0(
    bglm_honesty_humility_fsr$fsr_honesty_humility %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$fsr_honesty_humility %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$fsr_honesty_humility %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$fsr_honesty_humility %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  sincerity = paste0(
    bglm_honesty_humility_fsr$fsr_sincerity %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$fsr_sincerity %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$fsr_sincerity %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$fsr_sincerity %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  fairness = paste0(
    bglm_honesty_humility_fsr$fsr_fairness %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$fsr_fairness %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$fsr_fairness %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$fsr_fairness %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  greed_avoidance = paste0(
    bglm_honesty_humility_fsr$fsr_greed_avoidance %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$fsr_greed_avoidance %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$fsr_greed_avoidance %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$fsr_greed_avoidance %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  modesty = paste0(
    bglm_honesty_humility_fsr$fsr_modesty %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$fsr_modesty %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$fsr_modesty %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$fsr_modesty %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  )
)
# NOTE: These results for the posterior are reported in section 4.3 of the main paper.


# Compute BFSR prior summary
bfsr_prior_summary <- data.frame(
  honesty_humility = paste0(
    bglm_honesty_humility_fsr$prior_fsr_honesty_humility %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$prior_fsr_honesty_humility %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$prior_fsr_honesty_humility %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$prior_fsr_honesty_humility %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  sincerity = paste0(
    bglm_honesty_humility_fsr$prior_fsr_sincerity %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$prior_fsr_sincerity %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$prior_fsr_sincerity %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$prior_fsr_sincerity %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  fairness = paste0(
    bglm_honesty_humility_fsr$prior_fsr_fairness %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$prior_fsr_fairness %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$prior_fsr_fairness %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$prior_fsr_fairness %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  greed_avoidance = paste0(
    bglm_honesty_humility_fsr$prior_fsr_greed_avoidance %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$prior_fsr_greed_avoidance %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$prior_fsr_greed_avoidance %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$prior_fsr_greed_avoidance %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  modesty = paste0(
    bglm_honesty_humility_fsr$prior_fsr_modesty %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr$prior_fsr_modesty %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr$prior_fsr_modesty %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr$prior_fsr_modesty %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  )
)
# NOTE: These results for the prior are reported in section 4.3 of the main paper.


# Compute BFSR difference between posteriors and priors
bglm_honesty_humility_fsr_diff <- bglm_honesty_humility_fsr %>% dplyr::select(
  fsr_honesty_humility, fsr_sincerity, fsr_fairness, fsr_greed_avoidance, fsr_modesty
) - bglm_honesty_humility_fsr %>% dplyr::select(
  prior_fsr_honesty_humility, prior_fsr_sincerity, prior_fsr_fairness, prior_fsr_greed_avoidance, prior_fsr_modesty
)
bfsr_posterior_prior_difference_summary <- data.frame(
  honesty_humility = paste0(
    bglm_honesty_humility_fsr_diff$fsr_honesty_humility %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr_diff$fsr_honesty_humility %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr_diff$fsr_honesty_humility %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr_diff$fsr_honesty_humility %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  sincerity = paste0(
    bglm_honesty_humility_fsr_diff$fsr_sincerity %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr_diff$fsr_sincerity %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr_diff$fsr_sincerity %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr_diff$fsr_sincerity %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  fairness = paste0(
    bglm_honesty_humility_fsr_diff$fsr_fairness %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr_diff$fsr_fairness %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr_diff$fsr_fairness %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr_diff$fsr_fairness %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  greed_avoidance = paste0(
    bglm_honesty_humility_fsr_diff$fsr_greed_avoidance %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr_diff$fsr_greed_avoidance %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr_diff$fsr_greed_avoidance %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr_diff$fsr_greed_avoidance %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  modesty = paste0(
    bglm_honesty_humility_fsr_diff$fsr_modesty %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", bglm_honesty_humility_fsr_diff$fsr_modesty %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {bglm_honesty_humility_fsr_diff$fsr_modesty %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {bglm_honesty_humility_fsr_diff$fsr_modesty %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  )
)
# NOTE: These results for the difference between the posterior and prior are reported in section 4.3 of the main paper.





# Summarize FSR index -----------------------------------------------------


# Compute FSR summary
fsr_summary <- data.frame(
  honesty_humility = paste0(
    glm_honesty_humility_fit$fsr_honesty_humility_adj %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", glm_honesty_humility_fit$fsr_honesty_humility_adj %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {glm_honesty_humility_fit$fsr_honesty_humility_adj %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {glm_honesty_humility_fit$fsr_honesty_humility_adj %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  sincerity = paste0(
    glm_honesty_humility_fit$fsr_sincerity_adj %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", glm_honesty_humility_fit$fsr_sincerity_adj %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {glm_honesty_humility_fit$fsr_sincerity_adj %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {glm_honesty_humility_fit$fsr_sincerity_adj %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  fairness = paste0(
    glm_honesty_humility_fit$fsr_fairness_adj %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", glm_honesty_humility_fit$fsr_fairness_adj %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {glm_honesty_humility_fit$fsr_fairness_adj %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {glm_honesty_humility_fit$fsr_fairness_adj %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  greed_avoidance = paste0(
    glm_honesty_humility_fit$fsr_greed_avoidance_adj %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", glm_honesty_humility_fit$fsr_greed_avoidance_adj %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {glm_honesty_humility_fit$fsr_greed_avoidance_adj %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {glm_honesty_humility_fit$fsr_greed_avoidance_adj %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  modesty = paste0(
    glm_honesty_humility_fit$fsr_modesty_adj %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", glm_honesty_humility_fit$fsr_modesty_adj %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {glm_honesty_humility_fit$fsr_modesty_adj %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {glm_honesty_humility_fit$fsr_modesty_adj %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  )
)
# NOTE: These results for the maximum-likelihood distribution are reported in section 4.3 of the main paper when assessing robustness.


# Compute FSR difference between Bayesian and frequentist approach
fsr_bayesian_frequentist_diff <- bglm_honesty_humility_fsr %>% dplyr::select(
  fsr_honesty_humility, fsr_sincerity, fsr_fairness, fsr_greed_avoidance, fsr_modesty
) - glm_honesty_humility_fit %>% dplyr::select(
  fsr_honesty_humility_adj, fsr_sincerity_adj, fsr_fairness_adj, fsr_greed_avoidance_adj, fsr_modesty_adj
)
fsr_bayesian_frequentist_diff_summary <- data.frame(
  honesty_humility = paste0(
    fsr_bayesian_frequentist_diff$fsr_honesty_humility %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", fsr_bayesian_frequentist_diff$fsr_honesty_humility %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {fsr_bayesian_frequentist_diff$fsr_honesty_humility %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {fsr_bayesian_frequentist_diff$fsr_honesty_humility %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  sincerity = paste0(
    fsr_bayesian_frequentist_diff$fsr_sincerity %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", fsr_bayesian_frequentist_diff$fsr_sincerity %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {fsr_bayesian_frequentist_diff$fsr_sincerity %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {fsr_bayesian_frequentist_diff$fsr_sincerity %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  fairness = paste0(
    fsr_bayesian_frequentist_diff$fsr_fairness %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", fsr_bayesian_frequentist_diff$fsr_fairness %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {fsr_bayesian_frequentist_diff$fsr_fairness %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {fsr_bayesian_frequentist_diff$fsr_fairness %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  greed_avoidance = paste0(
    fsr_bayesian_frequentist_diff$fsr_greed_avoidance %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", fsr_bayesian_frequentist_diff$fsr_greed_avoidance %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {fsr_bayesian_frequentist_diff$fsr_greed_avoidance %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {fsr_bayesian_frequentist_diff$fsr_greed_avoidance %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  ),
  modesty = paste0(
    fsr_bayesian_frequentist_diff$fsr_modesty %>% mean() %>% round_down2(digits = digits$text),
    " (SE = ", fsr_bayesian_frequentist_diff$fsr_modesty %>% sd2() %>% round_up2(digits = digits$text),
    "; 95% HDI[", {fsr_bayesian_frequentist_diff$fsr_modesty %>% bayestestR::hdi()}$CI_low %>% round_down2(digits = digits$text),
    "; ", {fsr_bayesian_frequentist_diff$fsr_modesty %>% bayestestR::hdi()}$CI_high %>% round_up2(digits = digits$text), "])"
  )
)
# NOTE: These results for the difference between the posterior and maximum-likelihood are reported in section 4.3 of the main paper when assessing robustness.





# BGLM summary fit table --------------------------------------------------


# Get matrices of the estimated parameters
parameter_matrices <- get_parameter_matrices(draws = bglm_honesty_humility_draws, latent = latent, overt = overt, prior = "no")
Lambda <- parameter_matrices$Lambda
Beta <- parameter_matrices$Beta
Omega <- parameter_matrices$Omega
Delta <- parameter_matrices$Delta
Psi_Y <- parameter_matrices$Psi_Y
Psi_X <- parameter_matrices$Psi_X
Psi_Y_X <- parameter_matrices$Psi_Y_X
Psi_X_Y <- aperm(parameter_matrices$Psi_Y_X, c(2, 1, 3))
alpha_Y <- parameter_matrices$alpha_Y
alpha_X <- parameter_matrices$alpha_X

# Compute additional parameters
I_ood <- array(diag(1, o, o), dim = c(o, o, d))
I_lld <- array(diag(1, l, l), dim = c(l, l, d))
Omega_star <- invd(I_ood %-d% Omega); dimnames(Omega_star)[[1]] <- overt; dimnames(Omega_star)[[2]] <- overt
Beta_star <- invd(I_lld %-d% Beta); dimnames(Beta_star)[[1]] <- latent; dimnames(Beta_star)[[2]] <- latent
G_star <- invd(I_lld %-d% Beta %-d% (Delta %*d% Omega_star %*d% Lambda)); dimnames(G_star)[[1]] <- latent; dimnames(G_star)[[2]] <- latent
H_star <- invd(I_ood %-d% Omega %-d% (Lambda %*d% Beta_star %*d% Delta)); dimnames(H_star)[[1]] <- overt; dimnames(H_star)[[2]] <- overt
G <- Omega_star %*d% Lambda %*d% G_star; dimnames(G)[[1]] <- overt; dimnames(G)[[2]] <- latent
H <- Beta_star %*d% Delta %*d% H_star; dimnames(H)[[1]] <- latent; dimnames(H)[[2]] <- overt
Sigma_Y <- G_star %*d% Psi_Y %*d% td(G_star) + H %*d% Psi_X %*d% td(H) + G_star %*d% Psi_Y_X %*d% td(H) + H %*d% Psi_X_Y %*d% td(G_star); dimnames(Sigma_Y)[[1]] <- latent; dimnames(Sigma_Y)[[2]] <- latent
Sigma_X <- G %*d% Psi_Y %*d% td(G) + H_star %*d% Psi_X %*d% td(H_star)  + G %*d% Psi_Y_X %*d% td(H_star) + H_star %*d% Psi_X_Y %*d% td(G); dimnames(Sigma_X)[[1]] <- overt; dimnames(Sigma_X)[[2]] <- overt
Sigma_Y_X <- G_star %*d% Psi_Y %*d% td(G) + H %*d% Psi_X %*d% td(H_star) + G_star %*d% Psi_Y_X %*d% td(H_star) + H %*d% Psi_X_Y %*d% td(G); dimnames(Sigma_Y_X)[[1]] <- latent; dimnames(Sigma_Y_X)[[2]] <- overt
Sigma_X_Y <- G %*d% Psi_Y %*d% td(G_star) + H_star %*d% Psi_X %*d% td(H) + G %*d% Psi_Y_X %*d% td(H) + H_star %*d% Psi_X_Y %*d% td(G_star); dimnames(Sigma_X_Y)[[1]] <- overt; dimnames(Sigma_X_Y)[[2]] <- latent

# Standardize estimated parameters
D_Sigma_Y <- sqrtd(Sigma_Y %od% I_lld)
D_Sigma_X <- sqrtd(Sigma_X %od% I_ood)
D_Sigma_Y_inv <- invd(D_Sigma_Y)
D_Sigma_X_inv <- invd(D_Sigma_X)
Lambda_std = D_Sigma_X_inv %*d% Lambda %*d% D_Sigma_Y
Beta_std <- D_Sigma_Y_inv %*d% Beta %*d% D_Sigma_Y
Delta_std <- D_Sigma_Y_inv %*d% Delta %*d% D_Sigma_X
Omega_std <- D_Sigma_X_inv %*d% Omega %*d% D_Sigma_X
Psi_Y_std <- D_Sigma_Y_inv %*d% Psi_Y %*d% D_Sigma_Y_inv; Psi_Y_sqrt_std <- sqrtd(Psi_Y_std)
Psi_X_std <- D_Sigma_X_inv %*d% Psi_X %*d% D_Sigma_X_inv; Psi_X_sqrt_std <- sqrtd(Psi_X_std)
alpha_Y_std <- alpha_Y %*d% D_Sigma_Y_inv
alpha_X_std <- alpha_X %*d% D_Sigma_X_inv


# Create summary table
n_parameters_per_predictor <- 3
additional_parameters <- 3
bglm_honesty_humility_fit_summary <- data.frame(
  Predictor = c(
    rep("Honesty-humility", length(facets)+additional_parameters),
    rep("Sincerity", n_parameters_per_predictor*length(sincerity_indicators)+additional_parameters),
    rep("Fairness", n_parameters_per_predictor*length(fairness_indicators)+additional_parameters),
    rep("Greed-avoidance", n_parameters_per_predictor*length(greed_avoidance_indicators)+additional_parameters),
    rep("Modesty", n_parameters_per_predictor*length(modesty_indicators)+additional_parameters)
  ),
  Outcome = c(
    c("Sincerity", "Fairness", "Greed-avoidance", "Modesty"), rep("", additional_parameters),
    rep(sincerity_indicators, each = n_parameters_per_predictor), rep("", additional_parameters),
    rep(fairness_indicators, each = n_parameters_per_predictor), rep("", additional_parameters),
    rep(greed_avoidance_indicators, each = n_parameters_per_predictor), rep("", additional_parameters),
    rep(modesty_indicators, each = n_parameters_per_predictor), rep("", additional_parameters)
  ),
  Parameter = c(
    rep("$\\hat{\\lambda}$", length(facets)), "$\\hat{\\alpha}_{Y}$", "$\\hat{\\psi}_{Y}$", "$\\hat{\\rho}^2$",
    rep(c(rep(c("$\\hat{\\lambda}$", "$\\hat{\\alpha}_{X}$", "$\\hat{\\psi}_{X}$"), length(sincerity_indicators)), "$\\hat{\\alpha}_{Y}$", "$\\hat{\\psi}_{Y}$", "$\\hat{\\rho}^2$"), length(facets))
  ),
  Unstd = c(
    # Honesty-humility -> Sincerity + Fairness + Greed-avoidance + Modesty
    get_summary_table(bglm_honesty_humility_draws[,c("bsp_sincerity_mihonesty_humility", "bsp_fairness_mihonesty_humility", "bsp_greedavoidance_mihonesty_humility", "bsp_modesty_mihonesty_humility", "b_honestyhumility_Intercept", "sigma_honestyhumility")]), "-",
    
    # Sincerity -> Indicators
    get_summary_table(bglm_honesty_humility_draws[,c(
      "bsp_HSincSinc1_misincerity", "b_HSincSinc1_Intercept", "sigma_HSincSinc1",
      "bsp_HSincSinc2_misincerity", "b_HSincSinc2_Intercept", "sigma_HSincSinc2",
      "bsp_HSincSinc3_misincerity", "b_HSincSinc3_Intercept", "sigma_HSincSinc3",
      "bsp_HSincSinc4_misincerity", "b_HSincSinc4_Intercept", "sigma_HSincSinc4",
      "bsp_HSincSinc5_misincerity", "b_HSincSinc5_Intercept", "sigma_HSincSinc5",
      "bsp_HSincSinc6_misincerity", "b_HSincSinc6_Intercept", "sigma_HSincSinc6",
      "bsp_HSincSinc7_misincerity", "b_HSincSinc7_Intercept", "sigma_HSincSinc7",
      "bsp_HSincSinc8_misincerity", "b_HSincSinc8_Intercept", "sigma_HSincSinc8",
      "bsp_HSincSinc9_misincerity", "b_HSincSinc9_Intercept", "sigma_HSincSinc9",
      "bsp_HSincSinc10_misincerity", "b_HSincSinc10_Intercept", "sigma_HSincSinc10",
      "b_sincerity_Intercept", "sigma_sincerity"
    )]), "-",
    
    # Fairness -> Indicators
    get_summary_table(bglm_honesty_humility_draws[,c(
      "bsp_HFairFair1_mifairness", "b_HFairFair1_Intercept", "sigma_HFairFair1",
      "bsp_HFairFair2_mifairness", "b_HFairFair2_Intercept", "sigma_HFairFair2",
      "bsp_HFairFair3_mifairness", "b_HFairFair3_Intercept", "sigma_HFairFair3",
      "bsp_HFairFair4_mifairness", "b_HFairFair4_Intercept", "sigma_HFairFair4",
      "bsp_HFairFair5_mifairness", "b_HFairFair5_Intercept", "sigma_HFairFair5",
      "bsp_HFairFair6_mifairness", "b_HFairFair6_Intercept", "sigma_HFairFair6",
      "bsp_HFairFair7_mifairness", "b_HFairFair7_Intercept", "sigma_HFairFair7",
      "bsp_HFairFair8_mifairness", "b_HFairFair8_Intercept", "sigma_HFairFair8",
      "bsp_HFairFair9_mifairness", "b_HFairFair9_Intercept", "sigma_HFairFair9",
      "bsp_HFairFair10_mifairness", "b_HFairFair10_Intercept", "sigma_HFairFair10",
      "b_fairness_Intercept", "sigma_fairness"
    )]), "-",
    
    # Greed-avoidance -> Indicators
    get_summary_table(bglm_honesty_humility_draws[,c(
      "bsp_HGreeGree1_migreed_avoidance", "b_HGreeGree1_Intercept", "sigma_HGreeGree1",
      "bsp_HGreeGree2_migreed_avoidance", "b_HGreeGree2_Intercept", "sigma_HGreeGree2",
      "bsp_HGreeGree3_migreed_avoidance", "b_HGreeGree3_Intercept", "sigma_HGreeGree3",
      "bsp_HGreeGree4_migreed_avoidance", "b_HGreeGree4_Intercept", "sigma_HGreeGree4",
      "bsp_HGreeGree5_migreed_avoidance", "b_HGreeGree5_Intercept", "sigma_HGreeGree5",
      "bsp_HGreeGree6_migreed_avoidance", "b_HGreeGree6_Intercept", "sigma_HGreeGree6",
      "bsp_HGreeGree7_migreed_avoidance", "b_HGreeGree7_Intercept", "sigma_HGreeGree7",
      "bsp_HGreeGree8_migreed_avoidance", "b_HGreeGree8_Intercept", "sigma_HGreeGree8",
      "bsp_HGreeGree9_migreed_avoidance", "b_HGreeGree9_Intercept", "sigma_HGreeGree9",
      "bsp_HGreeGree10_migreed_avoidance", "b_HGreeGree10_Intercept", "sigma_HGreeGree10",
      "b_greedavoidance_Intercept", "sigma_greedavoidance"
    )]), "-",
    
    # Modesty -> Indicators
    get_summary_table(bglm_honesty_humility_draws[,c(
      "bsp_HModeMode1_mimodesty", "b_HModeMode1_Intercept", "sigma_HModeMode1",
      "bsp_HModeMode2_mimodesty", "b_HModeMode2_Intercept", "sigma_HModeMode2",
      "bsp_HModeMode3_mimodesty", "b_HModeMode3_Intercept", "sigma_HModeMode3",
      "bsp_HModeMode4_mimodesty", "b_HModeMode4_Intercept", "sigma_HModeMode4",
      "bsp_HModeMode5_mimodesty", "b_HModeMode5_Intercept", "sigma_HModeMode5",
      "bsp_HModeMode6_mimodesty", "b_HModeMode6_Intercept", "sigma_HModeMode6",
      "bsp_HModeMode7_mimodesty", "b_HModeMode7_Intercept", "sigma_HModeMode7",
      "bsp_HModeMode8_mimodesty", "b_HModeMode8_Intercept", "sigma_HModeMode8",
      "bsp_HModeMode9_mimodesty", "b_HModeMode9_Intercept", "sigma_HModeMode9",
      "bsp_HModeMode10_mimodesty", "b_HModeMode10_Intercept", "sigma_HModeMode10",
      "b_modesty_Intercept", "sigma_modesty"
    )]), "-"
  ),
  Std = c(
    # Honesty-humility -> Sincerity + Fairness + Greed-avoidance + Modesty
    get_summary_table(data.frame(x = Beta_std["sincerity","honesty_humility",])),
    get_summary_table(data.frame(x = Beta_std["fairness","honesty_humility",])),
    get_summary_table(data.frame(x = Beta_std["greed_avoidance","honesty_humility",])),
    get_summary_table(data.frame(x = Beta_std["modesty","honesty_humility",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"honesty_humility",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["honesty_humility","honesty_humility",])),
    get_summary_table(data.frame(x = bglm_honesty_humility_fsr[,"fsr_honesty_humility"])),
    
    # Sincerity -> Indicators
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc1","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc1",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc1","H_Sinc_Sinc1",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc2","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc2",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc2","H_Sinc_Sinc2",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc3","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc3",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc3","H_Sinc_Sinc3",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc4","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc4",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc4","H_Sinc_Sinc4",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc5","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc5",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc5","H_Sinc_Sinc5",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc6","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc6",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc6","H_Sinc_Sinc6",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc7","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc7",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc7","H_Sinc_Sinc7",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc8","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc8",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc8","H_Sinc_Sinc8",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc9","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc9",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc9","H_Sinc_Sinc9",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc10","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc10",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc10","H_Sinc_Sinc10",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"sincerity",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["sincerity","sincerity",])),
    get_summary_table(data.frame(x = bglm_honesty_humility_fsr[,"fsr_sincerity"])),
    
    # Fairness -> Indicators
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair1","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair1",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair1","H_Fair_Fair1",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair2","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair2",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair2","H_Fair_Fair2",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair3","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair3",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair3","H_Fair_Fair3",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair4","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair4",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair4","H_Fair_Fair4",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair5","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair5",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair5","H_Fair_Fair5",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair6","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair6",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair6","H_Fair_Fair6",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair7","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair7",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair7","H_Fair_Fair7",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair8","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair8",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair8","H_Fair_Fair8",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair9","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair9",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair9","H_Fair_Fair9",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair10","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair10",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair10","H_Fair_Fair10",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"fairness",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["fairness","fairness",])),
    get_summary_table(data.frame(x = bglm_honesty_humility_fsr[,"fsr_fairness"])),
    
    # Greed-avoidance -> Indicators
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree1","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree1",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree1","H_Gree_Gree1",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree2","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree2",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree2","H_Gree_Gree2",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree3","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree3",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree3","H_Gree_Gree3",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree4","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree4",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree4","H_Gree_Gree4",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree5","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree5",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree5","H_Gree_Gree5",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree6","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree6",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree6","H_Gree_Gree6",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree7","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree7",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree7","H_Gree_Gree7",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree8","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree8",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree8","H_Gree_Gree8",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree9","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree9",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree9","H_Gree_Gree9",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree10","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree10",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree10","H_Gree_Gree10",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"greed_avoidance",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["greed_avoidance","greed_avoidance",])),
    get_summary_table(data.frame(x = bglm_honesty_humility_fsr[,"fsr_greed_avoidance"])),
    
    # Modesty -> Indicators
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode1","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode1",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode1","H_Mode_Mode1",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode2","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode2",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode2","H_Mode_Mode2",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode3","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode3",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode3","H_Mode_Mode3",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode4","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode4",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode4","H_Mode_Mode4",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode5","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode5",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode5","H_Mode_Mode5",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode6","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode6",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode6","H_Mode_Mode6",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode7","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode7",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode7","H_Mode_Mode7",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode8","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode8",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode8","H_Mode_Mode8",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode9","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode9",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode9","H_Mode_Mode9",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode10","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode10",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode10","H_Mode_Mode10",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"modesty",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["modesty","modesty",])),
    get_summary_table(data.frame(x = bglm_honesty_humility_fsr[,"fsr_modesty"]))
  ),
  Rhat = c(
    # Honesty-humility -> Sincerity + Fairness + Greed-avoidance + Modesty
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_sincerity_mihonesty_humility"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_greedavoidance_mihonesty_humility"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_modesty_mihonesty_humility"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_honestyhumility"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$rhat[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_honesty_humility"] %>% round_up2(digits = digits$table),
    
    # Sincerity -> Indicators
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc1_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc1_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc1"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc2_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc2_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc2"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc3_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc3_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc3"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc4_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc4_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc4"] %>% round_up2(digits = digits$table),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc5"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc6_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc6_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc6"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc7_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc7_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc7"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc8_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc8_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc8"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc9_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc9_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc9"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc10_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc10_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc10"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_sincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$rhat[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_sincerity"] %>% round_up2(digits = digits$table),
    
    # Fairness -> Indicators
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair1_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair1_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair1"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair2_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair2_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair2"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair3_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair3_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair3"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair4_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair4_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair4"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair5_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair5_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair5"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair6_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair6_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair6"] %>% round_up2(digits = digits$table),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair7"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair8_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair8_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair8"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair9_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair9_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair9"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair10_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair10_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair10"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_fairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$rhat[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_fairness"] %>% round_up2(digits = digits$table),
    
    # Greed-avoidance -> Indicators
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree1_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree1_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree1"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree2_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree2_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree2"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree3_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree3_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree3"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree4_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree4_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree4"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree5_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree5_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree5"] %>% round_up2(digits = digits$table),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree6"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree7_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree7_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree7"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree8_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree8_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree8"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree9_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree9_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree9"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree10_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree10_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree10"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_greedavoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$rhat[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_greed_avoidance"] %>% round_up2(digits = digits$table),
    
    # Modesty -> Indicators
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode1"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode2_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode2_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode2"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode3_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode3_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode3"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode4_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode4_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode4"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode5_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode5_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode5"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode6_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode6_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode6"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode7_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode7_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode7"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode8_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode8_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode8"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode9_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode9_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode9"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode10_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode10_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode10"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$rhat[bglm_honesty_humility_fit_diagnostics$variable == "sigma_modesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$rhat[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_modesty"] %>% round_up2(digits = digits$table)
  ),
  MCSE = c(
    # Honesty-humility -> Sincerity + Fairness + Greed-avoidance + Modesty
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_sincerity_mihonesty_humility"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_greedavoidance_mihonesty_humility"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_modesty_mihonesty_humility"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_honestyhumility"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$mcse_mean[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_honesty_humility"] %>% round_up2(digits = digits$table),
    
    # Sincerity -> Indicators
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc1_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc1_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc1"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc2_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc2_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc2"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc3_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc3_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc3"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc4_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc4_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc4"] %>% round_up2(digits = digits$table),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc5"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc6_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc6_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc6"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc7_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc7_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc7"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc8_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc8_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc8"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc9_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc9_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc9"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc10_misincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc10_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc10"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_sincerity"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$mcse_mean[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_sincerity"] %>% round_up2(digits = digits$table),
    
    # Fairness -> Indicators
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair1_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair1_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair1"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair2_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair2_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair2"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair3_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair3_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair3"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair4_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair4_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair4"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair5_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair5_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair5"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair6_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair6_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair6"] %>% round_up2(digits = digits$table),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair7"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair8_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair8_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair8"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair9_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair9_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair9"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair10_mifairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair10_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair10"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_fairness"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$mcse_mean[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_fairness"] %>% round_up2(digits = digits$table),
    
    # Greed-avoidance -> Indicators
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree1_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree1_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree1"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree2_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree2_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree2"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree3_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree3_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree3"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree4_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree4_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree4"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree5_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree5_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree5"] %>% round_up2(digits = digits$table),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree6"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree7_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree7_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree7"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree8_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree8_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree8"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree9_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree9_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree9"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree10_migreed_avoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree10_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree10"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_greedavoidance"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$mcse_mean[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_greed_avoidance"] %>% round_up2(digits = digits$table),
    
    # Modesty -> Indicators
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode1"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode2_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode2_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode2"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode3_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode3_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode3"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode4_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode4_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode4"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode5_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode5_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode5"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode6_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode6_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode6"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode7_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode7_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode7"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode8_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode8_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode8"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode9_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode9_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode9"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode10_mimodesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode10_Intercept"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode10"] %>% round_up2(digits = digits$table),
    "-",
    bglm_honesty_humility_fit_diagnostics$mcse_mean[bglm_honesty_humility_fit_diagnostics$variable == "sigma_modesty"] %>% round_up2(digits = digits$table),
    bglm_honesty_humility_fsr_diagnostics$mcse_mean[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_modesty"] %>% round_up2(digits = digits$table)
  ),
  ESS_bulk = c(
    # Honesty-humility -> Sincerity + Fairness + Greed-avoidance + Modesty
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_sincerity_mihonesty_humility"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_greedavoidance_mihonesty_humility"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_modesty_mihonesty_humility"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_honestyhumility"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_bulk[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_honesty_humility"] %>% qformat(),
    
    # Sincerity -> Indicators
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc1_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc1_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc1"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc2_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc2_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc2"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc3_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc3_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc3"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc4_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc4_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc4"] %>% qformat(),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc5"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc6_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc6_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc6"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc7_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc7_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc7"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc8_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc8_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc8"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc9_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc9_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc9"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc10_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc10_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc10"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_sincerity"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_bulk[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_sincerity"] %>% qformat(),
    
    # Fairness -> Indicators
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair1_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair1_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair1"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair2_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair2_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair2"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair3_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair3_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair3"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair4_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair4_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair4"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair5_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair5_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair5"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair6_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair6_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair6"] %>% qformat(),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair7"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair8_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair8_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair8"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair9_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair9_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair9"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair10_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair10_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair10"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_fairness"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_bulk[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_fairness"] %>% qformat(),
    
    # Greed-avoidance -> Indicators
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree1_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree1_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree1"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree2_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree2_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree2"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree3_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree3_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree3"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree4_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree4_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree4"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree5_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree5_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree5"] %>% qformat(),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree6"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree7_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree7_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree7"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree8_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree8_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree8"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree9_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree9_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree9"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree10_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree10_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree10"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_greedavoidance"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_bulk[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_greed_avoidance"] %>% qformat(),
    
    # Modesty -> Indicators
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode1"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode2_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode2_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode2"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode3_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode3_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode3"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode4_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode4_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode4"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode5_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode5_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode5"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode6_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode6_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode6"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode7_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode7_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode7"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode8_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode8_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode8"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode9_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode9_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode9"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode10_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode10_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode10"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_bulk[bglm_honesty_humility_fit_diagnostics$variable == "sigma_modesty"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_bulk[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_modesty"] %>% qformat()
  ),
  ESS_tail = c(
    # Honesty-humility -> Sincerity + Fairness + Greed-avoidance + Modesty
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_sincerity_mihonesty_humility"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_greedavoidance_mihonesty_humility"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_modesty_mihonesty_humility"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_honestyhumility"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_tail[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_honesty_humility"] %>% qformat(),
    
    # Sincerity -> Indicators
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc1_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc1_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc1"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc2_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc2_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc2"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc3_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc3_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc3"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc4_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc4_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc4"] %>% qformat(),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc5"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc6_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc6_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc6"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc7_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc7_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc7"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc8_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc8_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc8"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc9_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc9_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc9"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HSincSinc10_misincerity"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HSincSinc10_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HSincSinc10"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_sincerity"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_tail[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_sincerity"] %>% qformat(),
    
    # Fairness -> Indicators
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair1_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair1_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair1"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair2_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair2_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair2"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair3_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair3_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair3"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair4_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair4_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair4"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair5_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair5_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair5"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair6_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair6_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair6"] %>% qformat(),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair7"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair8_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair8_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair8"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair9_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair9_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair9"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HFairFair10_mifairness"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HFairFair10_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HFairFair10"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_fairness"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_tail[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_fairness"] %>% qformat(),
    
    # Greed-avoidance -> Indicators
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree1_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree1_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree1"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree2_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree2_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree2"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree3_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree3_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree3"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree4_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree4_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree4"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree5_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree5_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree5"] %>% qformat(),
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree6"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree7_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree7_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree7"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree8_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree8_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree8"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree9_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree9_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree9"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HGreeGree10_migreed_avoidance"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HGreeGree10_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HGreeGree10"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_greedavoidance"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_tail[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_greed_avoidance"] %>% qformat(),
    
    # Modesty -> Indicators
    "-",
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode1"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode2_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode2_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode2"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode3_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode3_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode3"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode4_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode4_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode4"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode5_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode5_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode5"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode6_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode6_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode6"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode7_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode7_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode7"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode8_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode8_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode8"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode9_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode9_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode9"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "bsp_HModeMode10_mimodesty"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "b_HModeMode10_Intercept"] %>% qformat(),
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_HModeMode10"] %>% qformat(),
    "-",
    bglm_honesty_humility_fit_diagnostics$ess_tail[bglm_honesty_humility_fit_diagnostics$variable == "sigma_modesty"] %>% qformat(),
    bglm_honesty_humility_fsr_diagnostics$ess_tail[bglm_honesty_humility_fsr_diagnostics$variable == "fsr_modesty"] %>% qformat()
  ),
  Prior = c(
    rep("$\\mathcal{N}(0, .304)$", length(facets)), "$\\delta(0)$", "$\\Gamma(1, 3.29)$", "-",
    rep(c(rep(c("$\\mathcal{N}(0, .304)$", "$\\mathcal{N}(0, 1)$", "$\\Gamma(12.7, 11.7)$"), length(sincerity_indicators)), "$\\delta(0)$", "$\\Gamma(1, 3.29)$", "-"), length(facets))
  )
)
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Honesty-humility" & bglm_honesty_humility_fit_summary$Outcome == "Fairness" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\lambda}$"] <- "$\\delta(1)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Honesty-humility" & bglm_honesty_humility_fit_summary$Outcome == "Fairness" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{X}$"] <- "$\\delta(0)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Sincerity" & bglm_honesty_humility_fit_summary$Outcome == "H_Sinc_Sinc5" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\lambda}$"] <- "$\\delta(1)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Sincerity" & bglm_honesty_humility_fit_summary$Outcome == "H_Sinc_Sinc5" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{X}$"] <- "$\\delta(0)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Fairness" & bglm_honesty_humility_fit_summary$Outcome == "H_Fair_Fair7" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\lambda}$"] <- "$\\delta(1)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Fairness" & bglm_honesty_humility_fit_summary$Outcome == "H_Fair_Fair7" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{X}$"] <- "$\\delta(0)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Greed-avoidance" & bglm_honesty_humility_fit_summary$Outcome == "H_Gree_Gree6" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\lambda}$"] <- "$\\delta(1)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Greed-avoidance" & bglm_honesty_humility_fit_summary$Outcome == "H_Gree_Gree6" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{X}$"] <- "$\\delta(0)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Modesty" & bglm_honesty_humility_fit_summary$Outcome == "H_Mode_Mode1" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\lambda}$"] <- "$\\delta(1)$"
bglm_honesty_humility_fit_summary$Prior[bglm_honesty_humility_fit_summary$Predictor == "Modesty" & bglm_honesty_humility_fit_summary$Outcome == "H_Mode_Mode1" & bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{X}$"] <- "$\\delta(0)$"
bglm_honesty_humility_fit_summary$Parameter[bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\lambda}$" & bglm_honesty_humility_fit_summary$Prior == "$\\delta(1)$"] <- "$\\hat{\\lambda}^{\\dagger}$"
bglm_honesty_humility_fit_summary$Parameter[bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{Y}$" & bglm_honesty_humility_fit_summary$Prior == "$\\delta(0)$"] <- "$\\hat{\\alpha}_{Y}^{\\dagger}$"
bglm_honesty_humility_fit_summary$Parameter[bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{X}$" & bglm_honesty_humility_fit_summary$Prior == "$\\delta(0)$"] <- "$\\hat{\\alpha}_{X}^{\\dagger}$"
bglm_honesty_humility_fit_summary$Unstd[bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\lambda}^{\\dagger}$" & bglm_honesty_humility_fit_summary$Prior == "$\\delta(1)$"] <- "1 (0) [1; 1]"
bglm_honesty_humility_fit_summary$Unstd[bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{Y}^{\\dagger}$" & bglm_honesty_humility_fit_summary$Prior == "$\\delta(0)$"] <- "0 (0) [0; 0]"
bglm_honesty_humility_fit_summary$Unstd[bglm_honesty_humility_fit_summary$Parameter == "$\\hat{\\alpha}_{X}^{\\dagger}$" & bglm_honesty_humility_fit_summary$Prior == "$\\delta(0)$"] <- "0 (0) [0; 0]"

# View summary table
bglm_honesty_humility_fit_summary %>%
  knitr::kable(
    # format = "latex",
    digits = 3,
    row.names = FALSE,
    # col.names = c("$\\text{Indicator}$", "$\\text{Parameter}$", "$\\text{Posterior - Unstandardized}$", "$\\text{Posterior - Standardized}$", "$\\hat{R}$", "$\\text{MCSE}$", "$\\text{ESS}_{\\text{Bulk}}$", "$\\text{ESS}_{\\text{Tail}}$", "$\\text{Prior}$"),
    format.args = list(big.mark = ",", scientific = FALSE),
    align = "lcccccccc",
    escape = FALSE
  ) %>%
  kableExtra::kable_styling(
    full_width = FALSE,
    position = "center",
    fixed_thead = TRUE,
    latex_options = c("striped", "repeat_header", "HOLD_position"),
    font_size = 12
    
  ) %>%
  kableExtra::row_spec(0, bold = TRUE) %>%
  kableExtra::column_spec(1, bold = TRUE)
# NOTE: This summary table of the results for the Bayesian model is reported as table 1 in the supplementary materials.



# Save summary table
openxlsx::write.xlsx(bglm_honesty_humility_fit_summary, "Demonstration/Results/bglm_honesty_humility_fit_summary.xlsx")





# GLM summary fit table ---------------------------------------------------


# Get matrices of the estimated parameters
parameter_matrices <- get_parameter_matrices(draws = glm_honesty_humility_fit, latent = latent, overt = overt, prior = "no")
Lambda <- parameter_matrices$Lambda
Beta <- parameter_matrices$Beta
Omega <- parameter_matrices$Omega
Delta <- parameter_matrices$Delta
Psi_Y <- parameter_matrices$Psi_Y; Psi_Y_sqrt <- sqrtd(Psi_Y)
Psi_X <- parameter_matrices$Psi_X; Psi_X_sqrt <- sqrtd(Psi_X)
Psi_Y_X <- parameter_matrices$Psi_Y_X
Psi_X_Y <- aperm(parameter_matrices$Psi_Y_X, c(2, 1, 3))
alpha_Y <- parameter_matrices$alpha_Y
alpha_X <- parameter_matrices$alpha_X

# Compute additional parameters
I_ood <- array(diag(1, o, o), dim = c(o, o, d))
I_lld <- array(diag(1, l, l), dim = c(l, l, d))
Omega_star <- invd(I_ood %-d% Omega); dimnames(Omega_star)[[1]] <- overt; dimnames(Omega_star)[[2]] <- overt
Beta_star <- invd(I_lld %-d% Beta); dimnames(Beta_star)[[1]] <- latent; dimnames(Beta_star)[[2]] <- latent
G_star <- invd(I_lld %-d% Beta %-d% (Delta %*d% Omega_star %*d% Lambda)); dimnames(G_star)[[1]] <- latent; dimnames(G_star)[[2]] <- latent
H_star <- invd(I_ood %-d% Omega %-d% (Lambda %*d% Beta_star %*d% Delta)); dimnames(H_star)[[1]] <- overt; dimnames(H_star)[[2]] <- overt
G <- Omega_star %*d% Lambda %*d% G_star; dimnames(G)[[1]] <- overt; dimnames(G)[[2]] <- latent
H <- Beta_star %*d% Delta %*d% H_star; dimnames(H)[[1]] <- latent; dimnames(H)[[2]] <- overt
Sigma_Y <- G_star %*d% Psi_Y %*d% td(G_star) + H %*d% Psi_X %*d% td(H) + G_star %*d% Psi_Y_X %*d% td(H) + H %*d% Psi_X_Y %*d% td(G_star); dimnames(Sigma_Y)[[1]] <- latent; dimnames(Sigma_Y)[[2]] <- latent
Sigma_X <- G %*d% Psi_Y %*d% td(G) + H_star %*d% Psi_X %*d% td(H_star)  + G %*d% Psi_Y_X %*d% td(H_star) + H_star %*d% Psi_X_Y %*d% td(G); dimnames(Sigma_X)[[1]] <- overt; dimnames(Sigma_X)[[2]] <- overt
Sigma_Y_X <- G_star %*d% Psi_Y %*d% td(G) + H %*d% Psi_X %*d% td(H_star) + G_star %*d% Psi_Y_X %*d% td(H_star) + H %*d% Psi_X_Y %*d% td(G); dimnames(Sigma_Y_X)[[1]] <- latent; dimnames(Sigma_Y_X)[[2]] <- overt
Sigma_X_Y <- G %*d% Psi_Y %*d% td(G_star) + H_star %*d% Psi_X %*d% td(H) + G %*d% Psi_Y_X %*d% td(H) + H_star %*d% Psi_X_Y %*d% td(G_star); dimnames(Sigma_X_Y)[[1]] <- overt; dimnames(Sigma_X_Y)[[2]] <- latent

# Standardize estimated parameters
D_Sigma_Y <- sqrtd(Sigma_Y %od% I_lld)
D_Sigma_X <- sqrtd(Sigma_X %od% I_ood)
D_Sigma_Y_inv <- invd(D_Sigma_Y)
D_Sigma_X_inv <- invd(D_Sigma_X)
Lambda_std = D_Sigma_X_inv %*d% Lambda %*d% D_Sigma_Y
Beta_std <- D_Sigma_Y_inv %*d% Beta %*d% D_Sigma_Y
Delta_std <- D_Sigma_Y_inv %*d% Delta %*d% D_Sigma_X
Omega_std <- D_Sigma_X_inv %*d% Omega %*d% D_Sigma_X
Psi_Y_std <- D_Sigma_Y_inv %*d% Psi_Y %*d% D_Sigma_Y_inv; Psi_Y_sqrt_std <- sqrtd(Psi_Y_std)
Psi_X_std <- D_Sigma_X_inv %*d% Psi_X %*d% D_Sigma_X_inv; Psi_X_sqrt_std <- sqrtd(Psi_X_std)
alpha_Y_std <- alpha_Y %*d% D_Sigma_Y_inv
alpha_X_std <- alpha_X %*d% D_Sigma_X_inv

# Create summary table
n_parameters_per_predictor <- 3
additional_parameters <- 3
glm_honesty_humility_fit_summary <- data.frame(
  Estimate = c(
    rep("Honesty-humility", length(facets)+additional_parameters),
    rep("Sincerity", n_parameters_per_predictor*length(sincerity_indicators)+additional_parameters),
    rep("Fairness", n_parameters_per_predictor*length(fairness_indicators)+additional_parameters),
    rep("Greed-avoidance", n_parameters_per_predictor*length(greed_avoidance_indicators)+additional_parameters),
    rep("Modesty", n_parameters_per_predictor*length(modesty_indicators)+additional_parameters)
  ),
  Outcome = c(
    c("Sincerity", "Fairness", "Greed-avoidance", "Modesty"), rep("", additional_parameters),
    rep(sincerity_indicators, each = n_parameters_per_predictor), rep("", additional_parameters),
    rep(fairness_indicators, each = n_parameters_per_predictor), rep("", additional_parameters),
    rep(greed_avoidance_indicators, each = n_parameters_per_predictor), rep("", additional_parameters),
    rep(modesty_indicators, each = n_parameters_per_predictor), rep("", additional_parameters)
  ),
  Parameter = c(
    rep("$\\hat{\\lambda}$", length(facets)), "$\\hat{\\alpha}_{Y}$", "$\\hat{\\psi}_{Y}$", "$\\hat{\\rho}^2$",
    rep(c(rep(c("$\\hat{\\lambda}$", "$\\hat{\\alpha}_{X}$", "$\\hat{\\psi}_{X}$"), length(sincerity_indicators)), "$\\hat{\\alpha}_{Y}$", "$\\hat{\\psi}_{Y}$", "$\\hat{\\rho}^2$"), length(facets))
  ),
  Unstd = c(
    # Honesty-humility -> Sincerity + Fairness + Greed-avoidance + Modesty
    get_summary_table(glm_honesty_humility_fit[,c("bsp_sincerity_mihonesty_humility", "bsp_fairness_mihonesty_humility", "bsp_greedavoidance_mihonesty_humility", "bsp_modesty_mihonesty_humility", "b_honestyhumility_Intercept", "sigma_honestyhumility")]), "-",
    
    # Sincerity -> Indicators
    get_summary_table(glm_honesty_humility_fit[,c(
      "bsp_HSincSinc1_misincerity", "b_HSincSinc1_Intercept", "sigma_HSincSinc1",
      "bsp_HSincSinc2_misincerity", "b_HSincSinc2_Intercept", "sigma_HSincSinc2",
      "bsp_HSincSinc3_misincerity", "b_HSincSinc3_Intercept", "sigma_HSincSinc3",
      "bsp_HSincSinc4_misincerity", "b_HSincSinc4_Intercept", "sigma_HSincSinc4",
      "bsp_HSincSinc5_misincerity", "b_HSincSinc5_Intercept", "sigma_HSincSinc5",
      "bsp_HSincSinc6_misincerity", "b_HSincSinc6_Intercept", "sigma_HSincSinc6",
      "bsp_HSincSinc7_misincerity", "b_HSincSinc7_Intercept", "sigma_HSincSinc7",
      "bsp_HSincSinc8_misincerity", "b_HSincSinc8_Intercept", "sigma_HSincSinc8",
      "bsp_HSincSinc9_misincerity", "b_HSincSinc9_Intercept", "sigma_HSincSinc9",
      "bsp_HSincSinc10_misincerity", "b_HSincSinc10_Intercept", "sigma_HSincSinc10",
      "b_sincerity_Intercept", "sigma_sincerity"
    )]), "-",
    
    # Fairness -> Indicators
    get_summary_table(glm_honesty_humility_fit[,c(
      "bsp_HFairFair1_mifairness", "b_HFairFair1_Intercept", "sigma_HFairFair1",
      "bsp_HFairFair2_mifairness", "b_HFairFair2_Intercept", "sigma_HFairFair2",
      "bsp_HFairFair3_mifairness", "b_HFairFair3_Intercept", "sigma_HFairFair3",
      "bsp_HFairFair4_mifairness", "b_HFairFair4_Intercept", "sigma_HFairFair4",
      "bsp_HFairFair5_mifairness", "b_HFairFair5_Intercept", "sigma_HFairFair5",
      "bsp_HFairFair6_mifairness", "b_HFairFair6_Intercept", "sigma_HFairFair6",
      "bsp_HFairFair7_mifairness", "b_HFairFair7_Intercept", "sigma_HFairFair7",
      "bsp_HFairFair8_mifairness", "b_HFairFair8_Intercept", "sigma_HFairFair8",
      "bsp_HFairFair9_mifairness", "b_HFairFair9_Intercept", "sigma_HFairFair9",
      "bsp_HFairFair10_mifairness", "b_HFairFair10_Intercept", "sigma_HFairFair10",
      "b_fairness_Intercept", "sigma_fairness"
    )]), "-",
    
    # Greed-avoidance -> Indicators
    get_summary_table(glm_honesty_humility_fit[,c(
      "bsp_HGreeGree1_migreed_avoidance", "b_HGreeGree1_Intercept", "sigma_HGreeGree1",
      "bsp_HGreeGree2_migreed_avoidance", "b_HGreeGree2_Intercept", "sigma_HGreeGree2",
      "bsp_HGreeGree3_migreed_avoidance", "b_HGreeGree3_Intercept", "sigma_HGreeGree3",
      "bsp_HGreeGree4_migreed_avoidance", "b_HGreeGree4_Intercept", "sigma_HGreeGree4",
      "bsp_HGreeGree5_migreed_avoidance", "b_HGreeGree5_Intercept", "sigma_HGreeGree5",
      "bsp_HGreeGree6_migreed_avoidance", "b_HGreeGree6_Intercept", "sigma_HGreeGree6",
      "bsp_HGreeGree7_migreed_avoidance", "b_HGreeGree7_Intercept", "sigma_HGreeGree7",
      "bsp_HGreeGree8_migreed_avoidance", "b_HGreeGree8_Intercept", "sigma_HGreeGree8",
      "bsp_HGreeGree9_migreed_avoidance", "b_HGreeGree9_Intercept", "sigma_HGreeGree9",
      "bsp_HGreeGree10_migreed_avoidance", "b_HGreeGree10_Intercept", "sigma_HGreeGree10",
      "b_greedavoidance_Intercept", "sigma_greedavoidance"
    )]), "-",
    
    # Modesty -> Indicators
    get_summary_table(glm_honesty_humility_fit[,c(
      "bsp_HModeMode1_mimodesty", "b_HModeMode1_Intercept", "sigma_HModeMode1",
      "bsp_HModeMode2_mimodesty", "b_HModeMode2_Intercept", "sigma_HModeMode2",
      "bsp_HModeMode3_mimodesty", "b_HModeMode3_Intercept", "sigma_HModeMode3",
      "bsp_HModeMode4_mimodesty", "b_HModeMode4_Intercept", "sigma_HModeMode4",
      "bsp_HModeMode5_mimodesty", "b_HModeMode5_Intercept", "sigma_HModeMode5",
      "bsp_HModeMode6_mimodesty", "b_HModeMode6_Intercept", "sigma_HModeMode6",
      "bsp_HModeMode7_mimodesty", "b_HModeMode7_Intercept", "sigma_HModeMode7",
      "bsp_HModeMode8_mimodesty", "b_HModeMode8_Intercept", "sigma_HModeMode8",
      "bsp_HModeMode9_mimodesty", "b_HModeMode9_Intercept", "sigma_HModeMode9",
      "bsp_HModeMode10_mimodesty", "b_HModeMode10_Intercept", "sigma_HModeMode10",
      "b_modesty_Intercept", "sigma_modesty"
    )]), "-"
  ),
  Std = c(
    # Honesty-humility -> Sincerity + Fairness + Greed-avoidance + Modesty
    get_summary_table(data.frame(x = Beta_std["sincerity","honesty_humility",])),
    get_summary_table(data.frame(x = Beta_std["fairness","honesty_humility",])),
    get_summary_table(data.frame(x = Beta_std["greed_avoidance","honesty_humility",])),
    get_summary_table(data.frame(x = Beta_std["modesty","honesty_humility",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"honesty_humility",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["honesty_humility","honesty_humility",])),
    get_summary_table(data.frame(x = glm_honesty_humility_fit[,"fsr_honesty_humility_adj"])),
    
    # Sincerity -> Indicators
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc1","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc1",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc1","H_Sinc_Sinc1",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc2","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc2",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc2","H_Sinc_Sinc2",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc3","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc3",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc3","H_Sinc_Sinc3",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc4","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc4",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc4","H_Sinc_Sinc4",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc5","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc5",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc5","H_Sinc_Sinc5",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc6","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc6",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc6","H_Sinc_Sinc6",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc7","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc7",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc7","H_Sinc_Sinc7",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc8","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc8",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc8","H_Sinc_Sinc8",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc9","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc9",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc9","H_Sinc_Sinc9",])),
    get_summary_table(data.frame(x = Lambda_std["H_Sinc_Sinc10","sincerity",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Sinc_Sinc10",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Sinc_Sinc10","H_Sinc_Sinc10",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"sincerity",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["sincerity","sincerity",])),
    get_summary_table(data.frame(x = glm_honesty_humility_fit[,"fsr_sincerity_adj"])),
    
    # Fairness -> Indicators
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair1","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair1",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair1","H_Fair_Fair1",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair2","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair2",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair2","H_Fair_Fair2",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair3","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair3",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair3","H_Fair_Fair3",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair4","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair4",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair4","H_Fair_Fair4",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair5","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair5",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair5","H_Fair_Fair5",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair6","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair6",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair6","H_Fair_Fair6",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair7","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair7",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair7","H_Fair_Fair7",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair8","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair8",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair8","H_Fair_Fair8",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair9","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair9",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair9","H_Fair_Fair9",])),
    get_summary_table(data.frame(x = Lambda_std["H_Fair_Fair10","fairness",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Fair_Fair10",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Fair_Fair10","H_Fair_Fair10",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"fairness",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["fairness","fairness",])),
    get_summary_table(data.frame(x = glm_honesty_humility_fit[,"fsr_fairness_adj"])),
    
    # Greed-avoidance -> Indicators
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree1","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree1",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree1","H_Gree_Gree1",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree2","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree2",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree2","H_Gree_Gree2",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree3","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree3",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree3","H_Gree_Gree3",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree4","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree4",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree4","H_Gree_Gree4",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree5","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree5",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree5","H_Gree_Gree5",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree6","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree6",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree6","H_Gree_Gree6",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree7","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree7",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree7","H_Gree_Gree7",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree8","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree8",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree8","H_Gree_Gree8",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree9","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree9",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree9","H_Gree_Gree9",])),
    get_summary_table(data.frame(x = Lambda_std["H_Gree_Gree10","greed_avoidance",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Gree_Gree10",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Gree_Gree10","H_Gree_Gree10",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"greed_avoidance",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["greed_avoidance","greed_avoidance",])),
    get_summary_table(data.frame(x = glm_honesty_humility_fit[,"fsr_greed_avoidance_adj"])),
    
    # Modesty -> Indicators
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode1","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode1",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode1","H_Mode_Mode1",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode2","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode2",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode2","H_Mode_Mode2",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode3","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode3",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode3","H_Mode_Mode3",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode4","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode4",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode4","H_Mode_Mode4",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode5","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode5",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode5","H_Mode_Mode5",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode6","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode6",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode6","H_Mode_Mode6",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode7","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode7",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode7","H_Mode_Mode7",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode8","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode8",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode8","H_Mode_Mode8",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode9","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode9",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode9","H_Mode_Mode9",])),
    get_summary_table(data.frame(x = Lambda_std["H_Mode_Mode10","modesty",])),
    get_summary_table(data.frame(x = alpha_X_std[,"H_Mode_Mode10",])),
    get_summary_table(data.frame(x = Psi_X_sqrt_std["H_Mode_Mode10","H_Mode_Mode10",])),
    get_summary_table(data.frame(x = alpha_Y_std[,"modesty",])),
    get_summary_table(data.frame(x = Psi_Y_sqrt_std["modesty","modesty",])),
    get_summary_table(data.frame(x = glm_honesty_humility_fit[,"fsr_modesty_adj"]))
  )
)
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\lambda$" & glm_honesty_humility_fit_summary$Outcome == "Fairness"] <- "$\\lambda^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\alpha_{X}$" & glm_honesty_humility_fit_summary$Outcome == "Fairness"] <- "$\\alpha_{X}^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\lambda$" & glm_honesty_humility_fit_summary$Outcome == "H_Sinc_Sinc5"] <- "$\\lambda^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\alpha_{X}$" & glm_honesty_humility_fit_summary$Outcome == "H_Sinc_Sinc5"] <- "$\\alpha_{X}^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\lambda$" & glm_honesty_humility_fit_summary$Outcome == "H_Fair_Fair7"] <- "$\\lambda^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\alpha_{X}$" & glm_honesty_humility_fit_summary$Outcome == "H_Fair_Fair7"] <- "$\\alpha_{X}^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\lambda$" & glm_honesty_humility_fit_summary$Outcome == "H_Gree_Gree6"] <- "$\\lambda^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\alpha_{X}$" & glm_honesty_humility_fit_summary$Outcome == "H_Gree_Gree6"] <- "$\\alpha_{X}^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\lambda$" & glm_honesty_humility_fit_summary$Outcome == "H_Mode_Mode1"] <- "$\\lambda^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\alpha_{X}$" & glm_honesty_humility_fit_summary$Outcome == "H_Mode_Mode1"] <- "$\\alpha_{X}^{\\dagger}$"
glm_honesty_humility_fit_summary$Parameter[glm_honesty_humility_fit_summary$Parameter == "$\\alpha_{Y}$"] <- "$\\alpha_{Y}^{\\dagger}$"
glm_honesty_humility_fit_summary$Unstd[glm_honesty_humility_fit_summary$Parameter == "$\\lambda^{\\dagger}$"] <- "1 (0) [1; 1]"
glm_honesty_humility_fit_summary$Unstd[glm_honesty_humility_fit_summary$Parameter == "$\\alpha_{Y}^{\\dagger}$"] <- "0 (0) [0; 0]"
glm_honesty_humility_fit_summary$Unstd[glm_honesty_humility_fit_summary$Parameter == "$\\alpha_{X}^{\\dagger}$"] <- "0 (0) [0; 0]"

# View summary table
glm_honesty_humility_fit_summary %>%
  knitr::kable(
    # format = "latex",
    digits = 3,
    row.names = FALSE,
    # col.names = c("$\\text{Indicator}$", "$\\text{Parameter}$", "$\\text{Posterior - Unstandardized}$", "$\\text{Posterior - Standardized}$", "$\\hat{R}$", "$\\text{MCSE}$", "$\\text{ESS}_{\\text{Bulk}}$", "$\\text{ESS}_{\\text{Tail}}$", "$\\text{Prior}$"),
    format.args = list(big.mark = ",", scientific = FALSE),
    align = "lcccccccc",
    escape = FALSE
  ) %>%
  kableExtra::kable_styling(
    full_width = FALSE,
    position = "center",
    fixed_thead = TRUE,
    latex_options = c("striped", "repeat_header", "HOLD_position"),
    font_size = 12
    
  ) %>%
  kableExtra::row_spec(0, bold = TRUE) %>%
  kableExtra::column_spec(1, bold = TRUE)
# NOTE: This summary table of the results for the frequentist model is reported as table 2 in the supplementary materials.


# Save summary table
openxlsx::write.xlsx(glm_honesty_humility_fit_summary, "Demonstration/Results/glm_honesty_humility_fit_summary.xlsx")





# Power analysis for simulation -------------------------------------------


# Power analysis for numerical experiment
sample_size <- 200
alpha <- 0.05
simulations <- 10000
indicators_per_latent <- 3
factor_loading <- 0.473
power_analysis <- foreach::foreach(i = 1:simulations, .combine = "rbind") %dopar% {
  require(lavaan)
  cfa_model <- paste0(
    'y1_imp =~ ', paste0("x", 1:indicators_per_latent, collapse = " + ")
  )
  while(TRUE){
    df_test <- data.frame(
      y1 = rnorm(sample_size)
    )
    for(i in 1:indicators_per_latent){
      df_test[,paste0("x", i)] <- rnorm(sample_size, factor_loading*df_test$y1)
    }
    tryCatch(
      expr = {
        cfa_fit <- lavaan::cfa(
          cfa_model,
          data = df_test,
          std.lv = FALSE,
          std.ov = FALSE
        )
      }
    )
    # Check to prevent non-nonsensical results
    if(!is.null(cfa_fit) & lavaan::lavInspect(cfa_fit, what = "converged") & lavInspect(cfa_fit, "post.check")){
      break
    }else{
      warning("\n Model did not converge")
    }
  }
  cfa_fit_summary <- lavaan::standardizedSolution(cfa_fit)
  output <- data.frame(
    false_negative = !all(cfa_fit_summary$pvalue[cfa_fit_summary$op == "=~"] <= alpha)
  )
  return(output)
}
paste0((1-mean(power_analysis$false_negative))*100, "%")
# NOTE: Statistical power = 80.38%. Reported in section 3.4 of the main paper.





# Numerical experiment ----------------------------------------------------


# CONDITION: Varying magnitude of factor loadings #
population_size <- 1000000
n_simulations <- 100
factor_loading <- c(0.473, 0.6, 0.7) # Varying conditions
smallest_sample_size <- 200
largest_sample_size <- 1000
sample_sizes <- ceiling(seq(smallest_sample_size, largest_sample_size, length.out = n_simulations)) # Sample sizes to vary
simulation_results <- data.frame(
  factor_loading = numeric(),
  sample_size = numeric(),
  fsr_bias_mean = numeric(),
  fsr_bias_sd = numeric(),
  fsr_bias_ci_low = numeric(),
  fsr_bias_ci_high = numeric(),
  bsfr_narrow_bias_mean = numeric(),
  bsfr_narrow_bias_sd = numeric(),
  bsfr_narrow_bias_ci_low = numeric(),
  bsfr_narrow_bias_ci_high = numeric(),
  bsfr_narrow_fsr_diff_mean = numeric(),
  bsfr_narrow_fsr_diff_sd = numeric(),
  bsfr_narrow_fsr_diff_ci_low = numeric(),
  bsfr_narrow_fsr_diff_ci_high = numeric(),
  bsfr_wide_bias_mean = numeric(),
  bsfr_wide_bias_sd = numeric(),
  bsfr_wide_bias_ci_low = numeric(),
  bsfr_wide_bias_ci_high = numeric(),
  bsfr_wide_fsr_diff_mean = numeric(),
  bsfr_wide_fsr_diff_sd = numeric(),
  bsfr_wide_fsr_diff_ci_low = numeric(),
  bsfr_wide_fsr_diff_ci_high = numeric(),
  bsfr_wide_narrow_diff_mean = numeric(),
  bsfr_wide_narrow_diff_sd = numeric(),
  bsfr_wide_narrow_diff_ci_low = numeric(),
  bsfr_wide_narrow_diff_ci_high = numeric(),
  fsr_mse = numeric(),
  bsfr_narrow_mse = numeric(),
  bsfr_wide_mse = numeric()
)
n_latent <- 1
indicators_per_latent <- 4
draws <- 1000
burnin <- 500
bcfa_priors_narrow <- bcfa_priors_wide <- blavaan::dpriors()
bcfa_priors_narrow["lambda"] <- "normal(0, 0.304)"
bcfa_priors_narrow["theta"] <- "gamma(12.7, 11.7)[sd]"
bcfa_priors_narrow["psi"] <- "gamma(1, 3.29)[sd]"
population_df <- data.frame(
  y1_true = rnorm(population_size)
)
latent <- "y1"

t0 <- Sys.time()
pb <- utils::txtProgressBar(min = 0, max = length(factor_loading), style = 3)
for(i in 1:length(factor_loading)){
  # Set conditions
  Lambda <- rep(factor_loading[i], indicators_per_latent)
  overt <- paste0("x", 1:indicators_per_latent)
  for(j in 1:indicators_per_latent){
    population_df[,paste0("x", j)] <- rnorm(population_size, Lambda[j]*population_df$y1_true)
  }
  population_df <- apply(population_df, 2, standardize2) %>% as.data.frame()
  cfa_model <- paste0(
    'y1_imp =~ ', paste0("x", 1:indicators_per_latent, collapse = " + ")
  )
  cfa_fit_population <- lavaan::cfa(
    cfa_model,
    data = population_df,
    std.lv = FALSE,
    std.ov = FALSE
  )
  # summary(cfa_fit_population, standardize = TRUE)
  # fsr_y1_true <- lavaan::lavInspect(cfa_fit_population, what = "fs.reliability")
  fsr_y1_true <- cor2(population_df$y1_true, lavaan::lavPredict(cfa_fit_population, method = "regression"))^2
  
  for(j in 1:n_simulations){
    # Draw sample
    sample_df <- population_df[sample(1:nrow(population_df), sample_sizes[j], replace = FALSE),]
    row.names(sample_df) <- NULL
    sample_df <- apply(sample_df, 2, standardize2) %>% as.data.frame()
    # Bayesian approach (Narrow priors)
    bcfa_fit_sample_narrow <- suppressWarnings(suppressMessages(invisible(blavaan::bcfa(
      cfa_model,
      dp = bcfa_priors_narrow,
      data = sample_df,
      std.lv = FALSE,
      std.ov = FALSE,
      save.lvs = FALSE,
      seed = SEED,
      inits = "simple",
      sample = draws,
      burnin = burnin,
      n.chains = 1,
      verbose  = FALSE,
      bcontrol = list(
        refresh = 0,
        open_progress = FALSE
      )
    ))))
    # Bayesian approach (Wide priors)
    bcfa_fit_sample_wide <- suppressWarnings(suppressMessages(invisible(blavaan::bcfa(
      cfa_model,
      dp = bcfa_priors_wide,
      data = sample_df,
      std.lv = FALSE,
      std.ov = FALSE,
      save.lvs = FALSE,
      seed = SEED,
      inits = "simple",
      sample = draws,
      burnin = burnin,
      n.chains = 1,
      verbose  = FALSE,
      bcontrol = list(
        refresh = 0,
        open_progress = FALSE
      )
    ))))
    # summary(bcfa_fit_sample_narrow, standardize = TRUE)
    # summary(bcfa_fit_sample_wide, standardize = TRUE)
    bcfa_fit_sample_estimates_narrow <- blavaan::standardizedPosterior(bcfa_fit_sample_narrow) %>% as.data.frame()
    bcfa_fit_sample_estimates_narrow <- bcfa_fit_sample_estimates_narrow %>% rename(
      "bsp_x1_miy1" = "y1_imp=~x1",
      "bsp_x2_miy1" = "y1_imp=~x2",
      "bsp_x3_miy1" = "y1_imp=~x3",
      "sigma_x1" = "x1~~x1",
      "sigma_x2" = "x2~~x2",
      "sigma_x3" = "x3~~x3",
      "sigma_y1" = "y1_imp~~y1_imp"
    )
    bfsr_narrow_est <- suppressWarnings(suppressMessages(invisible(fsr_bayes(draws = bcfa_fit_sample_estimates_narrow, latent = latent, overt = overt, prior = "no"))))
    
    bcfa_fit_sample_estimates_wide <- blavaan::standardizedPosterior(bcfa_fit_sample_wide) %>% as.data.frame()
    bcfa_fit_sample_estimates_wide <- bcfa_fit_sample_estimates_wide %>% rename(
      "bsp_x1_miy1" = "y1_imp=~x1",
      "bsp_x2_miy1" = "y1_imp=~x2",
      "bsp_x3_miy1" = "y1_imp=~x3",
      "sigma_x1" = "x1~~x1",
      "sigma_x2" = "x2~~x2",
      "sigma_x3" = "x3~~x3",
      "sigma_y1" = "y1_imp~~y1_imp"
    )
    bfsr_wide_est <- suppressWarnings(suppressMessages(invisible(fsr_bayes(draws = bcfa_fit_sample_estimates_wide, latent = latent, overt = overt, prior = "no"))))
    
    
    # Frequentist estimates
    sample_df <- apply(sample_df, 2, datawizard::standardize) %>% as.data.frame()
    frequentist_fsr_est <- foreach::foreach(i = 1:draws, .combine = "rbind") %dopar% {
      require(tidyverse)
      require(datawizard)
      require(lavaan)
      while(TRUE){
        sample_df_i <- sample_df[sample(1:nrow(sample_df), nrow(sample_df), replace = TRUE),]
        row.names(sample_df_i) <- NULL
        sample_df_i <- apply(sample_df_i, 2, datawizard::standardize) %>% as.data.frame()
        tryCatch(
          expr = {
            cfa_fit_sample_i <- lavaan::cfa(
              cfa_model,
              data = sample_df_i,
              std.lv = FALSE,
              std.ov = FALSE
            )
          }
        )
        # Check to prevent non-nonsensical results
        if(!is.null(cfa_fit_sample_i) & lavaan::lavInspect(cfa_fit_sample_i, what = "converged") & lavInspect(cfa_fit_sample_i, "post.check")){
          break
        }else{
          warning("\n Frequentist model did not converge")
        }
      }
      sample_size_i <- lavaan::lavInspect(cfa_fit_sample_i, what = "ntotal")
      npar_i <- lavaan::lavInspect(cfa_fit_sample_i, what = "npar")
      output <- data.frame(
        fsr_y1 = lavaan::lavInspect(cfa_fit_sample_i, what = "fs.reliability"),
        fsr_y1_adj = fsr_to_fsr_adj(lavaan::lavInspect(cfa_fit_sample_i, what = "fs.reliability"), n = sample_size_i, p = npar_i)
      )
      return(output)
    }
    row.names(frequentist_fsr_est) <- NULL
    
    # Compute performance measures
    fsr_bias <- frequentist_fsr_est$fsr_y1_adj - fsr_y1_true
    fsr_bias_hdi <- bayestestR::hdi(fsr_bias, ci = 1 - alpha)
    
    bsfr_narrow_bias <- bfsr_wide_est$fsr_y1 - fsr_y1_true
    bsfr_narrow_bias_hdi <- bayestestR::hdi(bsfr_narrow_bias, ci = 1 - alpha)
    
    bsfr_narrow_fsr_diff <- bfsr_wide_est$fsr_y1 - frequentist_fsr_est$fsr_y1_adj
    bsfr_narrow_fsr_diff_hdi <- bayestestR::hdi(bsfr_narrow_fsr_diff, ci = 1 - alpha)
    
    bsfr_wide_bias <- bfsr_wide_est$fsr_y1 - fsr_y1_true
    bsfr_wide_bias_hdi <- bayestestR::hdi(bsfr_wide_bias, ci = 1 - alpha)
    
    bsfr_wide_fsr_diff <- bfsr_wide_est$fsr_y1 - frequentist_fsr_est$fsr_y1_adj
    bsfr_wide_fsr_diff_hdi <- bayestestR::hdi(bsfr_wide_fsr_diff, ci = 1 - alpha)
    
    bsfr_wide_narrow_diff <- bfsr_wide_est$fsr_y1 - bfsr_narrow_est$fsr_y1
    bsfr_wide_narrow_diff_hdi <- bayestestR::hdi(bsfr_wide_narrow_diff, ci = 1 - alpha)
    
    simulation_results <- simulation_results %>% add_row(
      factor_loading = factor_loading[i],
      sample_size = sample_sizes[j],
      fsr_bias_mean = mean(fsr_bias),
      fsr_bias_sd = sd2(fsr_bias),
      fsr_bias_ci_low = fsr_bias_hdi$CI_low,
      fsr_bias_ci_high = fsr_bias_hdi$CI_high,
      bsfr_narrow_bias_mean = mean(bsfr_narrow_bias),
      bsfr_narrow_bias_sd = sd2(bsfr_narrow_bias),
      bsfr_narrow_bias_ci_low = bsfr_narrow_bias_hdi$CI_low,
      bsfr_narrow_bias_ci_high = bsfr_narrow_bias_hdi$CI_high,
      bsfr_narrow_fsr_diff_mean = mean(bsfr_narrow_fsr_diff),
      bsfr_narrow_fsr_diff_sd = sd2(bsfr_narrow_fsr_diff),
      bsfr_narrow_fsr_diff_ci_low = bsfr_narrow_fsr_diff_hdi$CI_low,
      bsfr_narrow_fsr_diff_ci_high = bsfr_narrow_fsr_diff_hdi$CI_high,
      bsfr_wide_bias_mean = mean(bsfr_wide_bias),
      bsfr_wide_bias_sd = sd2(bsfr_wide_bias),
      bsfr_wide_bias_ci_low = bsfr_wide_bias_hdi$CI_low,
      bsfr_wide_bias_ci_high = bsfr_wide_bias_hdi$CI_high,
      bsfr_wide_fsr_diff_mean = mean(bsfr_wide_fsr_diff),
      bsfr_wide_fsr_diff_sd = sd2(bsfr_wide_fsr_diff),
      bsfr_wide_fsr_diff_ci_low = bsfr_wide_fsr_diff_hdi$CI_low,
      bsfr_wide_fsr_diff_ci_high = bsfr_wide_fsr_diff_hdi$CI_high,
      bsfr_wide_narrow_diff_mean = mean(bsfr_wide_narrow_diff),
      bsfr_wide_narrow_diff_sd = sd2(bsfr_wide_narrow_diff),
      bsfr_wide_narrow_diff_ci_low = bsfr_wide_narrow_diff_hdi$CI_low,
      bsfr_wide_narrow_diff_ci_high = bsfr_wide_narrow_diff_hdi$CI_high,
      fsr_mse = mean((frequentist_fsr_est$fsr_y1_adj - fsr_y1_true)^2, na.rm = TRUE),
      bsfr_narrow_mse = mean((bfsr_narrow_est$fsr_y1 - fsr_y1_true)^2, na.rm = TRUE),
      bsfr_wide_mse = mean((bfsr_wide_est$fsr_y1 - fsr_y1_true)^2, na.rm = TRUE)
    )
  }
  utils::setTxtProgressBar(pb, i)
}
close(pb)
t1 <- Sys.time()
t1 - t0
# NOTE: Computation time was approximately 1.45 hours

write.csv(simulation_results, "Simulation/simulation_results.csv", row.names = FALSE)





# Analyze results from numerical experiment -------------------------------


# Load simulation results
simulation_results <- readr::read_csv("Simulation/simulation_results.csv") %>% as.data.frame()


# Difference between the Bayesian approach (Wide priors) and the frequentist approach
paste0(round_down2(mean(sign(simulation_results$bsfr_wide_fsr_diff_ci_low[simulation_results$factor_loading == 0.473]) != sign(simulation_results$bsfr_wide_fsr_diff_ci_high[simulation_results$factor_loading == 0.473]))*100, digits = 2), "%")
# NOTE: 100.00% of the credibility intervals overlap zero when factor loading = 0.473
paste0(round_down2(mean(sign(simulation_results$bsfr_wide_fsr_diff_ci_low[simulation_results$factor_loading == 0.6]) != sign(simulation_results$bsfr_wide_fsr_diff_ci_high[simulation_results$factor_loading == 0.6]))*100, digits = 2), "%")
# NOTE: 98.00% of the credibility intervals overlap zero when factor loading = 0.6
paste0(round_down2(mean(sign(simulation_results$bsfr_wide_fsr_diff_ci_low[simulation_results$factor_loading == 0.7]) != sign(simulation_results$bsfr_wide_fsr_diff_ci_high[simulation_results$factor_loading == 0.7]))*100, digits = 2), "%")
# NOTE: 99.00% of the credibility intervals overlap zero when factor loading = 0.7


# Difference between the Bayesian approach (Narrow priors) and the frequentist approach
paste0(round_down2(mean(sign(simulation_results$bsfr_narrow_fsr_diff_ci_low[simulation_results$factor_loading == 0.473]) != sign(simulation_results$bsfr_narrow_fsr_diff_ci_high[simulation_results$factor_loading == 0.473]))*100, digits = 2), "%")
# NOTE: 100.00% of the credibility intervals overlap zero when factor loading = 0.473
paste0(round_down2(mean(sign(simulation_results$bsfr_narrow_fsr_diff_ci_low[simulation_results$factor_loading == 0.6]) != sign(simulation_results$bsfr_narrow_fsr_diff_ci_high[simulation_results$factor_loading == 0.6]))*100, digits = 2), "%")
# NOTE: 98.00% of the credibility intervals overlap zero when factor loading = 0.6
paste0(round_down2(mean(sign(simulation_results$bsfr_narrow_fsr_diff_ci_low[simulation_results$factor_loading == 0.7]) != sign(simulation_results$bsfr_narrow_fsr_diff_ci_high[simulation_results$factor_loading == 0.7]))*100, digits = 2), "%")
# NOTE: 99.00% of the credibility intervals overlap zero when factor loading = 0.7


# Difference between the Bayesian approach (Wide priors) and the Bayesian approach (Wide priors)
paste0(round_down2(mean(sign(simulation_results$bsfr_wide_narrow_diff_ci_low[simulation_results$factor_loading == 0.473]) != sign(simulation_results$bsfr_wide_narrow_diff_ci_high[simulation_results$factor_loading == 0.473]))*100, digits = 2), "%")
# NOTE: 100.00% of the credibility intervals overlap zero when factor loading = 0.473
paste0(round_down2(mean(sign(simulation_results$bsfr_wide_narrow_diff_ci_low[simulation_results$factor_loading == 0.6]) != sign(simulation_results$bsfr_wide_narrow_diff_ci_high[simulation_results$factor_loading == 0.6]))*100, digits = 2), "%")
# NOTE: 99.00% of the credibility intervals overlap zero when factor loading = 0.6
paste0(round_down2(mean(sign(simulation_results$bsfr_wide_narrow_diff_ci_low[simulation_results$factor_loading == 0.7]) != sign(simulation_results$bsfr_wide_narrow_diff_ci_high[simulation_results$factor_loading == 0.7]))*100, digits = 2), "%")
# NOTE: 100.00% of the credibility intervals overlap zero when factor loading = 0.7


# NOTE: The above results are reported in section 3.4 of the main paper.


# Plot: RMSE Comparisons for varying priors and sample sizes
fig_width <- 7
fig_height <- fig_width/golden_ratio
plot_simulation <- data.frame(
  approach = c(
    rep("Frequentist", nrow(simulation_results)),
    rep("Bayesian (Narrow priors)", nrow(simulation_results)),
    rep("Bayesian (Wide priors)", nrow(simulation_results))
  ) %>% factor(levels = c("Frequentist", "Bayesian (Narrow priors)", "Bayesian (Wide priors)")),
  rmse = c(
    simulation_results$fsr_mse,
    simulation_results$bsfr_narrow_mse,
    simulation_results$bsfr_wide_mse
  ) %>% sqrt(),
  sample_size = c(
    simulation_results$sample_size,
    simulation_results$sample_size,
    simulation_results$sample_size
  ),
  simulation = c(
    simulation_results$factor_loading,
    simulation_results$factor_loading,
    simulation_results$factor_loading
  ) %>% as.factor()
) %>% ggplot(aes(x = sample_size, y = rmse, group = approach, color = approach, linetype = approach, shape = approach)) +
  facet_wrap(~simulation, labeller = as_labeller(c(
    "0.473" = "$\\lambda = 0.473$",
    "0.6" = "$\\lambda = 0.6$",
    "0.7" = "$\\lambda = 0.7$"
  ))) +
  geom_line(linewidth = 1, alpha = 0.8) +
  scale_x_continuous(breaks = seq(100, 1000, 150)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  labs(
    x = "Sample size",
    y = "RMSE",
    group = "Approach",
    color = "Approach",
    linetype = "Approach"
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8, linetype = 2, alpha = 0.4) +
  theme(
    axis.title.x = element_text(margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 10))
  ) +
  theme_classic(); plot_simulation
ggsave("Simulation/Plots/plot_simulation.png", plot = plot_simulation, width = fig_width, height = fig_height, units = "in", dpi = fig_dpi)
ggsave("Simulation/Plots/plot_simulation.pdf", plot = plot_simulation, width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
tikz("Simulation/Plots/plot_simulation.tikz", standAlone = FALSE, width = 7, height = 4, sanitize = FALSE)
print(plot_simulation)
dev.off()
# NOTE: This plot is reported as figure 1 of the main paper.





#################
# END OF SCRIPT #
#################