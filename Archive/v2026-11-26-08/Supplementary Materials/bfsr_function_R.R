# Script containing a simple implementation of the Bayesian factor score reliability (BFSR) index function written in R (v4.4.3 [2025-02-28 ucrt]; R Core Team, 2025), as conceptualized by the paper ’Bayesian Factor Score Reliability: Conceptualizing a Bayesian Measure of Factor Score Indeterminacy’ and detailed in its supplementary materials.

# Author: Emil Niclas Meyer-Hansen, Department of Political Science, Aarhus University (email: emil098meyerhansen@gmail.com)
# Date: 2026/11/26
# Copyright: 2026 The author. License: Except where otherwise indicated, all contents of this document and associated files are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). All data and software, including but not necessarily limited to, source code, executable code, code snippets, code chunks, algorithms, and/or scripts, attributable to this document and/or its associated files are expressly excluded from the foregoing license, and unless otherwise indicated, are instead licensed under the GNU General Public License, version 3 (GPL-3.0) but permission is granted to re-license these materials under later versions. The processed data sampled from the Open-Source Psychometrics Project (2014; provided in Petras & Dantlgraber, 2022) is excluded from the foregoing licenses and is instead licensed under the GNU General Public License, version 2 (GPL-2.0). By engaging with this document and/or any associated files, which include, but are not necessarily limited to, downloading, using, viewing, and/or distributing any of them, in parts of whole, you agree to comply with the applicable license terms for the respective content types.

# Description: The BFSR function computes the posterior factor score reliability of plausible scores using the posterior samples associated from a factor analysis model.
# This function solely supports plausible scores generated using the regression-based (‘Thomson-Thurstone’) approach and posteriors from models fit with the brms R package (v2.22.0; Bürkner, 2017, 2018).
# It relies on a series of utility functions that implement the d-wise notation described in the paper, also defined below.





# Setup -------------------------------------------------------------------

# Load (or install) required packages
# install.packages("pacman")
pacman::p_load(
  einsum,
  MASS,
  stringr,
  posterior,
  tidyverse,
  datawizard,
  parallel,
  brms
)





# Utility functions --------------------------------------------------------

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





# Bayesian factor score reliability index function ------------------------

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
    bfsr <- one1ld %*d% (W_TT_star %od% I_lld) %*d% invd(Sigma_Y %od% I_lld); dimnames(bfsr)[[2]] <- paste0("rho_", latent)
    bfsr <- posterior::as_draws_df(aperm(bfsr, c(3, 1, 2))); colnames(bfsr)[1:l] <- paste0("rho_", latent)
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
    prior_bfsr <- one1ld %*d% (prior_W_TT_star %od% I_lld) %*d% invd(prior_Sigma_Y %od% I_lld); dimnames(prior_bfsr)[[2]] <- paste0("prior_", "rho_", latent)
    prior_bfsr <- posterior::as_draws_df(aperm(prior_bfsr, c(3, 1, 2))); colnames(prior_bfsr)[1:l] <- paste0("prior_", "rho_", latent)
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





# Example -----------------------------------------------------------------

# Specify cores for MCMC
cores <- chains <- max(c(1, parallel::detectCores()-1))

# Simulate data
sample_size <- 200
eta_true <- rnorm(sample_size) # Latent variable
lambda <- 2 # Factor loading
df <- data.frame(
  m1 = lambda*eta_true + rnorm(sample_size),
  m2 = lambda*eta_true + rnorm(sample_size),
  m3 = lambda*eta_true + rnorm(sample_size),
  m4 = lambda*eta_true + rnorm(sample_size)
)
df <- apply(df, 2, datawizard::standardize) %>% as.data.frame() # Standardize for simplicity

# Specify and fit Bayesian model for factor analysis
df$eta_imp <- as.numeric(NA)
bglm_model <- brms::bf(m1 ~ 0 + mi(eta_imp)) +
  brms::bf(m2 ~ 0 + mi(eta_imp)) +
  brms::bf(m3 ~ 0 + mi(eta_imp)) +
  brms::bf(m4 ~ 0 + mi(eta_imp)) +
  brms::bf(eta_imp | mi() ~ 0) +
  brms::set_rescor(rescor = FALSE)
brms::get_prior(bglm_model, data = df, family = gaussian())
bglm_prior <- c(
  # Specify factor loadings of indicators
  brms::set_prior("normal(0, 1)", class = "b", coef = "mieta_imp", resp = "m1"),
  brms::set_prior("normal(0, 1)", class = "b", coef = "mieta_imp", resp = "m2"),
  brms::set_prior("normal(0, 1)", class = "b", coef = "mieta_imp", resp = "m3"),
  brms::set_prior("normal(0, 1)", class = "b", coef = "mieta_imp", resp = "m4"),
  
  # Specify partial dispersions of indicators
  brms::set_prior("gamma(1, 1)", class = "sigma", resp = "m1", lb = 0),
  brms::set_prior("gamma(1, 1)", class = "sigma", resp = "m2", lb = 0),
  brms::set_prior("gamma(1, 1)", class = "sigma", resp = "m3", lb = 0),
  brms::set_prior("gamma(1, 1)", class = "sigma", resp = "m4", lb = 0),
  
  # Specify dispersion of latent variable
  brms::set_prior("gamma(1, 1)", class = "sigma", resp = "etaimp", lb = 0)
)
bglm_fit <- brms::brm(
  bglm_model,
  data = df,
  prior = bglm_prior,
  family = gaussian(),
  seed = 1,
  cores = cores,
  chains = chains,
  sample_prior = "yes", # Samples prior and posterior
  iter = 2000 # NOTE: Only used for brevity of example. Larger value generally needed for these models
)
# summary(bglm_fit)

# Compute prior and posterior Bayesian factor score reliability index
bglm_bfsr <- fsr_bayes(bglm_fit, overt = paste0("m", 1:4), latent = "eta_imp", prior = "yes",  ndraws = 100)
mean(bglm_bfsr$prior_rho_eta_imp); sd(bglm_bfsr$prior_rho_eta_imp) # Prior: Mean and SD
mean(bglm_bfsr$rho_eta_imp); sd(bglm_bfsr$rho_eta_imp) # Posterior: Mean and SD


#################
# END OF SCRIPT #
#################