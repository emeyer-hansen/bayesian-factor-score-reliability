# Script written in R (v4.4.3 [2025-02-28 ucrt]; R Core Team, 2025) for validating the derivations/proofs as detailed by supplementary materials for the paper ’Bayesian Factor Score Reliability: Conceptualizing a Bayesian Measure of Factor Score Indeterminacy’.

# Author: Emil Niclas Meyer-Hansen, Department of Political Science, Aarhus University (email: emil098meyerhansen@gmail.com)
# Date: 2026/11/26
# Copyright: 2026 The author. License: Except where otherwise indicated, all contents of this document and associated files are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). All data and software, including but not necessarily limited to, source code, executable code, code snippets, code chunks, algorithms, and/or scripts, attributable to this document and/or its associated files are expressly excluded from the foregoing license, and unless otherwise indicated, are instead licensed under the GNU General Public License, version 3 (GPL-3.0) but permission is granted to re-license these materials under later versions. The processed data sampled from the Open-Source Psychometrics Project (2014; provided in Petras & Dantlgraber, 2022) is excluded from the foregoing licenses and is instead licensed under the GNU General Public License, version 2 (GPL-2.0). By engaging with this document and/or any associated files, which include, but are not necessarily limited to, downloading, using, viewing, and/or distributing any of them, in parts of whole, you agree to comply with the applicable license terms for the respective content types.


# NOTE: This script contains the R code to reproduce the numerical evidence that corroborates the mathematical derivations and proof detailed in the supplementary files for the article 'Bayesian Factor Score Reliability: Conceptualizing a Bayesian Measure of Factor Score Indeterminacy'.





# Setup -------------------------------------------------------------------

# Load required R packages
library(pacman)
pacman::p_load(
  tidyverse,
  trialr
)

# Set document specifications
digits <- 3


# Define function for validating statements that accounts for numerical imprecision
validate <- function(x = NULL, y = NULL, criterion = TRUE){
  require(dplyr)
  txt1 <- "x and y are identical within machine precision"
  txt2 <- "x and y are not identical within machine precision"
  if(all(dplyr::near(x, y)) != criterion){
    if(criterion){
      output <- paste0("FALSE. ", txt2)
    }else{
      output <- paste0("FALSE. ", txt1)
    }
    stop(output)
  } 
  if(criterion){
    output <- paste0("TRUE. ", txt1)
  }else{
    output <- paste0("TRUE. ", txt2)
  }
  return(output)
}


# Set seed for replicable RNG
SEED <- 322047077 # NOTE: Same seed as used in the main paper
set.seed(SEED)





# General setup -----------------------------------


# NOTE: For details on the domain of the specifications below, see the main paper and supplementary materials.


n <- 1000 # Number of units
o <- 12 # Number of overt phenomena
l <- 3 # Number of latent phenomena
k <- l + o # Number of total phenomena

overt <- paste0("o", 1:o) # Names of overt phenomena
latent <- paste0("l", 1:l) # Names of latent phenomena
variables <- c(latent, overt) # Names of all phenomena


# Specify arithmetic mean vector
alpha <- matrix(rnorm(k), 1, k); colnames(alpha) <- variables # Mean vector for the entire data matrix
alpha_Y <- matrix(alpha[,latent], 1, l); colnames(alpha_Y) <- latent # Mean vector for the latent data matrix
alpha_X <- matrix(alpha[,overt], 1, o); colnames(alpha_X) <- overt # Mean vector for the overt data matrix


# Specify invertible block-matrix of coefficients
while(TRUE){
  # Simulate Beta
  while(TRUE){
    mu_Beta <- matrix(rnorm(l), 1, l)
    Sigma_Beta <- trialr::rlkjcorr(1, l, 1)
    Beta <- matrix(brms::rmulti_normal(l, mu_Beta, Sigma_Beta), l, l)
    diag(Beta) <- 0
    if(abs(det(diag(1, l, l) - Beta)) > 0){
      colnames(Beta) <- latent
      rownames(Beta) <- latent
      break
    } 
  }
  # Simulate Delta
  mu_Delta <- matrix(rnorm(o), 1, o)
  Sigma_Delta <- trialr::rlkjcorr(1, o, 1)
  Delta <- matrix(brms::rmulti_normal(l, mu_Delta, Sigma_Delta), l, o); colnames(Delta) <- overt; rownames(Delta) <- latent

  # Simulate Lambda
  mu_Lambda <- matrix(rnorm(l), 1, l)
  Sigma_Lambda <- trialr::rlkjcorr(1, l, 1)
  Lambda <- matrix(brms::rmulti_normal(o, mu_Lambda, Sigma_Lambda), o, l); colnames(Lambda) <- latent; rownames(Lambda) <- overt
  
  # Simulate Omega
  while(TRUE){
    mu_Omega <- matrix(rnorm(o), 1, o)
    Sigma_Omega <- trialr::rlkjcorr(1, o, 1)
    Omega <- matrix(brms::rmulti_normal(o, mu_Omega, Sigma_Omega), o, o)
    diag(Omega) <- 0
    if(abs(det(diag(1, o, o) - Omega)) > 0){
      colnames(Omega) <- overt
      rownames(Omega) <- overt
      break
    } 
  }
  # Combine into a coefficient block-matrix
  A <- rbind(cbind(Beta, Delta), cbind(Lambda, Omega))
  if(abs(det(diag(1, k) - A)) > 0){
    A_star <- solve(diag(1, k) - A)
    break
  }
}

# Compute components of (I - A)^-1
Beta_star <- solve(diag(1, l, l) - Beta)
Omega_star <- solve(diag(1, o, o) - Omega)
H_star <- solve(diag(1, o, o) - Omega - Lambda %*% Beta_star %*% Delta)
G_star <- solve(diag(1, l, l) - Beta - Delta %*% Omega_star %*% Lambda)
H <- Beta_star %*% Delta %*% H_star
G <- Omega_star %*% Lambda %*% G_star


# Simulate partial residuals (Xi)
mu_Xi <- matrix(rnorm(k), 1, k)
Psi <- trialr::rlkjcorr(1, k, 1); colnames(Psi) <- variables; rownames(Psi) <- variables
Xi <- matrix(brms::rmulti_normal(n, mu_Xi, Psi), n, k); colnames(Xi) <- variables
Xi_Y <- matrix(Xi[,latent], n, l); colnames(Xi_Y) <- latent
Xi_X <- matrix(Xi[,overt], n, o); colnames(Xi_X) <- overt

mu_Xi <- matrix(colMeans(Xi), 1, k); colnames(mu_Xi) <- variables
mu_Xi_Y <- matrix(mu_Xi[,latent], 1, l); colnames(mu_Xi_Y) <- latent
mu_Xi_X <- matrix(mu_Xi[,overt], 1, o); colnames(mu_Xi_X) <- overt

bar_Xi <- Xi - matrix(1, n, 1) %*% mu_Xi
bar_Xi_Y <- matrix(bar_Xi[,latent], n, l); colnames(bar_Xi_Y) <- latent
bar_Xi_X <- matrix(bar_Xi[,overt], n, o); colnames(bar_Xi_X) <- overt


# Simulate data matrices
df <- (matrix(1, n, 1) %*% alpha + Xi) %*% t(A_star)
df_Y <- matrix(df[,latent], n, l); colnames(df_Y) <- latent
df_X <- matrix(df[,overt], n, o); colnames(df_X) <- overt

mu <- matrix(colMeans(df), 1, k); colnames(mu) <- variables
mu_Y <- matrix(mu[,latent], 1, l); colnames(mu_Y) <- latent
mu_X <- matrix(mu[,overt], 1, o); colnames(mu_X) <- overt

bar_df <- df - matrix(1, n, 1) %*% mu
bar_Y <- matrix(bar_df[,latent], n, l); colnames(bar_Y) <- latent
bar_X <- matrix(bar_df[,overt], n, o); colnames(bar_X) <- overt


# Get residual matrices (E)
E <- df - matrix(1, n, 1) %*% mu
E_Y <- matrix(E[,latent], n, l); colnames(E_Y) <- latent
E_X <- matrix(E[,overt], n, o); colnames(E_X) <- overt

mu_E <- matrix(colMeans(E), 1, k); colnames(mu_E) <- variables
mu_E_Y <- matrix(mu_E[,latent], 1, l); colnames(mu_E_Y) <- latent
mu_E_X <- matrix(mu_E[,overt], 1, o); colnames(mu_E_X) <- overt

bar_E <- E - matrix(1, n, 1) %*% mu_E
bar_E_Y <- matrix(bar_E[,latent], n, l); colnames(bar_E_Y) <- latent
bar_E_X <- matrix(bar_E[,overt], n, o); colnames(bar_E_X) <- overt


# Get variance-covariance matrices
Sigma <- 1/n * t(bar_E) %*% bar_E
Sigma_Y <- matrix(Sigma[latent, latent], l, l); colnames(Sigma_Y) <- latent; rownames(Sigma_Y) <- latent
Sigma_X <- matrix(Sigma[overt, overt], o, o); colnames(Sigma_X) <- overt; rownames(Sigma_X) <- overt
Sigma_Y_X <- matrix(Sigma[latent, overt], l, o); colnames(Sigma_Y_X) <- overt; rownames(Sigma_Y_X) <- latent
Sigma_X_Y <- matrix(Sigma[overt, latent], o, l); colnames(Sigma_X_Y) <- latent; rownames(Sigma_X_Y) <- overt


# Get partial variance-covariance matrices
Psi <- 1/n * t(bar_Xi) %*% bar_Xi
Psi_Y <- matrix(Psi[latent, latent], l, l); colnames(Psi_Y) <- latent; rownames(Psi_Y) <- latent
Psi_X <- matrix(Psi[overt, overt], o, o); colnames(Psi_X) <- overt; rownames(Psi_X) <- overt
Psi_Y_X <- matrix(Psi[latent, overt], l, o); colnames(Psi_Y_X) <- overt; rownames(Psi_Y_X) <- latent
Psi_X_Y <- matrix(Psi[overt, latent], o, l); colnames(Psi_X_Y) <- latent; rownames(Psi_X_Y) <- overt


# Define transposition and diagonalization matrices
T_matrix <- rbind(cbind(matrix(0, l, o), diag(1, l, l)), cbind(diag(1, o, o), matrix(0, o, l)))
D_1 <- rbind(cbind(matrix(1, l, o), matrix(0, l, l)), cbind(matrix(0, o, o), matrix(1, o, l)))
D_2 <- rbind(cbind(matrix(1, o, o), matrix(0, o, l)), cbind(matrix(0, l, o), matrix(1, l, l)))


# Define scoring weights
W_tt <- (Sigma %*% T_matrix * D_1) %*% solve(t(T_matrix) %*% Sigma %*% T_matrix * D_2); rownames(W_tt) <- c(latent, overt); colnames(W_tt) <- c(overt, latent)
W_tt_Y <- matrix(W_tt[latent, overt], l, o); rownames(W_tt_Y) <- latent; colnames(W_tt_Y) <- overt
W_tt_X <- matrix(W_tt[overt, latent], o, l); rownames(W_tt_X) <- overt; colnames(W_tt_X) <- latent


# Define imputed values
df_imp <- matrix(1, n, 1) %*% mu + (df - matrix(1, n, 1) %*% mu) %*% T_matrix %*% t(W_tt)
df_imp_Y <- df_imp[,latent]
df_imp_X <- df_imp[,overt]


# Get adjusted factor scoring weights
W_tt_star <- W_tt %*% t(Sigma %*% T_matrix * D_1)
W_tt_star_Y <- matrix(W_tt_star[latent, latent], l, l); rownames(W_tt_star_Y) <- latent; colnames(W_tt_star_Y) <- latent
W_tt_star_X <- matrix(W_tt_star[overt, overt], o, o); rownames(W_tt_star_X) <- overt; colnames(W_tt_star_X) <- overt


# Define FSR indeX
rho_sq <- matrix(1, 1, k) %*% (W_tt_star * diag(1, k, k)) %*% solve(Sigma * diag(1, k, k))
rho_sq_Y <- matrix(rho_sq[,latent], 1, l); colnames(rho_sq_Y) <- latent
rho_sq_X <- matrix(rho_sq[,overt], 1, o); colnames(rho_sq_X) <- overt


# NOTE: Henceforth, only changes from these specifications will be denoted





# Generalized fundamental equations -----------------------------------------------

# Validate eqs. 1-2
validate(df_Y, matrix(1, n, 1) %*% alpha_Y + df_Y %*% t(Beta) + df_X %*% t(Delta) + Xi_Y) # TRUE
validate(df_X, matrix(1, n, 1) %*% alpha_X + df_Y %*% t(Lambda) + df_X %*% t(Omega) + Xi_X) # TRUE

# Validate eqs. 3-4
validate(cbind(df_Y, df_X), matrix(1, n, 1) %*% cbind(alpha_Y, alpha_X) + cbind(df_Y, df_X) %*% t(rbind(cbind(Beta, Delta), cbind(Lambda, Omega))) + cbind(Xi_Y, Xi_X)) # TRUE
validate(cbind(df_Y, df_X), df) # TRUE
validate(df, matrix(1, n, 1) %*% alpha + df %*% t(A) + Xi) # TRUE

# Validate eqs. 5-6
validate(df - df %*% t(A), matrix(1, n, 1) %*% alpha + Xi) # TRUE
validate(df - df %*% t(A), df %*% (diag(1, k, k) - t(A))) # TRUE

# Validate eq. 7
validate(diag(1, k, k) - A, rbind(cbind(diag(1, l, l) - Beta, -Delta), cbind(-Lambda, diag(1, o, o) - Omega))) # TRUE

# Validate eq. 8
validate(A_star, rbind(cbind(G_star, H), cbind(G, H_star))) # TRUE

# Validate eq. 9
validate(df, matrix(1, n, 1) %*% alpha %*% t(A_star) + Xi %*% t(A_star)) # TRUE

# Validate eqs. 10-12
validate(colMeans(df), colMeans(matrix(1, n, 1) %*% alpha %*% t(A_star) + Xi %*% t(A_star))) # TRUE
validate(colMeans(df), mu) # TRUE
validate(mu, alpha %*% t(A_star) + mu_Xi %*% t(A_star)) # TRUE
validate(mu, (alpha + mu_Xi) %*% t(A_star)) # TRUE

# Validate eqs. 13-15
validate(mu, cbind(mu_Y, mu_X)) # TRUE
validate(cbind(mu_Y, mu_X), (cbind(alpha_Y, alpha_X) + cbind(mu_Xi_Y, mu_Xi_X)) %*% rbind(cbind(t(G_star), t(G)), cbind(t(H), t(H_star)))) # TRUE
validate(mu_Y, (alpha_Y + mu_Xi_Y) %*% t(G_star) + (alpha_X + mu_Xi_X) %*% t(H)) # TRUE
validate(mu_X, (alpha_Y + mu_Xi_Y) %*% t(G) + (alpha_X + mu_Xi_X) %*% t(H_star)) # TRUE

# Validate eqs. 16-17
validate(df - matrix(1, n, 1) %*% mu, Xi %*% t(A_star) - matrix(1, n, 1) %*% mu_Xi %*% t(A_star)) # TRUE
validate(df - matrix(1, n, 1) %*% mu, E) # TRUE
validate(E, (Xi - matrix(1, n, 1) %*% mu_Xi) %*% t(A_star)) # TRUE

# Validate eq. 18
validate(Sigma, rbind(cbind(Sigma_Y, Sigma_Y_X), cbind(Sigma_X_Y, Sigma_X))) # TRUE

# Validate eq. 19
validate(Psi, rbind(cbind(Psi_Y, Psi_Y_X), cbind(Psi_X_Y, Psi_X))) # TRUE

# Validate eqs. 20-29
validate(Sigma, 1/n * t(bar_E) %*% bar_E) # TRUE
validate(Sigma, 1/n * t(bar_Xi %*% t(A_star)) %*% (bar_Xi %*% t(A_star))) # TRUE
validate(Sigma, 1/n * (A_star %*% t(bar_Xi))  %*% (bar_Xi %*% t(A_star))) # TRUE
validate(Sigma, A_star %*% (1/n * t(bar_Xi) %*% bar_Xi) %*% t(A_star)) # TRUE
validate(Sigma, A_star %*% Psi %*% t(A_star)) # TRUE
validate(Sigma, rbind(cbind(Sigma_Y, Sigma_Y_X), cbind(Sigma_X_Y, Sigma_X))) # TRUE
validate(rbind(cbind(Sigma_Y, Sigma_Y_X), cbind(Sigma_X_Y, Sigma_X)), rbind(cbind(G_star, H), cbind(G, H_star)) %*% rbind(cbind(Psi_Y, Psi_Y_X), cbind(Psi_X_Y, Psi_X)) %*% rbind(cbind(t(G_star), t(G)), cbind(t(H), t(H_star)))) # TRUE
validate(Sigma_Y, G_star %*% Psi_Y %*% t(G_star) + H %*% Psi_X %*% t(H) + G_star %*% Psi_Y_X %*% t(H) + H %*% Psi_X_Y %*% t(G_star)) # TRUE
validate(Sigma_X, G %*% Psi_Y %*% t(G) + H_star %*% Psi_X %*% t(H_star) + G %*% Psi_Y_X %*% t(H_star) + H_star %*% Psi_X_Y %*% t(G)) # TRUE
validate(Sigma_Y_X, G_star %*% Psi_Y %*% t(G) + H %*% Psi_X %*% t(H_star) + G_star %*% Psi_Y_X %*% t(H_star) + H %*% Psi_X_Y %*% t(G)) # TRUE
validate(Sigma_X_Y, G %*% Psi_Y %*% t(G_star) + H_star %*% Psi_X %*% t(H) + G %*% Psi_Y_X %*% t(H) + H_star %*% Psi_X_Y %*% t(G_star)) # TRUE


# End of derivation/proof #





# 1.2 Factor scores -------------------------------------------------------

# Validate eqs. 30-37
validate(
  W_tt,
  (Sigma %*% T_matrix * D_1) %*%
    solve(
      t(T_matrix) %*% Sigma %*% T_matrix * D_2
    )
) # TRUE
validate(
  W_tt,
  (rbind(cbind(Sigma_Y, Sigma_Y_X), cbind(Sigma_X_Y, Sigma_X)) %*% rbind(cbind(matrix(0, l, o), diag(1, l, l)), cbind(diag(1, o, o), matrix(0, o, l))) * rbind(cbind(matrix(1, l, o), matrix(0, l, l)), cbind(matrix(0, o, o), matrix(1, o, l)))) %*%
    solve(
      rbind(cbind(matrix(0, o, l), diag(1, o, o)), cbind(diag(1, l, l), matrix(0, l, o))) %*% rbind(cbind(Sigma_Y, Sigma_Y_X), cbind(Sigma_X_Y, Sigma_X)) %*% rbind(cbind(matrix(0, l, o), diag(1, l, l)), cbind(diag(1, o, o), matrix(0, o, l))) * rbind(cbind(matrix(1, o, o), matrix(0, o, l)), cbind(matrix(0, l, o), matrix(1, l, l)))
    )
) # TRUE
validate(
  W_tt, rbind(cbind(W_tt_Y, matrix(0, l, l)), cbind(matrix(0, o, o), W_tt_X))
) # TRUE
validate(
  W_tt, rbind(cbind(Sigma_Y_X, matrix(0, l, l)), cbind(matrix(0, o, o), Sigma_X_Y)) %*% rbind(cbind(solve(Sigma_X), matrix(0, o, l)), cbind(matrix(0, l, o), solve(Sigma_Y)))
) # TRUE
validate(
  W_tt, rbind(cbind(Sigma_Y_X %*% solve(Sigma_X), matrix(0, l, l)), cbind(matrix(0, o, o), Sigma_X_Y %*% solve(Sigma_Y)))
) # TRUE
validate(W_tt_Y, Sigma_Y_X %*% solve(Sigma_X)) # TRUE
validate(W_tt_Y, t(bar_Y) %*% bar_X %*% solve(t(bar_X) %*% bar_X)) # TRUE
validate(W_tt_X, Sigma_X_Y %*% solve(Sigma_Y)) # TRUE
validate(W_tt_X, t(bar_X) %*% bar_Y %*% solve(t(bar_Y) %*% bar_Y)) # TRUE

# Validate eqs. 38-41
validate(df_imp, matrix(1, n, 1) %*% mu + (df - matrix(1, n, 1) %*% mu) %*% T_matrix %*% t(W_tt)) # TRUE
validate(df_imp, cbind(df_imp_Y, df_imp_X)) # TRUE
validate(
  cbind(df_imp_Y, df_imp_X),
  matrix(1, n, 1) %*% cbind(mu_Y, mu_X) + (
    cbind(df_Y, df_X) - matrix(1, n, 1) %*% cbind(mu_Y, mu_X)
  )
  %*%
    rbind(cbind(matrix(0, l, o), diag(1, l, l)), cbind(diag(1, o, o), matrix(0, o, l)))
  %*%
    rbind(cbind(t(Sigma_Y_X %*% solve(Sigma_X)), matrix(0, o, o)), cbind(matrix(0, l, l), t(Sigma_X_Y %*% solve(Sigma_Y))))
) # TRUE
validate(df_imp_Y, matrix(1, n, 1) %*% mu_Y + (df_X - matrix(1, n, 1) %*% mu_X) %*% t(W_tt_Y)) # TRUE
validate(df_imp_X, matrix(1, n, 1) %*% mu_X + (df_Y - matrix(1, n, 1) %*% mu_Y) %*% t(W_tt_X)) # TRUE


# Validate that the imputed scores do not equal the true scores
validate(df, df_imp, FALSE) # TRUE


# End of derivation/proof #





# 1.3 Factor score reliability ------------------------------------------------

# Validate eqs. 42-47
validate(
  W_tt_star,
  W_tt %*% t(Sigma %*% T_matrix * D_1)
) # TRUE
validate(
  W_tt_star,
  rbind(cbind(W_tt_star_Y, matrix(0, l, o)), cbind(matrix(0, o, l), W_tt_star_X))
) # TRUE
validate(
  W_tt_star,
  rbind(
    cbind(Sigma_Y_X %*% solve(Sigma_X), matrix(0, l, l)),
    cbind(matrix(0, o, o), Sigma_X_Y %*% solve(Sigma_Y))
  ) %*%
    rbind(
      cbind(t(Sigma_Y_X), matrix(0, o, o)),
      cbind(matrix(0, l, l), t(Sigma_X_Y))
    )
) # TRUE
validate(W_tt_star_Y, W_tt_Y %*% t(Sigma_Y_X)) # TRUE
validate(W_tt_star_Y, W_tt_Y %*% Sigma_X_Y) # TRUE
validate(W_tt_star_X, W_tt_X %*% t(Sigma_X_Y)) # TRUE
validate(W_tt_star_X, W_tt_X %*% Sigma_Y_X) # TRUE


# Validate eqs. 48-51
validate(rho_sq, matrix(1, 1, k) %*% (W_tt_star * diag(1, k, k)) %*% solve(Sigma * diag(1, k, k)))
validate(rho_sq, cbind(rho_sq_Y, rho_sq_X)) # TRUE
validate(rho_sq,
         matrix(1, 1, k) %*% (
           rbind(cbind(W_tt_star_Y, matrix(0, l, o)), cbind(matrix(0, o, l), W_tt_star_X)) * diag(1, k, k)
          ) %*%
           solve(
             rbind(cbind(Sigma_Y, Sigma_Y_X), cbind(Sigma_X_Y, Sigma_X)) * diag(1, k, k)
           )
) # TRUE
validate(rho_sq_Y, matrix(1, 1, l) %*% (W_tt_star_Y * diag(1, l, l)) %*% solve(Sigma_Y * diag(1, l, l))) # TRUE
validate(rho_sq_X, matrix(1, 1, o) %*% (W_tt_star_X * diag(1, o, o)) %*% solve(Sigma_X * diag(1, o, o))) # TRUE


# End of derivation/proof #





# Factor score indeterminacy ----------------------------------------------


# First special exception

# Show that the imputed scores does not equal the true scores
validate(df, df_imp, FALSE) # TRUE

# Replace the 'T_matrix %*% t(W_tt)'-term with an identity matrix
df_imp <- matrix(1, n, 1) %*% mu + (df - matrix(1, n, 1) %*% mu) %*% diag(1, o+l, l+o)

# Show that the imputed scores now equal the true scores
validate(df, df_imp) # TRUE



# Second special exception
n <- 1000
o <- 12
l <- 3
k <- l + o

overt <- paste0("o", 1:o)
latent <- paste0("l", 1:l)
variables <- c(latent, overt)

alpha <- matrix(rnorm(k), 1, k); colnames(alpha) <- variables
alpha_Y <- matrix(alpha[,latent], 1, l); colnames(alpha_Y) <- latent
alpha_X <- matrix(alpha[,overt], 1, o); colnames(alpha_X) <- overt

while(TRUE){
  # Simulate Beta
  while(TRUE){
    mu_Beta <- matrix(rnorm(l), 1, l)
    Sigma_Beta <- trialr::rlkjcorr(1, l, 1)
    Beta <- matrix(brms::rmulti_normal(l, mu_Beta, Sigma_Beta), l, l)
    diag(Beta) <- 0
    if(abs(det(diag(1, l, l) - Beta)) > 0){
      colnames(Beta) <- latent
      rownames(Beta) <- latent
      break
    } 
  }
  # Simulate Delta
  mu_Delta <- matrix(rnorm(o), 1, o)
  Sigma_Delta <- trialr::rlkjcorr(1, o, 1)
  Delta <- matrix(brms::rmulti_normal(l, mu_Delta, Sigma_Delta), l, o); colnames(Delta) <- overt; rownames(Delta) <- latent
  
  # Simulate Lambda
  mu_Lambda <- matrix(rnorm(l), 1, l)
  Sigma_Lambda <- trialr::rlkjcorr(1, l, 1)
  Lambda <- matrix(brms::rmulti_normal(o, mu_Lambda, Sigma_Lambda), o, l); colnames(Lambda) <- latent; rownames(Lambda) <- overt
  
  # Simulate Omega
  while(TRUE){
    mu_Omega <- matrix(rnorm(o), 1, o)
    Sigma_Omega <- trialr::rlkjcorr(1, o, 1)
    Omega <- matrix(brms::rmulti_normal(o, mu_Omega, Sigma_Omega), o, o)
    diag(Omega) <- 0
    if(abs(det(diag(1, o, o) - Omega)) > 0){
      colnames(Omega) <- overt
      rownames(Omega) <- overt
      break
    } 
  }
  A <- rbind(cbind(Beta, Delta), cbind(Lambda, Omega))
  if(abs(det(diag(1, k) - A)) > 0){
    A_star <- solve(diag(1, k) - A)
    break
  }
}

Beta_star <- solve(diag(1, l, l) - Beta)
Omega_star <- solve(diag(1, o, o) - Omega)
H_star <- solve(diag(1, o, o) - Omega - Lambda %*% Beta_star %*% Delta)
G_star <- solve(diag(1, l, l) - Beta - Delta %*% Omega_star %*% Lambda)
H <- Beta_star %*% Delta %*% H_star
G <- Omega_star %*% Lambda %*% G_star

mu_Xi <- matrix(rnorm(k), 1, k) # NOTE: Unused
Psi <- trialr::rlkjcorr(1, k, 1); colnames(Psi) <- variables; rownames(Psi) <- variables # NOTE: Unused
Xi <- matrix(0, n, k); colnames(Xi) <- variables # NOTE: Has been changed to a null matrix
Xi_Y <- matrix(Xi[,latent], n, l); colnames(Xi_Y) <- latent
Xi_X <- matrix(Xi[,overt], n, o); colnames(Xi_X) <- overt

mu_Xi <- matrix(colMeans(Xi), 1, k); colnames(mu_Xi) <- variables
mu_Xi_Y <- matrix(mu_Xi[,latent], 1, l); colnames(mu_Xi_Y) <- latent
mu_Xi_X <- matrix(mu_Xi[,overt], 1, o); colnames(mu_Xi_X) <- overt

bar_Xi <- Xi - matrix(1, n, 1) %*% mu_Xi
bar_Xi_Y <- matrix(bar_Xi[,latent], n, l); colnames(bar_Xi_Y) <- latent
bar_Xi_X <- matrix(bar_Xi[,overt], n, o); colnames(bar_Xi_X) <- overt

df <- (matrix(1, n, 1) %*% alpha + Xi) %*% t(A_star)
df_Y <- matrix(df[,latent], n, l); colnames(df_Y) <- latent
df_X <- matrix(df[,overt], n, o); colnames(df_X) <- overt

mu <- matrix(colMeans(df), 1, k); colnames(mu) <- variables
mu_Y <- matrix(mu[,latent], 1, l); colnames(mu_Y) <- latent
mu_X <- matrix(mu[,overt], 1, o); colnames(mu_X) <- overt

bar_df <- df - matrix(1, n, 1) %*% mu
bar_Y <- matrix(bar_df[,latent], n, l); colnames(bar_Y) <- latent
bar_X <- matrix(bar_df[,overt], n, o); colnames(bar_X) <- overt

E <- df - matrix(1, n, 1) %*% mu
E_Y <- matrix(E[,latent], n, l); colnames(E_Y) <- latent
E_X <- matrix(E[,overt], n, o); colnames(E_X) <- overt

mu_E <- matrix(colMeans(E), 1, k); colnames(mu_E) <- variables
mu_E_Y <- matrix(mu_E[,latent], 1, l); colnames(mu_E_Y) <- latent
mu_E_X <- matrix(mu_E[,overt], 1, o); colnames(mu_E_X) <- overt

bar_E <- E - matrix(1, n, 1) %*% mu_E
bar_E_Y <- matrix(bar_E[,latent], n, l); colnames(bar_E_Y) <- latent
bar_E_X <- matrix(bar_E[,overt], n, o); colnames(bar_E_X) <- overt

Sigma <- 1/n * t(bar_E) %*% bar_E
Sigma_Y <- matrix(Sigma[latent, latent], l, l); colnames(Sigma_Y) <- latent; rownames(Sigma_Y) <- latent
Sigma_X <- matrix(Sigma[overt, overt], o, o); colnames(Sigma_X) <- overt; rownames(Sigma_X) <- overt
Sigma_Y_X <- matrix(Sigma[latent, overt], l, o); colnames(Sigma_Y_X) <- overt; rownames(Sigma_Y_X) <- latent
Sigma_X_Y <- matrix(Sigma[overt, latent], o, l); colnames(Sigma_X_Y) <- latent; rownames(Sigma_X_Y) <- overt

Psi <- 1/n * t(bar_Xi) %*% bar_Xi
Psi_Y <- matrix(Psi[latent, latent], l, l); colnames(Psi_Y) <- latent; rownames(Psi_Y) <- latent
Psi_X <- matrix(Psi[overt, overt], o, o); colnames(Psi_X) <- overt; rownames(Psi_X) <- overt
Psi_Y_X <- matrix(Psi[latent, overt], l, o); colnames(Psi_Y_X) <- overt; rownames(Psi_Y_X) <- latent
Psi_X_Y <- matrix(Psi[overt, latent], o, l); colnames(Psi_X_Y) <- latent; rownames(Psi_X_Y) <- overt

T_matrix <- rbind(cbind(matrix(0, l, o), diag(1, l, l)), cbind(diag(1, o, o), matrix(0, o, l)))
D_1 <- rbind(cbind(matrix(1, l, o), matrix(0, l, l)), cbind(matrix(0, o, o), matrix(1, o, l)))
D_2 <- rbind(cbind(matrix(1, o, o), matrix(0, o, l)), cbind(matrix(0, l, o), matrix(1, l, l)))

W_tt <- (Sigma %*% T_matrix * D_1) %*% MASS::ginv(t(T_matrix) %*% Sigma %*% T_matrix * D_2); rownames(W_tt) <- c(latent, overt); colnames(W_tt) <- c(overt, latent) # NOTE: Has been changed to invoke the generalized inverse
W_tt_Y <- matrix(W_tt[latent, overt], l, o); rownames(W_tt_Y) <- latent; colnames(W_tt_Y) <- overt
W_tt_X <- matrix(W_tt[overt, latent], o, l); rownames(W_tt_X) <- overt; colnames(W_tt_X) <- latent

df_imp <- matrix(1, n, 1) %*% mu + (df - matrix(1, n, 1) %*% mu) %*% T_matrix %*% t(W_tt)
df_imp_Y <- df_imp[,latent]
df_imp_X <- df_imp[,overt]


# Show that the true scores and imputed scores are now identical
validate(df, df_imp) # TRUE
validate(df, matrix(1, n, 1) %*% mu) # TRUE


# End of derivation/proof #





# Classical Fundamental Equations/Theorems of Factor Analysis -------------


n <- 1000
o <- 12
l_c <- 3 # NOTE: New, positive finite
l_u <- o # NOTE: New, must equal o
l <- l_c + l_u # NOTE: Updated, must sum to l
k <- l + o

overt <- paste0("o", 1:o)
latent_c <- paste0("l_c", 1:l_c) # NOTE: New
latent_u <- paste0("l_u", 1:l_u) # NOTE: New
latent <- c(latent_c, latent_u) # NOTE: Updated
variables <- c(latent, overt)

alpha <- matrix(rnorm(k), 1, k); colnames(alpha) <- variables
alpha_Y <- matrix(alpha[,latent], 1, l); colnames(alpha_Y) <- latent
alpha_Y_c <- matrix(alpha_Y[,latent_c], 1, l_c); colnames(alpha_Y_c) <- latent_c # NOTE: New
alpha_Y_u <- matrix(alpha_Y[,latent_u], 1, l_u); colnames(alpha_Y_u) <- latent_u # NOTE: New
alpha_X <- matrix(0, 1, o); colnames(alpha_X) <- overt # NOTE: alpha_X is assumed to be a 1 by o null matrix
alpha <- cbind(alpha_Y, alpha_X) # NOTE: alpha is updated

while(TRUE){
  # Simulate Beta
  while(TRUE){
    mu_Beta <- matrix(rnorm(l), 1, l) # NOTE: Unused
    Sigma_Beta <- trialr::rlkjcorr(1, l, 1) # NOTE: Unused
    Beta <- matrix(0, l, l) # NOTE: Beta is assumed to be an l by l null matrix
    diag(Beta) <- 0
    if(abs(det(diag(1, l, l) - Beta)) > 0){
      colnames(Beta) <- latent
      rownames(Beta) <- latent
      break
    } 
  }
  # Simulate Delta
  mu_Delta <- matrix(rnorm(o), 1, o) # NOTE: Unused
  Sigma_Delta <- trialr::rlkjcorr(1, o, 1) # NOTE: Unused
  Delta <- matrix(0, l, o); colnames(Delta) <- overt; rownames(Delta) <- latent # NOTE: Changed to an l by o null matrix

  
  # Simulate Lambda
  mu_Lambda <- matrix(rnorm(l), 1, l)
  Sigma_Lambda <- trialr::rlkjcorr(1, l, 1)
  Lambda <- matrix(brms::rmulti_normal(o, mu_Lambda, Sigma_Lambda), o, l); colnames(Lambda) <- latent; rownames(Lambda) <- overt
  Lambda_c <- matrix(Lambda[overt, latent_c], o, l_c); colnames(Lambda_c) <- latent_c; rownames(Lambda_c) <- overt # NOTE: New
  Lambda_u <- diag(diag(Lambda[overt, latent_u]), o, l_u); colnames(Lambda_u) <- latent_u; rownames(Lambda_u) <- overt # NOTE: Lambda_u is conceptualized to be diagonal
  Lambda <- cbind(Lambda_c, Lambda_u) # NOTE: Update Lambda
  
  # Simulate Omega
  while(TRUE){
    mu_Omega <- matrix(rnorm(o), 1, o) # NOTE: Unused
    Sigma_Omega <- trialr::rlkjcorr(1, o, 1) # NOTE: Unused
    Omega <- matrix(0, o, o)   # NOTE: Changed to an o by o null matrix
    diag(Omega) <- 0
    if(abs(det(diag(1, o, o) - Omega)) > 0){
      colnames(Omega) <- overt
      rownames(Omega) <- overt
      break
    } 
  }
  A <- rbind(cbind(Beta, Delta), cbind(Lambda, Omega))
  if(abs(det(diag(1, k) - A)) > 0){
    A_star <- solve(diag(1, k) - A)
    break
  }
}

Beta_star <- solve(diag(1, l, l) - Beta)
Omega_star <- solve(diag(1, o, o) - Omega)
H_star <- solve(diag(1, o, o) - Omega - Lambda %*% Beta_star %*% Delta)
G_star <- solve(diag(1, l, l) - Beta - Delta %*% Omega_star %*% Lambda)
H <- Beta_star %*% Delta %*% H_star
G <- Omega_star %*% Lambda %*% G_star

mu_Xi <- matrix(rnorm(k), 1, k)
Psi <- trialr::rlkjcorr(1, k, 1); colnames(Psi) <- variables; rownames(Psi) <- variables # NOTE: Unused
Xi <- matrix(rnorm(n*k, mu_Xi, diag(Psi)), n, k); colnames(Xi) <- variables # NOTE: Xi_Y (incl. Xi_Y_u and Xi_Y_u) and Xi_Y are assumed linearly independent
Xi[,latent_u] <- qr.resid(qr(cbind(1, Xi[,latent_c])), Xi[,latent_u]) # NOTE: Assume unit variance for Xi_Y_u and help enforce linear independence between Xi_Y_c and Xi_Y_u
Xi[,latent_u] <- Xi[,latent_u] %*% solve(chol(1/n * t(Xi[,latent_u]) %*% Xi[,latent_u])) # Help enforce Psi_Y_u to be an identity matrix
Xi[,overt] <- 0  # NOTE: Xi_X is assumed a null matrix
Xi_Y <- matrix(Xi[,latent], n, l); colnames(Xi_Y) <- latent
Xi_Y_c <- matrix(Xi_Y[,latent_c], n, l_c); colnames(Xi_Y_c) <- latent_c
Xi_Y_u <- matrix(Xi_Y[,latent_u], n, l_u); colnames(Xi_Y_u) <- latent_u
Xi_X <- matrix(Xi[,overt], n, o); colnames(Xi_X) <- overt

mu_Xi <- matrix(colMeans(Xi), 1, k); colnames(mu_Xi) <- variables
mu_Xi_Y <- matrix(mu_Xi[,latent], 1, l); colnames(mu_Xi_Y) <- latent
mu_Xi_Y_c <- matrix(mu_Xi_Y[,latent_c], 1, l_c); colnames(mu_Xi_Y_c) <- latent_c
mu_Xi_Y_u <- matrix(mu_Xi[,latent_u], 1, l_u); colnames(mu_Xi_Y_u) <- latent_u
mu_Xi_X <- matrix(mu_Xi[,overt], 1, o); colnames(mu_Xi_X) <- overt

bar_Xi <- Xi - matrix(1, n, 1) %*% mu_Xi
bar_Xi_Y <- matrix(bar_Xi[,latent], n, l); colnames(bar_Xi_Y) <- latent
bar_Xi_Y_c <- matrix(bar_Xi_Y[,latent_c], n, l_c); colnames(bar_Xi_Y_c) <- latent_c
bar_Xi_Y_u <- matrix(bar_Xi_Y[,latent_u], n, l_u); colnames(bar_Xi_Y_u) <- latent_u
bar_Xi_X <- matrix(bar_Xi[,overt], n, o); colnames(bar_Xi_X) <- overt

df <- (matrix(1, n, 1) %*% alpha + Xi) %*% t(A_star)
df_Y <- matrix(df[,latent], n, l); colnames(df_Y) <- latent
df_Y_c <- matrix(df_Y[,latent_c], n, l_c); colnames(df_Y_c) <- latent_c
df_Y_u <- matrix(df_Y[,latent_u], n, l_u); colnames(df_Y_u) <- latent_u
df_X <- matrix(df[,overt], n, o); colnames(df_X) <- overt

mu <- matrix(colMeans(df), 1, k); colnames(mu) <- variables
mu_Y <- matrix(mu[,latent], 1, l); colnames(mu_Y) <- latent
mu_Y_c <- matrix(mu_Y[,latent_c], 1, l_c); colnames(mu_Y_c) <- latent_c
mu_Y_u <- matrix(mu_Y[,latent_u], 1, l_u); colnames(mu_Y_u) <- latent_u
mu_X <- matrix(mu[,overt], 1, o); colnames(mu_X) <- overt

bar_df <- df - matrix(1, n, 1) %*% mu
bar_Y <- matrix(bar_df[,latent], n, l); colnames(bar_Y) <- latent
bar_Y_c <- matrix(bar_Y[,latent_c], n, l_c); colnames(bar_Y_c) <- latent_c
bar_Y_u <- matrix(bar_Y[,latent_u], n, l_u); colnames(bar_Y_u) <- latent_u
bar_X <- matrix(bar_df[,overt], n, o); colnames(bar_X) <- overt

E <- df - matrix(1, n, 1) %*% mu
E_Y <- matrix(E[,latent], n, l); colnames(E_Y) <- latent
E_Y_c <- matrix(E_Y[,latent_c], n, l_c); colnames(E_Y_c) <- latent_c
E_Y_u <- matrix(E_Y[,latent_u], n, l_u); colnames(E_Y_u) <- latent_u
E_X <- matrix(E[,overt], n, o); colnames(E_X) <- overt

mu_E <- matrix(colMeans(E), 1, k); colnames(mu_E) <- variables
mu_E_Y <- matrix(mu_E[,latent], 1, l); colnames(mu_E_Y) <- latent
mu_E_Y_c <- matrix(mu_E_Y[,latent_c], 1, l_c); colnames(mu_E_Y_c) <- latent_c
mu_E_Y_u <- matrix(mu_E_Y[,latent_u], 1, l_u); colnames(mu_E_Y_u) <- latent_u
mu_E_X <- matrix(mu_E[,overt], 1, o); colnames(mu_E_X) <- overt

bar_E <- E - matrix(1, n, 1) %*% mu_E
bar_E_Y <- matrix(bar_E[,latent], n, l); colnames(bar_E_Y) <- latent
bar_E_Y_c <- matrix(bar_E_Y[,latent_c], n, l_c); colnames(bar_E_Y_c) <- latent_c
bar_E_Y_u <- matrix(bar_E_Y[,latent_u], n, l_u); colnames(bar_E_Y_u) <- latent_u
bar_E_X <- matrix(bar_E[,overt], n, o); colnames(bar_E_X) <- overt

Sigma <- 1/n * t(bar_E) %*% bar_E
Sigma_Y <- matrix(Sigma[latent, latent], l, l); colnames(Sigma_Y) <- latent; rownames(Sigma_Y) <- latent
Sigma_X <- matrix(Sigma[overt, overt], o, o); colnames(Sigma_X) <- overt; rownames(Sigma_X) <- overt
Sigma_Y_X <- matrix(Sigma[latent, overt], l, o); colnames(Sigma_Y_X) <- overt; rownames(Sigma_Y_X) <- latent
Sigma_X_Y <- matrix(Sigma[overt, latent], o, l); colnames(Sigma_X_Y) <- latent; rownames(Sigma_X_Y) <- overt

Sigma_Y_c <- matrix(Sigma_Y[latent_c, latent_c], l_c, l_c); colnames(Sigma_Y_c) <- latent_c; rownames(Sigma_Y_c) <- latent_c
Sigma_Y_u <- matrix(Sigma_Y[latent_u, latent_u], l_u, l_u); colnames(Sigma_Y_u) <- latent_u; rownames(Sigma_Y_u) <- latent_u
Sigma_Y_c_u <- matrix(Sigma_Y[latent_c, latent_u], l_c, l_u); colnames(Sigma_Y_c_u) <- latent_u; rownames(Sigma_Y_c_u) <- latent_c
Sigma_Y_u_c <- matrix(Sigma_Y[latent_u, latent_c], l_u, l_c); colnames(Sigma_Y_u_c) <- latent_c; rownames(Sigma_Y_u_c) <- latent_u

Psi <- 1/n * t(bar_Xi) %*% bar_Xi
Psi_Y <- matrix(Psi[latent, latent], l, l); colnames(Psi_Y) <- latent; rownames(Psi_Y) <- latent
Psi_X <- matrix(Psi[overt, overt], o, o); colnames(Psi_X) <- overt; rownames(Psi_X) <- overt
Psi_Y_X <- matrix(Psi[latent, overt], l, o); colnames(Psi_Y_X) <- overt; rownames(Psi_Y_X) <- latent
Psi_X_Y <- matrix(Psi[overt, latent], o, l); colnames(Psi_X_Y) <- latent; rownames(Psi_X_Y) <- overt

Psi_Y_c <- matrix(Psi_Y[latent_c, latent_c], l_c, l_c); colnames(Psi_Y_c) <- latent_c; rownames(Psi_Y_c) <- latent_c
Psi_Y_u <- matrix(Psi_Y[latent_u, latent_u], l_u, l_u); colnames(Psi_Y_u) <- latent_u; rownames(Psi_Y_u) <- latent_u
Psi_Y_c_u <- matrix(Psi_Y[latent_c, latent_u], l_c, l_u); colnames(Psi_Y_c_u) <- latent_u; rownames(Psi_Y_c_u) <- latent_c
Psi_Y_u_c <- matrix(Psi_Y[latent_u, latent_c], l_u, l_c); colnames(Psi_Y_u_c) <- latent_c; rownames(Psi_Y_u_c) <- latent_u


# Validate eq. 51
validate(df_Y, cbind(df_Y_c, df_Y_u)) # TRUE
validate(Lambda, cbind(Lambda_c, Lambda_u)) # TRUE
validate(df_X, df_Y_c %*% t(Lambda_c) + df_Y_u %*% t(Lambda_u)) # TRUE


# Validate eqs. 52-53
validate(Beta_star, diag(1, l, l)) # TRUE
validate(Omega_star, diag(1, o, o)) # TRUE
validate(G_star, diag(1, l, l)) # TRUE
validate(H_star, diag(1, o, o)) # TRUE
validate(H, matrix(0, l, o)) # TRUE
validate(G, rbind(cbind(Lambda_c, Lambda_u))) # TRUE
validate(Sigma_X, G %*% Psi_Y %*% t(G) + H_star %*% Psi_X %*% t(H_star) + G %*% Psi_Y_X %*% t(H_star) + H_star %*% Psi_X_Y %*% t(G)) # TRUE
validate(Sigma_X, G %*% Psi_Y %*% t(G) + G %*% Psi_Y_X %*% t(H_star) + H_star %*% Psi_X_Y %*% t(G)) # TRUE
validate(Sigma_X, G %*% Psi_Y %*% t(G) + G %*% Psi_Y_X %*% t(H_star)) # TRUE
validate(Sigma_X, G %*% Psi_Y %*% t(G)) # TRUE
validate(Sigma_X, cbind(Lambda_c, Lambda_u) %*% rbind(cbind(Psi_Y_c, Psi_Y_c_u), cbind(Psi_Y_u_c, Psi_Y_u)) %*% t(cbind(Lambda_c, Lambda_u))) # TRUE
validate(Sigma_X, Lambda_c %*% Psi_Y_c %*% t(Lambda_c) + Lambda_u %*% Psi_Y_u %*% t(Lambda_u)) # TRUE
validate(Sigma_X, Lambda_c %*% Psi_Y_c %*% t(Lambda_c) + Lambda_u %*% t(Lambda_u)) # TRUE


# End of derivation/proof #





# Classical Fundamental Equation of Component Analysis --------------------

n <- 1000
o <- 12
l <- o # NOTE: Assume as many components as overt variables
k <- l + o

overt <- paste0("o", 1:o)
latent <- paste0("l", 1:l)
variables <- c(latent, overt)

alpha <- matrix(rnorm(k), 1, k); colnames(alpha) <- variables
alpha_Y <- matrix(0, 1, l); colnames(alpha_Y) <- latent # NOTE: Changed to a one by l null matrix
alpha_X <- matrix(alpha[,overt], 1, o); colnames(alpha_X) <- overt
alpha <- cbind(alpha_Y, alpha_X)

while(TRUE){
  # Simulate Beta
  while(TRUE){
    mu_Beta <- matrix(rnorm(l), 1, l) # NOTE: Unused
    Sigma_Beta <- trialr::rlkjcorr(1, l, 1) # NOTE: Unused
    Beta <- matrix(0, l, l) # NOTE: Changed to an l by l null matrix
    diag(Beta) <- 0
    if(abs(det(diag(1, l, l) - Beta)) > 0){
      colnames(Beta) <- latent
      rownames(Beta) <- latent
      break
    } 
  }
  # Simulate Delta
  mu_Delta <- matrix(rnorm(o), 1, o)
  Sigma_Delta <- trialr::rlkjcorr(1, o, 1)
  Delta <- matrix(brms::rmulti_normal(l, mu_Delta, Sigma_Delta), l, o); colnames(Delta) <- overt; rownames(Delta) <- latent
  
  # Simulate Lambda
  mu_Lambda <- matrix(rnorm(l), 1, l) # NOTE: Unused
  Sigma_Lambda <- trialr::rlkjcorr(1, l, 1) # NOTE: Unused
  Lambda <- matrix(0, o, l); colnames(Lambda) <- latent; rownames(Lambda) <- overt  # NOTE: Changed to an o by o null matrix
  
  # Simulate Omega
  while(TRUE){
    mu_Omega <- matrix(rnorm(o), 1, o) # NOTE: Unused
    Sigma_Omega <- trialr::rlkjcorr(1, o, 1) # NOTE: Unused
    Omega <- matrix(0, o, o) # NOTE: Changed to an o by o null matrix
    diag(Omega) <- 0
    if(abs(det(diag(1, o, o) - Omega)) > 0){
      colnames(Omega) <- overt
      rownames(Omega) <- overt
      break
    } 
  }
  A <- rbind(cbind(Beta, Delta), cbind(Lambda, Omega))
  if(abs(det(diag(1, k) - A)) > 0){
    A_star <- solve(diag(1, k) - A)
    break
  }
}

Beta_star <- solve(diag(1, l, l) - Beta)
Omega_star <- solve(diag(1, o, o) - Omega)
H_star <- solve(diag(1, o, o) - Omega - Lambda %*% Beta_star %*% Delta)
G_star <- solve(diag(1, l, l) - Beta - Delta %*% Omega_star %*% Lambda)
H <- Beta_star %*% Delta %*% H_star
G <- Omega_star %*% Lambda %*% G_star

mu_Xi <- matrix(rnorm(k), 1, k)
Psi <- trialr::rlkjcorr(1, k, 1); colnames(Psi) <- variables; rownames(Psi) <- variables
Xi <- matrix(brms::rmulti_normal(n, mu_Xi, Psi), n, k); colnames(Xi) <- variables
Xi_Y <- matrix(0, n, l); colnames(Xi_Y) <- latent # NOTE: Changed to a null matrix
Xi_X <- matrix(Xi[,overt], n, o); colnames(Xi_X) <- overt
Xi <- cbind(Xi_Y, Xi_X) # Update Xi

mu_Xi <- matrix(colMeans(Xi), 1, k); colnames(mu_Xi) <- variables
mu_Xi_Y <- matrix(mu_Xi[,latent], 1, l); colnames(mu_Xi_Y) <- latent
mu_Xi_X <- matrix(mu_Xi[,overt], 1, o); colnames(mu_Xi_X) <- overt

bar_Xi <- Xi - matrix(1, n, 1) %*% mu_Xi
bar_Xi_Y <- matrix(bar_Xi[,latent], n, l); colnames(bar_Xi_Y) <- latent
bar_Xi_X <- matrix(bar_Xi[,overt], n, o); colnames(bar_Xi_X) <- overt

df <- (matrix(1, n, 1) %*% alpha + Xi) %*% t(A_star)
df_Y <- matrix(df[,latent], n, l); colnames(df_Y) <- latent
df_X <- matrix(df[,overt], n, o); colnames(df_X) <- overt

mu <- matrix(colMeans(df), 1, k); colnames(mu) <- variables
mu_Y <- matrix(mu[,latent], 1, l); colnames(mu_Y) <- latent
mu_X <- matrix(mu[,overt], 1, o); colnames(mu_X) <- overt

bar_df <- df - matrix(1, n, 1) %*% mu
bar_Y <- matrix(bar_df[,latent], n, l); colnames(bar_Y) <- latent
bar_X <- matrix(bar_df[,overt], n, o); colnames(bar_X) <- overt

E <- df - matrix(1, n, 1) %*% mu
E_Y <- matrix(E[,latent], n, l); colnames(E_Y) <- latent
E_X <- matrix(E[,overt], n, o); colnames(E_X) <- overt

mu_E <- matrix(colMeans(E), 1, k); colnames(mu_E) <- variables
mu_E_Y <- matrix(mu_E[,latent], 1, l); colnames(mu_E_Y) <- latent
mu_E_X <- matrix(mu_E[,overt], 1, o); colnames(mu_E_X) <- overt

bar_E <- E - matrix(1, n, 1) %*% mu_E
bar_E_Y <- matrix(bar_E[,latent], n, l); colnames(bar_E_Y) <- latent
bar_E_X <- matrix(bar_E[,overt], n, o); colnames(bar_E_X) <- overt

Sigma <- 1/n * t(bar_E) %*% bar_E
Sigma_Y <- matrix(Sigma[latent, latent], l, l); colnames(Sigma_Y) <- latent; rownames(Sigma_Y) <- latent
Sigma_X <- matrix(Sigma[overt, overt], o, o); colnames(Sigma_X) <- overt; rownames(Sigma_X) <- overt
Sigma_Y_X <- matrix(Sigma[latent, overt], l, o); colnames(Sigma_Y_X) <- overt; rownames(Sigma_Y_X) <- latent
Sigma_X_Y <- matrix(Sigma[overt, latent], o, l); colnames(Sigma_X_Y) <- latent; rownames(Sigma_X_Y) <- overt

Psi <- 1/n * t(bar_Xi) %*% bar_Xi
Psi_Y <- matrix(Psi[latent, latent], l, l); colnames(Psi_Y) <- latent; rownames(Psi_Y) <- latent
Psi_X <- matrix(Psi[overt, overt], o, o); colnames(Psi_X) <- overt; rownames(Psi_X) <- overt
Psi_Y_X <- matrix(Psi[latent, overt], l, o); colnames(Psi_Y_X) <- overt; rownames(Psi_Y_X) <- latent
Psi_X_Y <- matrix(Psi[overt, latent], o, l); colnames(Psi_X_Y) <- latent; rownames(Psi_X_Y) <- overt

# Validate eq. 54
validate(df_Y, df_X %*% t(Delta)) # TRUE


# End of derivation/proof #





# Classical Standard Linear Model -----------------------------------------

n <- 1000
o_eta <- 6
o_chi <- 6
o <- o_eta + o_chi
l <- 3
k <- l + o

overt_eta <- paste0("o_eta", 1:o_eta)
overt_chi <- paste0("o_chi", 1:o_chi)
overt <- c(overt_eta, overt_chi)
latent <- paste0("l", 1:l)
variables <- c(latent, overt)

alpha <- matrix(0, 1, k); colnames(alpha) <- variables  # Changed to a 1 by k null matrix
alpha_Y <- matrix(alpha[,latent], 1, l); colnames(alpha_Y) <- latent
alpha_X <- matrix(alpha[,overt], 1, o); colnames(alpha_X) <- overt
alpha_X_eta <- matrix(alpha_X[,overt_eta], 1, o_eta); colnames(alpha_X_eta) <- overt_eta
alpha_X_chi <- matrix(alpha_X[,overt_chi], 1, o_chi); colnames(alpha_X_chi) <- overt_chi

while(TRUE){
  # Simulate Beta
  while(TRUE){
    mu_Beta <- matrix(rnorm(l), 1, l) # NOTE: Unused here
    Sigma_Beta <- trialr::rlkjcorr(1, l, 1) # NOTE: Unused here
    Beta <- matrix(0, l, l) # NOTE: Changed to an l by l null matrix
    diag(Beta) <- 0
    if(abs(det(diag(1, l, l) - Beta)) > 0){
      colnames(Beta) <- latent
      rownames(Beta) <- latent
      break
    } 
  }
  # Simulate Delta
  mu_Delta <- matrix(rnorm(o), 1, o) # NOTE: Unused here
  Sigma_Delta <- trialr::rlkjcorr(1, o, 1) # NOTE: Unused here
  Delta <- matrix(0, l, o); colnames(Delta) <- overt; rownames(Delta) <- latent # NOTE: Changed to an l by o null matrix
  
  # Simulate Lambda
  mu_Lambda <- matrix(rnorm(l), 1, l) # NOTE: Unused here
  Sigma_Lambda <- trialr::rlkjcorr(1, l, 1) # NOTE: Unused here
  Lambda <- matrix(0, o, l); colnames(Lambda) <- latent; rownames(Lambda) <- overt # NOTE: Changed to an o by l null matrix
  
  # Simulate Omega
  while(TRUE){
    mu_Omega <- matrix(rnorm(o), 1, o)
    Sigma_Omega <- trialr::rlkjcorr(1, o, 1)
    Omega <- matrix(brms::rmulti_normal(o, mu_Omega, Sigma_Omega), o, o)
    Omega[lower.tri(Omega)] <- 0 # NOTE: New, Omega is made strictly recursive. Conceptual, not strictly necessary
    diag(Omega) <- 0
    if(abs(det(diag(1, o, o) - Omega)) > 0){
      colnames(Omega) <- overt
      rownames(Omega) <- overt
      Omega_eta <- Omega[overt_eta, overt_eta] # NOTE: New
      Omega_chi <- Omega[overt_chi, overt_chi] # NOTE: New
      Omega_eta_chi <- Omega[overt_eta, overt_chi] # NOTE: New
      Omega_chi_eta <- Omega[overt_chi, overt_eta] # NOTE: New
      validate(Omega, rbind(cbind(Omega_eta, Omega_eta_chi), cbind(Omega_chi_eta, Omega_chi))) # NOTE: New
      Omega_chi_eta <- matrix(0, o_chi, o_eta) # NOTE: Changed so that Omega_chi_eta is a null matrix. Conceptual, not strictly necessary
      Omega_eta <- matrix(0, o_eta, o_eta) # NOTE: Changed so that Omega_eta is a null matrix
      Omega_chi <- matrix(0, o_chi, o_chi) # NOTE: Changed so that Omega_chi is a null matrix. Conceptual, not strictly necessary
      Omega <- rbind(cbind(Omega_eta, Omega_eta_chi), cbind(Omega_chi_eta, Omega_chi)) # NOTE: Omega is updated
      break
    } 
  }
  A <- rbind(cbind(Beta, Delta), cbind(Lambda, Omega))
  if(abs(det(diag(1, k) - A)) > 0){
    A_star <- solve(diag(1, k) - A)
    break
  }
}

Beta_star <- solve(diag(1, l, l) - Beta)
Omega_star <- solve(diag(1, o, o) - Omega)
H_star <- solve(diag(1, o, o) - Omega - Lambda %*% Beta_star %*% Delta)
G_star <- solve(diag(1, l, l) - Beta - Delta %*% Omega_star %*% Lambda)
H <- Beta_star %*% Delta %*% H_star
G <- Omega_star %*% Lambda %*% G_star

mu_Xi <- matrix(rnorm(k), 1, k)
Psi <- trialr::rlkjcorr(1, k, 1); colnames(Psi) <- variables; rownames(Psi) <- variables
Xi <- matrix(brms::rmulti_normal(n, mu_Xi, Psi), n, k); colnames(Xi) <- variables
Xi[,latent] <- 0 # New, Xi_Y is assumed a null matrix
Xi_Y <- matrix(Xi[,latent], n, l); colnames(Xi_Y) <- latent
Xi_X <- matrix(Xi[,overt], n, o); colnames(Xi_X) <- overt
Xi_X_eta <- matrix(Xi_X[,overt_eta], n, o_eta); colnames(Xi_X_eta) <- overt_eta
Xi_X_chi <- matrix(Xi_X[,overt_chi], n, o_chi); colnames(Xi_X_chi) <- overt_chi

mu_Xi <- matrix(colMeans(Xi), 1, k); colnames(mu_Xi) <- variables
mu_Xi_Y <- matrix(mu_Xi[,latent], 1, l); colnames(mu_Xi_Y) <- latent
mu_Xi_X <- matrix(mu_Xi[,overt], 1, o); colnames(mu_Xi_X) <- overt
mu_Xi_X_eta <- matrix(mu_Xi_X[,overt_eta], 1, o_eta); colnames(mu_Xi_X_eta) <- overt_eta
mu_Xi_X_chi <- matrix(mu_Xi_X[,overt_chi], 1, o_chi); colnames(mu_Xi_X_chi) <- overt_chi

bar_Xi <- Xi - matrix(1, n, 1) %*% mu_Xi
bar_Xi_Y <- matrix(bar_Xi[,latent], n, l); colnames(bar_Xi_Y) <- latent
bar_Xi_X <- matrix(bar_Xi[,overt], n, o); colnames(bar_Xi_X) <- overt
bar_Xi_X_eta <- matrix(bar_Xi_X[,overt_eta], n, o_eta); colnames(bar_Xi_X_eta) <- overt_eta
bar_Xi_X_chi <- matrix(bar_Xi_X[,overt_chi], n, o_chi); colnames(bar_Xi_X_chi) <- overt_chi

df <- (matrix(1, n, 1) %*% alpha + Xi) %*% t(A_star)
df_Y <- matrix(df[,latent], n, l); colnames(df_Y) <- latent
df_X <- matrix(df[,overt], n, o); colnames(df_X) <- overt
df_X_eta <- matrix(df_X[,overt_eta], n, o_eta); colnames(df_X_eta) <- overt_eta
df_X_chi <- matrix(df_X[,overt_chi], n, o_chi); colnames(df_X_chi) <- overt_chi

mu <- matrix(colMeans(df), 1, k); colnames(mu) <- variables
mu_Y <- matrix(mu[,latent], 1, l); colnames(mu_Y) <- latent
mu_X <- matrix(mu[,overt], 1, o); colnames(mu_X) <- overt
mu_X_eta <- matrix(mu_X[,overt_eta], 1, o_eta); colnames(mu_X_eta) <- overt_eta
mu_X_chi <- matrix(mu_X[,overt_chi], 1, o_chi); colnames(mu_X_chi) <- overt_chi

bar_df <- df - matrix(1, n, 1) %*% mu
bar_Y <- matrix(bar_df[,latent], n, l); colnames(bar_Y) <- latent
bar_X <- matrix(bar_df[,overt], n, o); colnames(bar_X) <- overt
bar_X_eta <- matrix(bar_X[,overt_eta], n, o_eta); colnames(bar_X_eta) <- overt_eta
bar_X_chi <- matrix(bar_X[,overt_chi], n, o_chi); colnames(bar_X_chi) <- overt_chi

E <- df - matrix(1, n, 1) %*% mu
E_Y <- matrix(E[,latent], n, l); colnames(E_Y) <- latent
E_X <- matrix(E[,overt], n, o); colnames(E_X) <- overt
E_X_eta <- matrix(E_X[,overt_eta], n, o_eta); colnames(E_X_eta) <- overt_eta
E_X_chi <- matrix(E_X[,overt_chi], n, o_chi); colnames(E_X_chi) <- overt_chi

mu_E <- matrix(colMeans(E), 1, k); colnames(mu_E) <- variables
mu_E_Y <- matrix(mu_E[,latent], 1, l); colnames(mu_E_Y) <- latent
mu_E_X <- matrix(mu_E[,overt], 1, o); colnames(mu_E_X) <- overt
mu_E_X_eta <- matrix(mu_E_X[,overt_eta], 1, o_eta); colnames(mu_E_X_eta) <- overt_eta
mu_E_X_chi <- matrix(mu_E_X[,overt_chi], 1, o_chi); colnames(mu_E_X_chi) <- overt_chi

bar_E <- E - matrix(1, n, 1) %*% mu_E
bar_E_Y <- matrix(bar_E[,latent], n, l); colnames(bar_E_Y) <- latent
bar_E_X <- matrix(bar_E[,overt], n, o); colnames(bar_E_X) <- overt
bar_E_X_eta <- matrix(bar_E_X[,overt_eta], n, o_eta); colnames(bar_E_X_eta) <- overt_eta
bar_E_X_chi <- matrix(bar_E_X[,overt_chi], n, o_chi); colnames(bar_E_X_chi) <- overt_chi

Sigma <- 1/n * t(bar_E) %*% bar_E
Sigma_Y <- matrix(Sigma[latent, latent], l, l); colnames(Sigma_Y) <- latent; rownames(Sigma_Y) <- latent
Sigma_X <- matrix(Sigma[overt, overt], o, o); colnames(Sigma_X) <- overt; rownames(Sigma_X) <- overt
Sigma_Y_X <- matrix(Sigma[latent, overt], l, o); colnames(Sigma_Y_X) <- overt; rownames(Sigma_Y_X) <- latent
Sigma_X_Y <- matrix(Sigma[overt, latent], o, l); colnames(Sigma_X_Y) <- latent; rownames(Sigma_X_Y) <- overt

Sigma_X_eta <- matrix(Sigma_X[overt_eta, overt_eta], o_eta, o_eta); colnames(Sigma_X_eta) <- overt_eta; rownames(Sigma_X_eta) <- overt_eta
Sigma_X_chi <- matrix(Sigma_X[overt_chi, overt_chi], o_chi, o_chi); colnames(Sigma_X_chi) <- overt_chi; rownames(Sigma_X_chi) <- overt_chi
Sigma_X_eta_chi <- matrix(Sigma_X[overt_eta, overt_chi], o_eta, o_chi); colnames(Sigma_X_eta_chi) <- overt_chi; rownames(Sigma_X_eta_chi) <- overt_eta
Sigma_X_chi_eta <- matrix(Sigma_X[overt_chi, overt_eta], o_chi, o_eta); colnames(Sigma_X_chi_eta) <- overt_eta; rownames(Sigma_X_chi_eta) <- overt_chi

Psi <- 1/n * t(bar_Xi) %*% bar_Xi
Psi_Y <- matrix(Psi[latent, latent], l, l); colnames(Psi_Y) <- latent; rownames(Psi_Y) <- latent
Psi_X <- matrix(Psi[overt, overt], o, o); colnames(Psi_X) <- overt; rownames(Psi_X) <- overt
Psi_Y_X <- matrix(Psi[latent, overt], l, o); colnames(Psi_Y_X) <- overt; rownames(Psi_Y_X) <- latent
Psi_X_Y <- matrix(Psi[overt, latent], o, l); colnames(Psi_X_Y) <- latent; rownames(Psi_X_Y) <- overt

Psi_X_eta <- matrix(Psi_X[overt_eta, overt_eta], o_eta, o_eta); colnames(Psi_X_eta) <- overt_eta; rownames(Psi_X_eta) <- overt_eta
Psi_X_chi <- matrix(Psi_X[overt_chi, overt_chi], o_chi, o_chi); colnames(Psi_X_chi) <- overt_chi; rownames(Psi_X_chi) <- overt_chi
Psi_X_eta_chi <- matrix(Psi_X[overt_eta, overt_chi], o_eta, o_chi); colnames(Psi_X_eta_chi) <- overt_chi; rownames(Psi_X_eta_chi) <- overt_eta
Psi_X_chi_eta <- matrix(Psi_X[overt_chi, overt_eta], o_chi, o_eta); colnames(Psi_X_chi_eta) <- overt_eta; rownames(Psi_X_chi_eta) <- overt_chi

# Validate special assumptions
validate(df_X, cbind(df_X_eta, df_X_chi)) # TRUE
validate(Xi_X, cbind(Xi_X_eta, Xi_X_chi)) # TRUE
validate(Omega, rbind(cbind(Omega_eta, Omega_eta_chi), cbind(Omega_chi_eta, Omega_chi))) # TRUE
validate(Omega_chi_eta, matrix(0, o_chi, o_eta)) # TRUE
validate(Omega_eta, matrix(0, o_eta, o_eta)) # TRUE

# Validate eq. 55
validate(df_X_eta, df_X_chi %*% t(Omega_eta_chi) + Xi_X_eta) # TRUE



# End of derivation/proof #





#################
# END OF SCRIPT #
#################