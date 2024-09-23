#install.packages("telescope")
#installed.packages("extraDistr")
#
rm(list = ls())
#
library("telescope")
library("stats")
library("extraDistr")       # Package for the "dbnbinom" probability mass function

source("./sampleUniNormMixture.R")
source("./priorOnK_spec.R")
source("./prior_alpha_e0.R")
source("./sample_e0_alpha.R")
source("./sampleK.R")
source("./identifyMixture.R")


############################
## (1) Set the Parameters ##
############################
set.seed(123)
TT <- 3             # Number of time series observations
N <- 50            # Number of cross-section observations
MC <- 3               # Number of Monte Carlo iterations
#
sigma_Z <- 1        # standard deviation of the covariates Z
beta <- 0           # value of the common parameter
## Parameters of the mixture:
K <- 3              # Number of the mixture components
w_true <- c(1/3,1/3,1/3) # Vector of weights of each component of the mixture
alpha_true <- c(-2,0,1.5)  # True values of the incidental parameter for the K mixture components
var_u <- c(1,1,1)          # True values of the variance of the model error for the K mixture components
## Parameters of the MCMC sampler:
Mmax <- 100   # the maximum number of iterations
thin <- 1       # the thinning imposed to reduce auto-correlation in the chain by only recording every `thined' observation
burnin <- 100   # the number of burn-in iterations not recorded
M <- Mmax/thin  # the number of recorded observations.

Kmax <- 50      # the maximum number of components possible during sampling
Kinit <- 10     # the initial number of filled components

## Initial values:
# we use initial values corresponding to equal component sizes and half of the value of `C0` for the variances.

eta_0 <- rep(1/Kinit, Kinit)
r <- 1    # dimension
c0 <- 2 + (r-1)/2
g0 <- 0.2 + (r-1) / 2

# Initial value for beta
beta_0 <- beta

###############
## (2) Prior ##
###############
# Static specification for the weights with a fixed prior on $e_0$ where the value is set to 1.
G <- "MixStatic" 

priorOn_n0 <- priorOnE0_spec("e0const", 1)        # Prior on the hyperparameter of the Dirichlet prior

priorOn_K <- priorOnK_spec("BNB_143", 30)          # Prior on K. This gives the function to compute log(PMF)-log(.5).

# Prior on the common coefficient beta: N(beta_prior,Omega_prior)
beta_prior <- 0
Omega_prior <- 1

# Initialization of matrices to store the results
theta_alpha <- array(,dim=c(M,Kmax,MC))
theta_Sigma2 <- array(,dim=c(M,Kmax,MC))
Beta <- matrix(,M,MC)
Eta <- array(,dim=c(M,Kmax,MC))
S_label <- array(,dim=c(M,N,MC))
Nk <- array(,dim=c(M,Kmax,MC))
K_MC <- matrix(,M,MC)
#
Kplus <- matrix(,M,MC) 
Kplus_stat <- matrix(,MC,4)
K_stat <- matrix(,MC,5)
#
# Estimators
theta_alpha_hat <- matrix(,MC,Kmax)
theta_Sigma2_hat <- matrix(,MC,Kmax)
Eta_hat <- matrix(,MC,Kmax)
Beta_hat <- matrix(,MC,3)

#####################
## (3) Monte Carlo ##
#####################
for (mc in 1:MC){
  #print(mc)
  cat("mc = ",mc)
  ## (3.1) Simulate the Z and Y 
  #
  ZZ <- rnorm(N*TT, mean = 0, sd = sigma_Z)
  Z <- matrix(ZZ,TT,N)
  
  rm(ZZ)
  S_unif <- runif(N, min = 0, max = 1)
  S <- matrix(,N,1)               # S is the latent allocation variable 
  Y <- matrix(,TT,N)
  
  for (ii in 1:N){
    # Fill the latent allocation variable S
    if (S_unif[ii] >= 0 & S_unif[ii] < w_true[1]) {
      S[ii,1] <- 1
    } else if (S_unif[ii] >= w_true[1] & S_unif[ii] < (w_true[1] + w_true[2])) {
      S[ii,1] <- 2
    } else {
      S[ii,1] <- 3
    }
    u <- rnorm(TT, mean = 0, sd = sqrt(var_u[S[ii]]))
    Y[,ii] <- alpha_true[S[ii]] + beta*Z[,ii] + u
  }
  
  rm(S,S_unif,u)
  Y_mean <- colMeans(Y)
  
  ## (3.2) Initial values of the other parameters (not specified outside the loop):
  # Component-specific priors on \alpha_k (i.e. N(b0,B0)) and \sigma_k^2 (i.e. IG(c_0,C_0)) following Richardson and Green (1997).
  R <- diff(range(Y))
  C0 <- diag(c(0.02*(R^2)), nrow = r)
  sigma2_0 <- array(0, dim = c(1, 1, Kinit))
  sigma2_0[1, 1, ] <- 0.5 * C0
  G0 <- diag(10/(R^2), nrow = r)
  B0 <- diag((R^2), nrow = r)                       # prior variance of alpha_k
  b0 <- as.matrix((max(Y) + min(Y))/2, ncol = 1)    # prior mean of alpha_k 

  ## (3.3) Initial partition of the data and initial parameter values to start the MCMC. 
  # We use `kmeans()` to determine the initial partition $S_0$ and the initial component-specific means \mu_0.
  
  cl_y <- kmeans(Y_mean, centers = Kinit, nstart = 30)
  S_0 <- cl_y$cluster
  alpha_0 <- t(cl_y$centers)
  
  ## (3.3) MCMC sampling 
  # We draw samples from the posterior using the telescoping sampler of Fruhwirth-Schnatter. 
  estGibbs <- sampleUniNormMixture(
    Y, Z, S_0, alpha_0, sigma2_0, eta_0, beta_0,
    c0, g0, G0, C0, b0, B0, beta_prior, Omega_prior,
    M, burnin, thin, Kmax,
    G, priorOn_K, priorOn_n0)
  
  #The sampling function returns a named list where the sampled parameters and latent variables are contained. 
  theta_alpha[,,mc] <- estGibbs$Mu
  theta_Sigma2[,,mc] <- estGibbs$Sigma2
  Beta[,mc] <- estGibbs$Beta
  Eta[,,mc] <- estGibbs$Eta
  S_label[,,mc] <- estGibbs$S
  Nk[,,mc] <- estGibbs$Nk
  K_MC[,mc] <- estGibbs$K
  Kplus[,mc] <- estGibbs$Kp
  nonnormpost_mode_list <- estGibbs$nonnormpost_mode_list
  acc <- estGibbs$acc
  e0 <- estGibbs$e0
  alpha_Dir <- estGibbs$alpha  
  
  #########################################
  ## Identification of the mixture model ##
  #########################################
  ## Step 1: Estimating $K_+$ and $K$
  ## K_+ ##
  Nk_mc <- Nk[,,mc]
  Kplus_mc <- rowSums(Nk_mc != 0, na.rm = TRUE)  
  p_Kplus <- tabulate(Kplus_mc, nbins = max(Kplus_mc))/M  
  # The distribution of $K_+$ is characterized using the 1st and 3rd quartile as well as the median.
  Kplus_stat[mc,1:3] <- quantile(Kplus_mc, probs = c(0.25, 0.5, 0.75))
  
  # Point estimate for K_+ by taking the mode and determine the number of MCMC draws where exactly K_+ 
  # components were filled.
  Kplus_hat <- which.max(p_Kplus)
  Kplus_stat[mc,4] <- Kplus_hat
  M0 <- sum(Kplus_mc == Kplus_hat)
  
  ## K ##
  # We also determine the posterior distribution of the number of components K directly drawn using the 
  # telescoping sampler.
  K_mc <- K_MC[,mc]
  p_K <- tabulate(K_mc, nbins = max(K_mc, na.rm = TRUE))/M
  round(p_K[1:20], digits = 2)
  
  # The posterior mode can be determined as well as the posterior mean and quantiles of the posterior.
  K_hat <- which.max(tabulate(K_mc, nbins = max(K_mc)))
  K_stat[mc,4] <- K_hat
  K_stat[mc,5] <- mean(K_mc)
  K_stat[mc,1:3] <- quantile(K_mc, probs = c(0.25, 0.5, 0.75))
  
  ## Step 2: Extracting the draws with exactly \hat{K}_+ non-empty components
  # Select those draws where the number of filled components was exactly \hat{K}_+:
  index <- Kplus_mc == Kplus_hat
  Nk_mc[is.na(Nk_mc)] <- 0
  Nk_Kplus <- (Nk_mc[index, ] > 0)
  
  # We extract the cluster means, data cluster sizes and cluster assignments for the draws where exactly 
  # \hat{K}_+ components were filled.
  Mu_inter <- estGibbs$Mu[index, , , drop = FALSE]
  Mu_Kplus <- array(0, dim = c(M0, 1, Kplus_hat)) 
  Mu_Kplus[, 1, ] <- Mu_inter[, 1, ][Nk_Kplus]          # Here, 'Nk_Kplus' is used to subset the extracted 
  # elements from Mu_inter. It selects specific rows  
  # from the first column of Mu_inter based on the 
  # indices in Nk_Kplus.
  
  Sigma2_inter <- estGibbs$Sigma2[index, , , drop = FALSE]
  Sigma2_Kplus <- array(0, dim = c(M0, 1, Kplus_hat)) 
  Sigma2_Kplus[, 1, ] <- Sigma2_inter[, 1, ][Nk_Kplus]
  
  Eta_inter <- estGibbs$Eta[index, ]
  Eta_Kplus <- matrix(Eta_inter[Nk_Kplus], ncol = Kplus_hat)
  Eta_Kplus <- sweep(Eta_Kplus, 1, rowSums(Eta_Kplus), "/")     # it normalizes the weights eta
  
  w <- which(index)
  S_Kplus <- matrix(0, M0, ncol(Y))
  for (i in seq_along(w)) {
    m <- w[i]
    perm_S <- rep(0, Kmax)
    perm_S[Nk_mc[m, ] != 0] <- 1:Kplus_hat         # It assigns the sequence '1:Kplus_hat' to the positions in 
    # perm_S where Nk[m, ] is not zero.
    S_Kplus[i, ] <- perm_S[S_label[m, ,mc]]              # 'perm_S[S[m, ]]' indexes the vector 'perm_S' using the 
    # values in S[m, ]. This permutes or reassigns the values in 
    # 'S[m, ]' based on the mapping in 'perm_S'.
  }
  
  ## Step 3: Clustering and relabeling of the MCMC draws in the point process representation
  # For model identification, we cluster the draws of the means where exactly \hat{K}_+ components were 
  # filled in the point process representation using k-means clustering. 
  Func_init <- nonnormpost_mode_list[[Kplus_hat]]$mu  
  identified_Kplus <- identifyMixture(
    Mu_Kplus, Mu_Kplus, Eta_Kplus, S_Kplus, Sigma2_Kplus, Func_init)
  
  #A named list is returned which contains the proportion of draws where the clustering did not result in a 
  # permutation and hence no relabeling could be performed and the draws had to be omitted.
  identified_Kplus$non_perm_rate
  
  
  ## Step 4: Estimating data cluster specific parameters and determining the final partition 
  # The relabeled draws are also returned which can be used to determine posterior mean values for data 
  # cluster specific parameters.
  Mu_mean <- colMeans(identified_Kplus$Mu)
  theta_alpha_hat[mc,1:length(Mu_mean)] <- colMeans(identified_Kplus$Mu)
  theta_Sigma2_hat[mc,1:length(Mu_mean)] <- colMeans(identified_Kplus$Sigma2)
  Eta_hat[mc,1:length(colMeans(identified_Kplus$Eta))] <- colMeans(identified_Kplus$Eta)
  Beta_hat[mc,1] <- mean(Beta[,mc])
  Beta_hat[mc,2:3] <- quantile(Beta[,mc], probs = c(0.025,0.975), na.rm =TRUE)
  
  rm(R, C0, sigma2_0, G0, B0, b0, cl_y, S_0, alpha_0, estGibbs, nonnormpost_mode_list, acc, e0, alpha_Dir, 
     Kplus_mc, Nk_mc, p_Kplus, Kplus_hat, K_mc, index, Mu_inter, Mu_Kplus, Nk_Kplus, Sigma2_inter, 
     Sigma2_Kplus, Eta_inter, Eta_Kplus, w, S_Kplus, perm_S, Func_init, identified_Kplus, Mu_mean)
}

theta_alpha_hat2<-t(apply(theta_alpha_hat[,1:3],1,sort))
theta_alpha_hat_means <- colMeans(theta_alpha_hat2)
Beta_hat_mean <- colMeans(Beta_hat)
MSE_beta = mean((Beta_hat[,1] - beta)^2)
MSE_alpha = colMeans((theta_alpha_hat2[,1:3] - matrix(rep(sort(alpha_true),MC), ncol=3, byrow=TRUE))^2)
RMSE_alpha <- sqrt(MSE_alpha)

Eta_hat2 <- t(matrix(Eta_hat[order(row(theta_alpha_hat[,1:3]), theta_alpha_hat[,1:3])], byrow = FALSE, ncol = length(w_true)))
Eta_hat_means <- colMeans(Eta_hat2)
MSE_eta = colMeans((Eta_hat2[,1:3] - matrix(rep(w_true,MC), ncol=3, byrow=TRUE))^2)

theta_Sigma2_hat2 <- t(matrix(theta_Sigma2_hat[order(row(theta_alpha_hat[,1:3]), theta_alpha_hat[,1:3])], byrow = FALSE, ncol = length(var_u)))
theta_Sigma2_hat_means <- colMeans(theta_Sigma2_hat2)
Sigma2_sorted <- var_u[order(sort(alpha_true))]
MSE_Sigma2 = colMeans((theta_Sigma2_hat2[,1:3] - matrix(rep(Sigma2_sorted,MC), ncol=3, byrow=TRUE))^2)

#MSE_beta = mean(Beta_hat[,1] - beta)^2
#MSE_alpha = colMeans((theta_alpha_hat[,1:3] - matrix(rep(alpha_true,MC), ncol=3, byrow=TRUE))^2)
#MSE_eta = colMeans((Eta_hat[,1:3] - matrix(rep(w_true,MC), ncol=3, byrow=TRUE))^2)

save(theta_alpha_hat_means,Beta_hat_mean,Eta_hat_means,Eta_hat_means,alpha_true,w_true,var_u,TT,N,MC, 
     file = "./N50_T3.RData")

load("./N50_T3.RData")
