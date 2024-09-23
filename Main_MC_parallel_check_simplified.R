rm(list = ls())
#
library("telescope")
library("stats")
library("extraDistr")       # Package for the "dbnbinom" probability mass function
library(parallel)

source("./sampleUniNormMixture.R")
source("./priorOnK_spec.R")
source("./prior_alpha_e0.R")
source("./sample_e0_alpha.R")
source("./sampleK.R")
source("./identifyMixture.R")

# Determine the number of cores availale in my machine
#num_cores <- detectCores() - 1   # Use all but one core

num_cores <- 50

cat(num_cores)

############################
## (1) Set the Parameters ##
############################
set.seed(123)
TT <- 100  #5           # Number of time series observations
N <- 100   #5         # Number of cross-section observations
MC <- 50   #10             # Number of Monte Carlo iterations
#
sigma_Z <- 1        # standard deviation of the covariates Z
beta <- 0           # value of the common parameter
## Parameters of the mixture:
K <- 3              # Number of the mixture components
w_true <- c(1/3,1/3,1/3) # Vector of weights of each component of the mixture
alpha_true <- c(-5,0,5)  # True values of the incidental parameter for the K mixture components
var_u <- c(1,1,1)          # True values of the variance of the model error for the K mixture components
## Parameters of the MCMC sampler:
Mmax <- 10000 #100  # the maximum number of iterations
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


##########################################################################################################
# Initialize a cluster
cl <- makeCluster(num_cores)

# Load the necessary packages on all workers
clusterEvalQ(cl, library(extraDistr))

# Export necessary objects and functions to the cluster
clusterExport(cl, c("sampleUniNormMixture", "priorOnK_spec", "identifyMixture", "priorOnE0_spec",
                    "priorOnAlpha_spec", "sampleE0", "sampleK_e0",
                    "eta_0", "beta_0", "c0", "g0", "M", "burnin", "thin", "Kmax", "G", 
                    "priorOn_K", "priorOn_n0", "N", "TT", "w_true", "alpha_true", "var_u", "MC",
                    "sigma_Z", "beta", "K", "Mmax", 
                    "beta_prior", "Omega_prior", "Kinit","r"))

#####################
## (3) Monte Carlo ##
#####################
# Simplified version to identify the problem
results <- parLapply(cl, 1:MC, function(mc) {
  tryCatch({
    cat("mc = ", mc, "\n")
    
    ## (3.1) Simulate the Z and Y 
    ZZ <- rnorm(N * TT, mean = 0, sd = sigma_Z)
    Z <- matrix(ZZ, TT, N)
    S_unif <- runif(N, min = 0, max = 1)
    S <- matrix(, N, 1)
    Y <- matrix(, TT, N)
    
    for (ii in 1:N) {
      if (S_unif[ii] >= 0 & S_unif[ii] < w_true[1]) {
        S[ii, 1] <- 1
      } else if (S_unif[ii] >= w_true[1] & S_unif[ii] < (w_true[1] + w_true[2])) {
        S[ii, 1] <- 2
      } else {
        S[ii, 1] <- 3
      }
      u <- rnorm(TT, mean = 0, sd = sqrt(var_u[S[ii]]))
      Y[, ii] <- alpha_true[S[ii]] + beta * Z[, ii] + u
    }
    Y_mean <- colMeans(Y)
    
    ## (3.2) Initial values of the other parameters 
    R <- diff(range(Y))
    C0 <- diag(c(0.02 * (R^2)), nrow = r)
    sigma2_0 <- array(0, dim = c(1, 1, Kinit))
    sigma2_0[1, 1, ] <- 0.5 * C0
    G0 <- diag(10 / (R^2), nrow = r)
    B0 <- diag((R^2), nrow = r)
    b0 <- as.matrix((max(Y) + min(Y)) / 2, ncol = 1)
    
    cl_y <- kmeans(Y_mean, centers = Kinit, nstart = 30)
    S_0 <- cl_y$cluster
    alpha_0 <- t(cl_y$centers)
    
    ## (3.3) MCMC sampling 
    estGibbs <- sampleUniNormMixture(Y, Z, S_0, alpha_0, sigma2_0, eta_0, beta_0, c0, g0, G0, C0, b0, B0, 
                                     beta_prior, Omega_prior, M, burnin, thin, Kmax, G, priorOn_K, priorOn_n0)
    

    list(
      theta_alpha = estGibbs$Mu,
      theta_Sigma2 = estGibbs$Sigma2,
      Beta = estGibbs$Beta,
      Eta = estGibbs$Eta,
      S_label = estGibbs$S,
      Nk =  estGibbs$Nk,
      K_MC = estGibbs$K,
      Kplus = estGibbs$Kp,
      nonnormpost_mode_list = estGibbs$nonnormpost_mode_list
#      theta_alpha_hat = colMeans(identified_Kplus$Mu),
#      theta_Sigma2_hat = colMeans(identified_Kplus$Sigma2),
#      Eta_hat = colMeans(identified_Kplus$Eta),
#      Beta_hat = c(mean(estGibbs$Beta), quantile(estGibbs$Beta, probs = c(0.025, 0.975), na.rm = TRUE)),
#      Kplus_stat = c(quantile(Kplus_mc, probs = c(0.25, 0.5, 0.75)), Kplus_hat),
#      K_stat = c(quantile(K_mc, probs = c(0.25, 0.5, 0.75)), K_hat, mean(K_mc))
    )
    
#    print(c(theta_alpha[,,mc],theta_Sigma2[,,mc], Beta[,mc], Eta[,,mc], S_label[,,mc]))
#  }, error = function(e) {
#    list(error = e$message, trace = traceback())
  })
})

# Stop the cluster
stopCluster(cl)


# Initialization of matrices to store the results
theta_alpha <- array(,dim=c(M,Kmax,MC))
theta_Sigma2 <- array(,dim=c(M,Kmax,MC))
Beta <- matrix(,M,MC)
Eta <- array(,dim=c(M,Kmax,MC))
S_label <- array(,dim=c(M,N,MC))
Nk <- array(,dim=c(M,Kmax,MC))
K_MC <- matrix(,M,MC)
Kplus <- matrix(,M,MC) 
nonnormpost_mode_list <- array(list(), MC)

# Combine the results from the list
for (mc in 1:MC) {
  theta_alpha[,,mc] <- results[[mc]]$theta_alpha
  theta_Sigma2[,,mc] <- results[[mc]]$theta_Sigma2
  Beta[,mc] <- results[[mc]]$Beta
  Eta[,,mc] <- results[[mc]]$Eta
  S_label[,,mc] <- results[[mc]]$S_label
  Nk[,,mc] <- results[[mc]]$Nk
  K_MC[,mc] <- results[[mc]]$K_MC
  Kplus[,mc] <- results[[mc]]$Kplus
  nonnormpost_mode_list[[mc]] <- as.list(results[[mc]]$nonnormpost_mode_list)
}

# Post-processing
# Initialization of matrices
Kplus_stat <- matrix(,MC,4)
K_stat <- matrix(,MC,4)
#
# Estimators
theta_alpha_hat <- matrix(,MC,Kmax)
theta_Sigma2_hat <- matrix(,MC,Kmax)
Eta_hat <- matrix(,MC,Kmax)
Beta_hat <- matrix(,MC,3)

for (mc in 1:MC){
  Nk_mc <- Nk[,,mc]
  Kplus_mc <- rowSums(Nk_mc != 0, na.rm = TRUE)
  p_Kplus <- tabulate(Kplus_mc, nbins = max(Kplus_mc)) / M
  Kplus_hat <- which.max(p_Kplus)
  Kplus_stat[mc,1:4] <- c(Kplus_hat,quantile(Kplus_mc, probs = c(0.25, 0.5, 0.75)))
  M0 <- sum(Kplus_mc == Kplus_hat)
  
  K_mc <- K_MC[,mc]
  p_K <- tabulate(K_mc, nbins = max(K_mc, na.rm = TRUE)) / M
  K_hat <- which.max(tabulate(K_mc, nbins = max(K_mc)))
  K_stat[mc,1:4] <- c(K_hat, quantile(K_mc, probs = c(0.25, 0.5, 0.75)))
  
  index <- Kplus_mc == Kplus_hat
  Nk_mc[is.na(Nk_mc)] <- 0
  Nk_Kplus <- (Nk_mc[index, ] > 0)
  
  alpha_inter <- results[[mc]]$theta_alpha[index, , , drop = FALSE]
  #alpha_inter <- theta_alpha[index, ,mc , drop = FALSE]
  alpha_Kplus <- array(0, dim = c(M0, 1, Kplus_hat)) 
  alpha_Kplus[, 1, ] <- alpha_inter[, 1, ][Nk_Kplus]
  
  Sigma2_inter <- theta_Sigma2[index, ,mc, drop = FALSE]
  Sigma2_Kplus <- array(0, dim = c(M0, 1, Kplus_hat))
  Sigma2_Kplus[, 1, ] <- Sigma2_inter[, 1, ][Nk_Kplus]
  
  Eta_inter <- Eta[index,,mc]
  Eta_Kplus <- matrix(Eta_inter[Nk_Kplus], ncol = Kplus_hat)
  Eta_Kplus <- sweep(Eta_Kplus, 1, rowSums(Eta_Kplus), "/")
  
  w <- which(index)
  S_Kplus <- matrix(0, M0, N)
  for (i in seq_along(w)) {
    m <- w[i]
    perm_S <- rep(0, Kmax)
    perm_S[Nk_mc[m, ] != 0] <- 1:Kplus_hat
    S_Kplus[i, ] <- perm_S[S_label[m, ,mc]]
  }
  
  Func_init <- nonnormpost_mode_list[[mc]][[Kplus_hat]]$mu  
  identified_Kplus <- identifyMixture(alpha_Kplus, alpha_Kplus, Eta_Kplus, S_Kplus, Sigma2_Kplus, Func_init)
  
  alpha_mean <- colMeans(identified_Kplus$Mu)
  theta_alpha_hat[mc,1:length(alpha_mean)] <- colMeans(identified_Kplus$Mu)
  theta_Sigma2_hat[mc,1:length(alpha_mean)] <- colMeans(identified_Kplus$Sigma2)
  Eta_hat[mc,1:length(colMeans(identified_Kplus$Eta))] <- colMeans(identified_Kplus$Eta)
  Beta_hat[mc,1] <- mean(Beta[,mc])
  Beta_hat[mc,2:3] <- quantile(Beta[,mc], probs = c(0.025,0.975), na.rm =TRUE)
  
  #Mu_mean <- colMeans(identified_Kplus$Mu)
  
  rm(Nk_mc, Kplus_mc, p_Kplus, Kplus_hat, M0, K_mc, p_K, K_hat, index, Nk_Kplus, alpha_inter,
     alpha_Kplus, Sigma2_inter, Sigma2_Kplus, Eta_inter, Eta_Kplus, w, S_Kplus, perm_S, Func_init)
}
theta_alpha_hat2<-t(apply(theta_alpha_hat,1,sort))
theta_alpha_hat2<-t(apply(theta_alpha_hat[,1:3],1,sort))
theta_alpha_hat_means <- colMeans(theta_alpha_hat2)
MSE_beta = mean((Beta_hat[,1] - beta)^2)
MSE_alpha = colMeans((theta_alpha_hat2[,1:3] - matrix(rep(alpha_true,MC), ncol=3, byrow=TRUE))^2)
RMSE_alpha <- sqrt(MSE_alpha)

# order the rows of eta according to theta_alpha_hat2
Eta_hat2 <- matrix(Eta_hat[order(row(theta_alpha_hat), theta_alpha_hat)], byrow = TRUE, ncol = ncol(Eta_hat))
Eta_hat_means <- colMeans(Eta_hat2)
MSE_eta = colMeans((Eta_hat2[,1:3] - matrix(rep(w_true,MC), ncol=3, byrow=TRUE))^2)


