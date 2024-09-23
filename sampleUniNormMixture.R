#' Telescoping sampling of a Bayesian finite univariate Gaussian mixture where a prior 
#' on the number of components K is specified. 
#'
#' @description
#'   * The MCMC scheme is implemented as suggested in Frühwirth-Schnatter et al (2021).
#'   * The priors on the model parameters are specified as in Frühwirth-Schnatter et al (2021),
#'     see the vignette for details and notation.
#'   * The parametrizations of the gamma and inverse gamma distribution are used as in 
#'     Frühwirth-Schnatter et al (2021), see also the vignette.
#' 
#' @param y A numeric matrix; containing the data.
#' @param S A numeric matrix; containing the initial cluster
#'     assignments.
#' @param mu A numeric matrix; containing the initial cluster-specific
#'     mean values.
#' @param sigma2 A numeric matrix; containing the initial cluster-specific
#'     variance values.
#' @param eta A numeric vector; containing the initial cluster sizes.
#' @param c0 A numeric vector; hyperparameter of the prior on \eqn{\sigma^2_k}.
#' @param g0 A numeric vector; hyperparameter of the prior on \eqn{\sigma^2_k}
#' @param G0 A numeric vector; hyperparameter of the prior on \eqn{\sigma^2_k}
#' @param C0_0 A numeric vector; initial value of hyperparameter \eqn{C_0}.
#' @param b0 A numeric vector; hyperparameter of the prior on \eqn{\mu_k}.
#' @param B0 A numeric vector; hyperparameter of the prior on \eqn{\mu_k}.
#' @param M A numeric scalar; specifying the number of recorded
#'     iterations.
#' @param burnin A numeric scalar; specifying the number of burn-in
#'     iterations.
#' @param thin A numeric scalar; specifying the thinning used for the
#'     iterations.
#' @param Kmax A numeric scalar; the maximum number of components. 
#' @param G A character string; either `"MixDynamic"` or `"MixStatic"`.
#' @param priorOnK A named list; providing the prior on the number of
#'     components K, see [priorOnK_spec()].
#' @param priorOnWeights A named list; providing the prior on the mixture weights.
#' @param verbose A logical; indicating if some intermediate clustering
#'     results should be printed.
#' @return A named list containing:
#'   * `"Mu"`: sampled component means.
#'   * `"Eta"`: sampled weights.
#'   * `"S"`: sampled assignments.
#'   * `"Nk"`: number of observations assigned to the different components, for each iteration.
#'   * `"K"`: sampled number of components.
#'   * `"Kplus"`: number of filled, i.e., non-empty components, for each iteration.
#'   * `"e0"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{e_0} is random).
#'   * `"alpha"`: sampled Dirichlet parameter of the prior on the weights (if \eqn{\alpha} is random).
#'   * `"acc"`: logical vector indicating acceptance in the Metropolis-Hastings step when sampling either \eqn{e_0} or \eqn{\alpha}.
#'
#' @examples
#' if (requireNamespace("mclust", quietly = TRUE)) {
#'     data("acidity", package = "mclust")
#'     y <- acidity
#'     
#'     N <- length(y)
#'     r <- 1
#'     
#'     Mmax <- 200
#'     thin <- 1
#'     burnin <- 100
#'     M <- Mmax/thin
#'     Kmax <- 50  
#'     Kinit <- 10
#'     
#'     G <- "MixStatic" 
#'     priorOnE0 <- priorOnE0_spec("e0const", 0.01)
#'     priorOnK <- priorOnK_spec("Pois_1", 50)
#'     
#'     R <- diff(range(y))
#'     c0 <- 2 + (r-1)/2
#'     C0 <- diag(c(0.02*(R^2)), nrow = r)
#'     g0 <- 0.2 + (r-1) / 2
#'     G0 <- diag(10/(R^2), nrow = r)
#'     B0 <- diag((R^2), nrow = r)
#'     b0 <- as.matrix((max(y) + min(y))/2, ncol = 1)  
#'     
#'     cl_y <- kmeans(y, centers = Kinit, nstart = 100)
#'     S_0 <- cl_y$cluster
#'     mu_0 <- t(cl_y$centers)
#'     eta_0 <- rep(1/Kinit, Kinit)
#'     sigma2_0 <- array(0, dim = c(1, 1, Kinit))
#'     sigma2_0[1, 1, ] <- 0.5 * C0
#' 
#'     result <- sampleUniNormMixture(
#'         y, S_0, mu_0, sigma2_0, eta_0,
#'         c0, g0, G0, C0, b0, B0,
#'         M, burnin, thin, Kmax,
#'         G, priorOnK, priorOnE0)
#'     
#'     K <- result$K
#'     Kplus <- result$Kplus  
#'     
#'     plot(seq_along(K), K, type = "l", ylim = c(0, max(K)),
#'          xlab = "iteration", main = "",
#'          ylab = expression("K" ~ "/" ~ K["+"]), col = 1)
#'     lines(seq_along(Kplus), Kplus, col = 2)
#'     legend("topright", legend = c("K", expression(K["+"])),
#'            col = 1:2, lty = 1, box.lwd = 0)
#' }
#' 
sampleUniNormMixture <-
    function(y, Z, S, mu, sigma2, eta, beta, c0, g0, G0, C0_0, 
             b0, B0, beta_prior, Omega_prior, M, burnin, thin, Kmax,
             G = c("MixDynamic", "MixStatic"),
             priorOnK, priorOnWeights,
             verbose = FALSE) {
### A effacer en suite
#y <- Y
#S <- S_0
#mu <- alpha_0
#sigma2 <- sigma2_0
#eta <- eta_0
#beta <- beta_0
#C0_0 <- C0
#priorOnK <- priorOn_K
#priorOnWeights <- priorOn_n0
#verbose = FALSE
###

    y <- as.matrix(y)
    ## initial number of components
    K_j <- length(eta)

    ## prior on K and weights
    log_pK <- priorOnK$log_pK
    G <- match.arg(G)
    if (G == "MixDynamic") {
        log_pAlpha <- priorOnWeights$param$log_pAlpha  
        a_alpha <- priorOnWeights$param$a_alpha        
        b_alpha <- priorOnWeights$param$b_alpha        
        alpha <- priorOnWeights$param$alpha
        e0 <- alpha / K_j
    } else {
        e0 <- priorOnWeights$param$e0
        alpha <- e0 * K_j
        log_p_e0 <- priorOnWeights$log_p_e0
    }
    s0_proposal <- priorOnWeights$param$s0_proposal
    
    N <- ncol(y)   # number of cross-section observations
    TT <- nrow(y)  # number of time periods
    r <- 1         # dimension
    
    ## initializing current values
    beta_j <- beta
    invB0 <- chol2inv(chol(B0))     # "chol2inv" Invert a symmetric, positive definite square matrix from 
                                    # its Cholesky decomposition
    eta_j <- eta
    mu_j <- mu
    sigma2_j <- sigma2
    cholsigma2_j <- array(0, dim = c(r, r, K_j))
    invsigma2_j <- array(0, dim = c(r, r, K_j))
    for (k in 1:K_j) {
        cholsigma2_j[, , k] <- chol(sigma2_j[, , k])
        invsigma2_j[, , k] <- chol2inv(cholsigma2_j[, , k])
    }
    S_j <- S
    C0_j <- C0_0
    Nk_j <- tabulate(S_j, K_j)    # Nk_j counts how many observations are per component of the mixture
    if (verbose) {                # When verbose is FALSE or set to a lower level, the program suppresses 
                                  # these detailed messages, producing less output.
        cat("0 ", Nk_j)           # The cat() function concatenates and outputs its arguments
    }
    Kp_j <- sum(Nk_j != 0)        # number of nonempty components
    acc <- FALSE                  # What is "acc"?
    
    ## generating matrices for storing the draws:
    result <- list(Eta = matrix(NA_real_, M, Kmax),
                   Beta = matrix(NA_real_, M, 1),
                   Mu = array(NA_real_, dim = c(M, r, Kmax)),
                   Sigma2 = array(NA_real_, dim = c(M, r, Kmax)),
                   Nk = matrix(NA_real_, M, Kmax),  # It stores the number of observations per component
                   S = matrix(NA_real_, M, N),
                   K = rep(NA_real_, M),            # Number of components in the mixture distribution
                   Kp = rep(0, M),                  # Number of clusters in the data
                   #
                   mixlik = rep(0, M),
                   mixprior = rep(0, M),
                   nonnormpost = rep(0, M),
                   nonnormpost_mode_list = vector("list", Kmax),
                   mixlik_mode_list = vector("list", Kmax),
                   # Hyperparameters:
                   C0 = array(NA_real_, dim = c(M, r, r)),
                   e0 = rep(NA_real_, M),
                   alpha = rep(NA_real_, M),
                   acc = rep(NA_real_, M))

    ## Storing the initial values
    result$Mu[1, , 1:K_j] <- mu_j
    result$Sigma2[1, , 1:K_j] <- sigma2_j
    result$Beta[1,1] <- beta_j
    result$Eta[1, 1:K_j] <- eta_j
    result$S[1, ] <- S_j
    result$K[1] <- K_j
    result$Kp[1] <- Kp_j
    result$Nk[1, 1:K_j] <- Nk_j
    result$e0[1] <- e0
    result$alpha[1] <- alpha
    result$acc[1] <- acc
    for (k in 1:Kmax) {
        result$nonnormpost_mode_list[[k]] <- list(nonnormpost = -(10)^18)
        result$mixlik_mode_list[[k]] <- list(mixlik = -(10)^18)
    }
    
    ##---------------------- simulation ----------------------------------------------
    
    m <- 2
    Mmax <- M * thin
    while (m <= Mmax || m <= burnin) {
        if (verbose && !(m%%500)) {             # 'verbose': checks if the verbose variable is TRUE.
                                                # '!(m%%500)': checks if m %% 500 is equal to 0. 'm %% 500' calculates 
                                                # the remainder when m is divided by 500. 
                                                # The ! operator negates the result. 
                                                # So, !(m %% 500) will be TRUE only when m is a multiple of 500.
            cat("\n", m, " ", Nk_j)             # the 'cat()' function prints a newline character (\n), 
                                                # followed by the current value of m, a space, and the 
                                                # value of Nk_j
        }
        
        if (m == burnin) {
            m <- 1
            burnin <- 0
        }
        
      print(m)
      #cat("m=", m)
        ## FIRST STEP: classify observations and determine new partition
        # (1.a) of Algorithm 2 of Fruhwirth-Schnatter et al. (2021)
        y_betaZ <- y - beta*Z
        mat <- matrix(,N,K_j)
        for (i in 1:N){
          mat[i,] <- sapply(1:K_j, function(k)
            eta_j[k] * mvtnorm::dmvnorm(t(y_betaZ[,i]), rep(mu_j[,k],1,TT), sigma2_j[,,k]*diag(TT)))
        }
        ## If I want to normalise the element of mat use the following commands:
        # mat_sum <- rowSums(mat)
        # mat_norm <- mat/mat_sum
        ## then use mat_norm instead of mat to compute S_j below.

# OLD CODE
#        mat <- sapply(1:K_j, function(k)
#            eta_j[k] * mvtnorm::dmvnorm(y, mu_j[, k] + beta_j*Z, as.matrix(sigma2_j[, , k])))

        S_j <- apply(mat, 1, function(x) sample(1:K_j, 1, prob = x, replace = TRUE))
        
        ## determine partition
        # (1.b) of Algorithm 2 of Fruhwirth-Schnatter et al. (2021)
        Nk_j <- tabulate(S_j, K_j)    # length(Nk_j)=K_j
        Kp_j <- sum(Nk_j != 0)        # Number of cluster
        
        ## reorder the components
        perm <- c(which(Nk_j > 0), which(Nk_j == 0))
        mu_j <- mu_j[, perm, drop = FALSE]  # length(mu_j[1, ])=K_j
        sigma2_j <- sigma2_j[, , perm, drop = FALSE]  # length(sigma_j[1, 1, ])=K_j
        S_ <- rep(FALSE, N)
        for (i in 1:length(perm)) {               # length(perm) = K_j
            S_[S_j == i] <- which(perm == i)
        }
        S_j <- S_
        Nk_j <- tabulate(S_j, Kp_j)  # length(Nk_j) = Kp_j
        
        ## SECOND STEP: update parameters conditional on partition C=(N_1, ..., N_K+):
        
        ## (2.a) update parameters of filled components (i): sample Sigma^{-} for filled components
        Ck <- array(0, dim = c(r, r, Kp_j))
        ck <- c0 + TT * Nk_j/2
        for (k in 1:Kp_j) {
            #Ck[, , k] <- C0_j + 0.5 * sum((y[,S_j == k] - mu_j[, k])^2)
            Ck[, , k] <- C0_j + 0.5 * sum((y_betaZ[,S_j == k] - mu_j[, k])^2)
            invsigma2_j[, , k] <- rgamma(1, shape = ck[k], rate = Ck[, , k])
            sigma2_j[, , k] <- 1/invsigma2_j[k]
        }
        
        ## (i): Sample mu_j for filled components
        mean_yk <- matrix(0, nrow = 1, ncol = Kp_j)
        for (t in 1:TT){
          mean_yk_t <- matrix(sapply(1:Kp_j, function(k) rowMeans(y_betaZ[t,S_j == k, drop = FALSE])), nrow = 1, ncol = Kp_j)
          mean_yk <- mean_yk + mean_yk_t
        }
        # mean_yk <- matrix(sapply(1:Kp_j, function(k) colMeans(y[S_j == k, , drop = FALSE])), ncol = Kp_j)
        Bk <- array(0, dim = c(r, r, Kp_j))
        bk <- matrix(0, r, Kp_j)
        for (k in 1:Kp_j) {
            Bk[, , k] <- 1/(invB0 + TT*Nk_j[k] * invsigma2_j[, , k])
            bk[, k] <- Bk[, , k] * (invB0 * b0 + invsigma2_j[, , k] * mean_yk[, k] * Nk_j[k])
            mu_j[, k] <- rnorm(1, mean = bk[, k], sd = sqrt(Bk[, , k]))
        }
        
        ## (2.b) sample hyperparameters conditional on P (i): sample C0
        gK <- g0 + Kp_j * c0            # Kp_j is the number of components of the mixture
        C0_j <- rgamma(1, shape = gK, rate = G0 + sum(invsigma2_j[, , 1:Kp_j]))
        
        ## (2.c): sample beta from the conditional posterior given (S, K_+ and theta)
        # This is an extra step that is not in the original Fruhwirth-Schnatter et al. (2021) algorithm.
        mu_S_j <- matrix(,N,1)
        sigma2_S_j <- matrix(,N,1)
        for (i in 1:N){
          mu_S_j[i,] <- mu_j[,S_j[i]]
          sigma2_S_j[i,] <- sigma2_j[,,S_j[i]]
        }
        # Posterior variance of the conditional posterior of beta
        Omega_post <- 1/(1/Omega_prior + sum(colSums(Z^2)/sigma2_S_j))
        # Posterior mean of the conditional posterior of beta
        beta_post <- Omega_post*(beta_prior/Omega_prior + sum(colSums((y - rep(t(mu_S_j),TT,1))*Z)/sigma2_S_j) )
        # Draw a new beta:
        beta_j <- rnorm(1, mean = beta_post, sd = sqrt(Omega_post))
        
        
        ## THIRD STEP: sample K and alpha (or e0) conditional on partition
        
        if (G == "MixDynamic") {
            ## (3a) Sample K, if e0=alpha/K (=dependent on K)
            K_j <- sampleK_alpha(Kp_j, Kmax, Nk_j, alpha, log_pK)
            
            ## (3b) Sample alpha, if alpha~p(a_alpha, b_alpha)
            value <- sampleAlpha(N, Nk_j, K_j, alpha, s0_proposal, log_pAlpha)
            alpha <- value$alpha
            e0 <- alpha / K_j
            acc <- value$acc
        } else {
            ## (3a*) Sample K, if e0 fixed or e0~G(a_e, b_e) (independent of K):
            K_j <- sampleK_e0(Kp_j, Kmax, log_pK, log_p_e0, e0, N)
            
            ## (3b*) Sample e0, if e0~G(a_e, b_e) (independent of K)
            value <- sampleE0(K_j, Kp_j, N, Nk_j, s0_proposal, e0, log_p_e0)
            e0 <- value$e0
            alpha <- e0 * K_j
            acc <- value$acc
        }
        
        ## FOURTH step: add empty components conditional on K
        
        ## (4a) Add/remove empty components
        if (K_j > Kp_j) {
            Nk_j <- c(Nk_j[1:Kp_j], rep(0, (K_j - Kp_j)))  # length(Nk_j)=K_j
            invsigma2_e <- array(NA_real_, dim = c(r, r, K_j - Kp_j))
            invsigma2_j <- abind::abind(invsigma2_j[, , 1:Kp_j, drop = FALSE], invsigma2_e, along = 3)
            
            sigma2_e <- array(NA_real_, dim = c(r, r, K_j - Kp_j))
            sigma2_j <- abind::abind(sigma2_j[, , 1:Kp_j, drop = FALSE], sigma2_e, along = 3)
            
            mu_j <- cbind(mu_j[, 1:Kp_j, drop = FALSE], matrix(0, r, K_j - Kp_j))
            
            for (k in (Kp_j + 1):K_j) {
                invsigma2_j[, , k] <- rgamma(1, shape = c0, rate = C0_j)
                sigma2_j[, , k] <- 1/invsigma2_j[, , k]
                mu_j[, k] <- rnorm(1, mean = b0, sd = sqrt(B0))
            }
        } else {
            invsigma2_j <- invsigma2_j[, , 1:K_j, drop = FALSE]
            sigma2_j <- sigma2_j[, , 1:K_j, drop = FALSE]
            mu_j <- mu_j[, 1:K_j, drop = FALSE]
        }
        
        ## (4b): Sample eta_j:
        ek <- e0 + Nk_j
        eta_j <- bayesm::rdirichlet(ek)
        
        ## FIFTH STEP: evaluating the mixture likelihood and storing the current values
        
        ## evaluating the mixture likelihood:
        y_betaZ <- y - beta*Z
        mat_neu <- matrix(,N,K_j)
        for (i in 1:N){
          mat_neu[i,] <- sapply(1:K_j, function(k)
            eta_j[k] * mvtnorm::dmvnorm(t(y_betaZ[,i]), rep(mu_j[,k],1,TT), sigma2_j[,,k]*diag(TT)))
        }
        ## OLD CODE
        # mat_neu <- sapply(1:K_j, function(k)
        #     eta_j[k] * mvtnorm::dmvnorm(y, mu_j[, k], as.matrix(sigma2_j[, , k])))
        mixlik_j <- sum(log(rowSums(mat_neu)))
  
        ## evaluating the mixture prior
        mixprior_j <- log(MCMCpack::ddirichlet(as.vector(eta_j), rep(e0, K_j))) +
            sum(mvtnorm::dmvnorm(t(mu_j), b0, B0, log = TRUE)) +
            sum(dgamma(invsigma2_j, shape = c0, rate = C0_j, log = TRUE)) +
            dgamma(C0_j, shape = g0, rate = G0) + log_pK(K_j)
        
        if ((burnin == 0) && !(m%%thin)) {
            result$mixlik[m] <- mixlik_j
            result$mixprior[m] <- mixprior_j
            result$nonnormpost[m] <- result$mixlik[m] + result$mixprior[m]
        }

        ## storing the nonnormalized posterior for having good
        ## starting points when clustering the draws in the point
        ## process repres.
        if (((burnin == 0) && !(m%%thin)) && (result$nonnormpost[m] > result$nonnormpost_mode_list[[Kp_j]]$nonnormpost)) {
            result$nonnormpost_mode_list[[Kp_j]] <- list(nonnormpost = result$nonnormpost[m],
                                                         mu = mu_j[, Nk_j != 0],
                                                         S = S_j,
                                                         bk = bk,
                                                         Bk = Bk,
                                                         eta = eta_j,
                                                         sigma2_j = sigma2_j[, , Nk_j != 0],
                                                         ck = ck,
                                                         Ck = Ck)
        }
        
        ## intermediate step: storing the results for given S_j, Nk, K_j
        if ((burnin == 0) && !(m%%thin)) {
            ## storing the new values
            result$Mu[m/thin, , 1:K_j] <- mu_j
            result$Sigma2[m/thin, , 1:K_j] <- sigma2_j
            result$Beta[m/thin,1] <- beta_j
            result$Eta[m/thin, 1:K_j] <- eta_j
            result$S[m/thin, ] <- S_j
            result$Nk[m/thin, 1:K_j] <- Nk_j
            result$K[m/thin] <- K_j
            result$Kp[m/thin] <- Kp_j
            result$C0[m/thin, , ] <- C0_j
            result$e0[m/thin] <- e0
            result$alpha[m/thin] <- alpha
            result$acc[m/thin] <- acc
        }
        m <- m + 1
    }
    return(result)
}

