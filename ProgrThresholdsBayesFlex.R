

### discrete and continuous

fit.thresholdFlex <- function(data,
                          F,
                          f,
                          g,
                          gprime = NULL,
                          q = 1,
                          dalpha,
                          alpha.par = NULL,
                          mu.delta0 = 0,
                          sd.delta0 = 2,
                          mu.lambda = 0,
                          sd.lambda = 1,
                          niter = 10000,
                          burnin = 2000,
                          thin = 1,
                          proposal.sd = 0.15,
                          start = NULL,
                          discrete = FALSE,
                          verbose = TRUE) {
  
  ## ------------------------------------------------------------
  ## Data
  ## ------------------------------------------------------------
  
  Y <- as.matrix(data)
  
  P <- nrow(Y)
  I <- ncol(Y)
  
  if (length(g) != I)
    stop("g must contain one function for each item.")
  
  if (!discrete && length(gprime) != I)
    stop("gprime must contain one function for each item.")
  
  if (!is.function(F))
    stop("F must be a function.")
  
  if (!is.function(f))
    stop("f must be the derivative of F.")
  
  if (!is.function(dalpha))
    stop("dalpha must be a density function for alpha.")
  
  if (q < 1 || q != as.integer(q))
    stop("q must be a positive integer.")
  
  
  ## ------------------------------------------------------------
  ## Checks for discrete responses
  ## ------------------------------------------------------------
  
  if (discrete) {
    
    ## All non-missing responses must be
    ## non-negative integers
    
    if (any(!is.na(Y) &
            (Y < 0 | Y != floor(Y))))
      stop("For discrete = TRUE, responses must be non-negative integers.")
    
    ## Check g at all observed response values
    ymax <- max(Y, na.rm = TRUE)
    
    for (i in seq_len(I)) {
      
      if (!is.function(g[[i]]))
        stop("Each element of g must be a function.")
      
      gy <- sapply(0:ymax, g[[i]])
      
      if (any(!is.finite(gy)))
        stop("g must be finite for all observed response categories.")
      
      ## g must be strictly increasing over the observed categories
      if (length(gy) > 1 && any(diff(gy) <= 0))
        stop("g must be strictly increasing for discrete responses.")
    }
  }
  
  
  ## ------------------------------------------------------------
  ## Starting values
  ## ------------------------------------------------------------
  
  if (is.null(start)) {
    
    theta <- matrix(rnorm(P * q), P, q)
    
    alpha <- matrix(rnorm(I * q), I, q)
    
    delta0 <- rep(0, I)
    
    lambda <- rep(0, I)
    
  } else {
    
    theta  <- as.matrix(start$theta)
    alpha  <- as.matrix(start$alpha)
    delta0 <- as.numeric(start$delta0)
    lambda <- as.numeric(start$lambda)
    
    if (!all(dim(theta) == c(P, q)))
      stop("start$theta must have dimensions P x q.")
    
    if (!all(dim(alpha) == c(I, q)))
      stop("start$alpha must have dimensions I x q.")
    
    if (length(delta0) != I)
      stop("start$delta0 must have length equal to the number of items.")
    
    if (length(lambda) != I)
      stop("start$lambda must have length equal to the number of items.")
  }
  
  
  ## ------------------------------------------------------------
  ## Log-likelihood
  ## ------------------------------------------------------------
  
  loglik <- function(theta, alpha, delta0, lambda) {
    
    delta <- exp(lambda)
    
    ll <- 0
    
    for (p in seq_len(P)) {
      
      for (i in seq_len(I)) {
        
        y <- Y[p, i]
        
        if (is.na(y))
          next
        
        
        ## ========================================================
        ## DISCRETE RESPONSE
        ## ========================================================
        
        if (discrete) {
          
          eta <- sum(theta[p, ] * alpha[i, ])
          
          
          ## ------------------------------------------------------
          ## y = 0
          ## P(Y = 0) = 1 - P(Y > 0)
          ## ------------------------------------------------------
          
          if (y == 0) {
            
            gy <- g[[i]](0)
            
            z <- eta -
              delta0[i] -
              delta[i] * gy
            
            Fy <- F(z)
            
            prob <- 1 - Fy
            
            
            ## ------------------------------------------------------
            ## y >= 1
            ##
            ## P(Y = y)
            ## =
            ## P(Y > y-1) - P(Y > y)
            ## ------------------------------------------------------
            
          } else {
            
            gy.lower <- g[[i]](y - 1)
            gy.upper <- g[[i]](y)
            
            z.lower <- eta -
              delta0[i] -
              delta[i] * gy.lower
            
            z.upper <- eta -
              delta0[i] -
              delta[i] * gy.upper
            
            F.lower <- F(z.lower)
            F.upper <- F(z.upper)
            
            prob <- F.lower - F.upper
          }
          
          
          ## ------------------------------------------------------
          ## Probability check
          ## ------------------------------------------------------
          
          if (!is.finite(prob) || prob <= 0)
            return(-Inf)
          
          ll <- ll + log(prob)
          
          
          ## ========================================================
          ## CONTINUOUS RESPONSE
          ## ========================================================
          
        } else {
          
          gy <- g[[i]](y)
          
          gpy <- gprime[[i]](y)
          
          if (!is.finite(gy) ||
              !is.finite(gpy) ||
              gpy <= 0)
            return(-Inf)
          
          eta <- sum(theta[p, ] * alpha[i, ])
          
          z <- eta -
            delta0[i] -
            delta[i] * gy
          
          lf <- log(f(z))
          
          if (!is.finite(lf))
            return(-Inf)
          
          logfy <- lambda[i] +
            log(gpy) +
            lf
          
          if (!is.finite(logfy))
            return(-Inf)
          
          ll <- ll + logfy
        }
      }
    }
    
    ll
  }
  
  
  ## ------------------------------------------------------------
  ## Log-prior
  ## ------------------------------------------------------------
  
  logprior <- function(theta, alpha, delta0, lambda) {
    
    lp <- 0
    
    ## theta_p ~ N_q(0,I)
    lp <- lp +
      sum(dnorm(theta,
                mean = 0,
                sd = 1,
                log = TRUE))
    
    
    ## User-specified distribution of alpha
    for (i in seq_len(I)) {
      
      da <- dalpha(alpha[i, ], alpha.par)
      
      if (!is.finite(da) || da <= 0)
        return(-Inf)
      
      lp <- lp + log(da)
    }
    
    
    ## delta0
    lp <- lp +
      sum(dnorm(delta0,
                mean = mu.delta0,
                sd = sd.delta0,
                log = TRUE))
    
    
    ## lambda = log(delta)
    lp <- lp +
      sum(dnorm(lambda,
                mean = mu.lambda,
                sd = sd.lambda,
                log = TRUE))
    
    lp
  }
  
  
  ## ------------------------------------------------------------
  ## Log-posterior
  ## ------------------------------------------------------------
  
  logpost <- function(theta, alpha, delta0, lambda) {
    
    ll <- loglik(theta,
                 alpha,
                 delta0,
                 lambda)
    
    if (!is.finite(ll))
      return(-Inf)
    
    lp <- logprior(theta,
                   alpha,
                   delta0,
                   lambda)
    
    if (!is.finite(lp))
      return(-Inf)
    
    ll + lp
  }
  
  
  ## ------------------------------------------------------------
  ## Initial value
  ## ------------------------------------------------------------
  
  current <- logpost(theta,
                     alpha,
                     delta0,
                     lambda)
  
  if (!is.finite(current))
    stop("Initial parameter values have non-finite posterior.")
  
  
  ## ------------------------------------------------------------
  ## Storage
  ## ------------------------------------------------------------
  
  nkeep <- floor((niter - burnin) / thin)
  
  if (nkeep <= 0)
    stop("niter must be larger than burnin.")
  
  theta.save <- array(
    NA_real_,
    dim = c(nkeep, P, q)
  )
  
  alpha.save <- array(
    NA_real_,
    dim = c(nkeep, I, q)
  )
  
  delta0.save <- matrix(
    NA_real_,
    nrow = nkeep,
    ncol = I
  )
  
  lambda.save <- matrix(
    NA_real_,
    nrow = nkeep,
    ncol = I
  )
  
  accept <- 0
  keep <- 0
  
  
  ## ------------------------------------------------------------
  ## MCMC
  ## ------------------------------------------------------------
  
  for (iter in seq_len(niter)) {
    
    ## ----------------------------------------------------------
    ## theta
    ## ----------------------------------------------------------
    
    for (p in seq_len(P)) {
      
      for (j in seq_len(q)) {
        
        proposal <- theta
        
        proposal[p, j] <-
          proposal[p, j] +
          rnorm(1, 0, proposal.sd)
        
        proposed <-
          logpost(proposal,
                  alpha,
                  delta0,
                  lambda)
        
        if (log(runif(1)) < proposed - current) {
          
          theta <- proposal
          current <- proposed
          accept <- accept + 1
        }
      }
    }
    
    
    ## ----------------------------------------------------------
    ## alpha
    ## ----------------------------------------------------------
    
    for (i in seq_len(I)) {
      
      for (j in seq_len(q)) {
        
        proposal <- alpha
        
        proposal[i, j] <-
          proposal[i, j] +
          rnorm(1, 0, proposal.sd)
        
        proposed <-
          logpost(theta,
                  proposal,
                  delta0,
                  lambda)
        
        if (log(runif(1)) < proposed - current) {
          
          alpha <- proposal
          current <- proposed
          accept <- accept + 1
        }
      }
    }
    
    
    ## ----------------------------------------------------------
    ## delta0
    ## ----------------------------------------------------------
    
    for (i in seq_len(I)) {
      
      proposal <- delta0
      
      proposal[i] <-
        proposal[i] +
        rnorm(1, 0, proposal.sd)
      
      proposed <-
        logpost(theta,
                alpha,
                proposal,
                lambda)
      
      if (log(runif(1)) < proposed - current) {
        
        delta0 <- proposal
        current <- proposed
        accept <- accept + 1
      }
    }
    
    
    ## ----------------------------------------------------------
    ## lambda
    ## ----------------------------------------------------------
    
    for (i in seq_len(I)) {
      
      proposal <- lambda
      
      proposal[i] <-
        proposal[i] +
        rnorm(1, 0, proposal.sd)
      
      proposed <-
        logpost(theta,
                alpha,
                delta0,
                proposal)
      
      if (log(runif(1)) < proposed - current) {
        
        lambda <- proposal
        current <- proposed
        accept <- accept + 1
      }
    }
    
    
    ## ----------------------------------------------------------
    ## Save
    ## ----------------------------------------------------------
    
    if (iter > burnin &&
        (iter - burnin) %% thin == 0) {
      
      keep <- keep + 1
      
      theta.save[keep, , ] <-
        array(theta, dim = c(1, P, q))
      
      alpha.save[keep, , ] <-
        array(alpha, dim = c(1, I, q))
      
      delta0.save[keep, ] <- delta0
      
      lambda.save[keep, ] <- lambda
    }
    
    
    ## ----------------------------------------------------------
    ## Progress
    ## ----------------------------------------------------------
    
    if (verbose && iter %% 1000 == 0)
      cat("Iteration:", iter, "\n")
  }
  
  
  ## ------------------------------------------------------------
  ## Posterior means
  ## ------------------------------------------------------------
  
  theta.hat <- apply(
    theta.save,
    c(2, 3),
    mean
  )
  
  alpha.hat <- apply(
    alpha.save,
    c(2, 3),
    mean
  )
  
  delta0.hat <- colMeans(delta0.save)
  
  lambda.hat <- colMeans(lambda.save)
  
  delta.hat <- colMeans(exp(lambda.save))
  
  
  ## ------------------------------------------------------------
  ## Acceptance rate
  ## ------------------------------------------------------------
  
  acceptance.rate <-
    accept /
    (niter * (P * q + I * q + I + I))
  
  
  ## ------------------------------------------------------------
  ## Return
  ## ------------------------------------------------------------
  
  list(
    
    estimates = list(
      theta = theta.hat,
      alpha = alpha.hat,
      delta0 = delta0.hat,
      lambda = lambda.hat,
      delta = delta.hat
    ),
    
    samples = list(
      theta = theta.save,
      alpha = alpha.save,
      delta0 = delta0.save,
      lambda = lambda.save
    ),
    
    acceptance.rate = acceptance.rate,
    
    call = match.call()
  )
}

##############################################

loglik.thresholdFlex <- function(data,
                             theta,
                             alpha,
                             delta0,
                             lambda,
                             F,
                             f,
                             g,
                             gprime = NULL,
                             discrete = FALSE) {
  
  Y <- as.matrix(data)
  
  P <- nrow(Y)
  I <- ncol(Y)
  
  theta <- as.matrix(theta)
  alpha <- as.matrix(alpha)
  
  if (nrow(theta) != P)
    stop("Incorrect dimensions of theta.")
  
  if (nrow(alpha) != I)
    stop("alpha must have one row for each item.")
  
  if (length(delta0) != I)
    stop("delta0 must have length equal to the number of items.")
  
  if (length(lambda) != I)
    stop("lambda must have length equal to the number of items.")
  
  if (length(g) != I)
    stop("g must contain one function for each item.")
  
  delta <- exp(lambda)
  
  ll <- 0
  
  for (p in seq_len(P)) {
    
    for (i in seq_len(I)) {
      
      y <- Y[p, i]
      
      if (is.na(y))
        next
      
      eta <- sum(theta[p, ] * alpha[i, ])
      
      if (discrete) {
        
        ## -------------------------------------------------------
        ## Discrete nonnegative integer response
        ## -------------------------------------------------------
        
        if (!is.finite(y) ||
            y < 0 ||
            y != floor(y))
          return(-Inf)
        
        if (y == 0) {
          
          gy <- g[[i]](0)
          
          if (!is.finite(gy))
            return(-Inf)
          
          z <- eta - delta0[i] - delta[i] * gy
          
          prob <- 1 - F(z)
          
        } else {
          
          gy0 <- g[[i]](y - 1)
          gy1 <- g[[i]](y)
          
          if (!is.finite(gy0) ||
              !is.finite(gy1) ||
              gy1 <= gy0)
            return(-Inf)
          
          z0 <- eta - delta0[i] - delta[i] * gy0
          z1 <- eta - delta0[i] - delta[i] * gy1
          
          prob <- F(z0) - F(z1)
        }
        
        if (!is.finite(prob) || prob <= 0)
          return(-Inf)
        
        ll <- ll + log(prob)
        
      } else {
        
        ## -------------------------------------------------------
        ## Continuous response
        ## -------------------------------------------------------
        
        if (length(gprime) != I)
          stop("gprime must contain one function for each item.")
        
        gy <- g[[i]](y)
        gpy <- gprime[[i]](y)
        
        if (!is.finite(gy) ||
            !is.finite(gpy) ||
            gpy <= 0)
          return(-Inf)
        
        z <- eta - delta0[i] - delta[i] * gy
        
        lf <- log(f(z))
        
        if (!is.finite(lf))
          return(-Inf)
        
        logfy <- lambda[i] + log(gpy) + lf
        
        if (!is.finite(logfy))
          return(-Inf)
        
        ll <- ll + logfy
      }
    }
  }
  
  return(ll)
}

################################################
fit.threshold.binary <- function(data, F,
                                 q = 1,
                                 dalpha,
                                 alpha.par = NULL,
                                 mu.delta0 = 0,
                                 sd.delta0 = 2,
                                 niter = 10000,
                                 burnin = 2000,
                                 thin = 1,
                                 proposal.sd = 0.15,
                                 start = NULL,
                                 verbose = TRUE) {
  
  Y <- as.matrix(data)
  P <- nrow(Y)
  I <- ncol(Y)
  
  ## Check binary responses
  if (!all(Y %in% c(0, 1))) {
    stop("data must contain only 0 and 1")
  }
  
  ## --------------------------------------------------
  ## Initial values
  ## --------------------------------------------------
  
  if (is.null(start)) {
    
    ## theta ~ N(0,I)
    theta <- matrix(rnorm(P * q), P, q)
    
    ## alpha
    if (!is.null(alpha.par)) {
      
      alpha <- matrix(NA_real_, I, q)
      
      for (i in 1:I) {
        alpha[i, ] <- MASS::mvrnorm(
          1,
          mu = alpha.par$mu,
          Sigma = alpha.par$Sigma
        )
      }
      
    } else {
      
      ## Reasonable starting values
      alpha <- matrix(1, I, q)
      
    }
    
    ## delta0
    delta0 <- rnorm(I, mu.delta0, sd.delta0)
    
  } else {
    
    theta <- start$theta
    alpha <- start$alpha
    delta0 <- start$delta0
  }
  
  
  ## --------------------------------------------------
  ## Log likelihood
  ## --------------------------------------------------
  
  loglik <- function(theta, alpha, delta0) {
    
    eta <- theta %*% t(alpha)
    
    eta <- sweep(eta, 2, delta0, "-")
    
    prob <- F(eta)
    
    prob <- pmin(pmax(prob, 1e-12), 1 - 1e-12)
    
    sum(
      Y * log(prob) +
        (1 - Y) * log(1 - prob)
    )
  }
  
  
  ## --------------------------------------------------
  ## Log priors
  ## --------------------------------------------------
  
  logprior.theta <- function(theta) {
    
    ## theta_p ~ N_q(0,I_q)
    -0.5 * sum(theta^2)
  }
  
  
  logprior.alpha <- function(alpha) {
    
    if (is.null(alpha.par)) {
      
      ## No alpha prior specified
      0
      
    } else {
      
      val <- 0
      
      for (i in 1:I) {
        
        val <- val +
          mvtnorm::dmvnorm(
            alpha[i, ],
            mean = alpha.par$mu,
            sigma = alpha.par$Sigma,
            log = TRUE
          )
      }
      
      val
    }
  }
  
  
  logprior.delta0 <- function(delta0) {
    
    sum(
      dnorm(
        delta0,
        mean = mu.delta0,
        sd = sd.delta0,
        log = TRUE
      )
    )
  }
  
  
  ## --------------------------------------------------
  ## Complete log posterior
  ## --------------------------------------------------
  
  logpost <- function(theta, alpha, delta0) {
    
    loglik(theta, alpha, delta0) +
      logprior.theta(theta) +
      logprior.alpha(alpha) +
      logprior.delta0(delta0)
  }
  
  
  ## --------------------------------------------------
  ## Storage
  ## --------------------------------------------------
  
  nkeep <- floor((niter - burnin) / thin)
  
  theta.save <- array(
    NA_real_,
    dim = c(nkeep, P, q)
  )
  
  alpha.save <- array(
    NA_real_,
    dim = c(nkeep, I, q)
  )
  
  delta0.save <- matrix(
    NA_real_,
    nkeep,
    I
  )
  
  logpost.save <- numeric(nkeep)
  
  
  ## Current posterior
  lp <- logpost(theta, alpha, delta0)
  
  keep <- 0
  
  
  ## --------------------------------------------------
  ## MCMC
  ## --------------------------------------------------
  
  for (iter in 1:niter) {
    
    ## ----------------------------------------------
    ## Update theta
    ## ----------------------------------------------
    
    for (p in 1:P) {
      
      theta.old <- theta[p, ]
      
      theta.new <- theta.old +
        rnorm(q, 0, proposal.sd)
      
      theta[p, ] <- theta.new
      
      lp.new <- logpost(theta, alpha, delta0)
      
      if (log(runif(1)) < lp.new - lp) {
        lp <- lp.new
      } else {
        theta[p, ] <- theta.old
      }
    }
    
    
    ## ----------------------------------------------
    ## Update alpha
    ## ----------------------------------------------
    
    for (i in 1:I) {
      
      alpha.old <- alpha[i, ]
      
      alpha.new <- alpha.old +
        rnorm(q, 0, proposal.sd)
      
      alpha[i, ] <- alpha.new
      
      lp.new <- logpost(theta, alpha, delta0)
      
      if (log(runif(1)) < lp.new - lp) {
        lp <- lp.new
      } else {
        alpha[i, ] <- alpha.old
      }
    }
    
    
    ## ----------------------------------------------
    ## Update delta0
    ## ----------------------------------------------
    
    for (i in 1:I) {
      
      delta0.old <- delta0[i]
      
      delta0.new <- delta0.old +
        rnorm(1, 0, proposal.sd)
      
      delta0[i] <- delta0.new
      
      lp.new <- logpost(theta, alpha, delta0)
      
      if (log(runif(1)) < lp.new - lp) {
        lp <- lp.new
      } else {
        delta0[i] <- delta0.old
      }
    }
    
    
    ## ----------------------------------------------
    ## Store
    ## ----------------------------------------------
    
    if (iter > burnin &&
        ((iter - burnin) %% thin == 0)) {
      
      keep <- keep + 1
      
      theta.save[keep, , ] <- theta
      alpha.save[keep, , ] <- alpha
      delta0.save[keep, ] <- delta0
      logpost.save[keep] <- lp
    }
    
    
    if (verbose && iter %% 1000 == 0) {
      
      cat(
        "Iteration:", iter,
        " log-posterior:", lp, "\n"
      )
    }
  }
  
  
  ## --------------------------------------------------
  ## Posterior means
  ## --------------------------------------------------
  
  theta.mean <- apply(theta.save, c(2, 3), mean)
  alpha.mean <- apply(alpha.save, c(2, 3), mean)
  delta0.mean <- colMeans(delta0.save)
  
  
  ## --------------------------------------------------
  ## Return
  ## --------------------------------------------------
  
  list(
    theta = theta.save,
    alpha = alpha.save,
    delta0 = delta0.save,
    theta.mean = theta.mean,
    alpha.mean = alpha.mean,
    delta0.mean = delta0.mean,
    logpost = logpost.save
  )
 ##################
  list(
    estimates = list(
      theta = theta.mean,
      alpha = alpha.mean,
      delta0 = delta0.mean
    ),
    samples = list(
      theta = theta.save,
      alpha = alpha.save,
      delta0 = delta0.save
    ),
    logpost = logpost.save
  )
  
  }