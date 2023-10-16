print.spemix <- function(x, ...) {
  #print fitted object
  cat(ncol(x$bestmod$z), "groups and scale structure", x$bestmod$gpar$model, "was selected by the BIC.", "\n", "BIC was", x$BIC[x$bestBIC[1],x$bestBIC[2]])
}
marginal_m_l <- function(my, mu, Sig, bval, p, log.d = FALSE) {
  #needed for calculating density of SPE
  val <- log(bval) + log(gamma(p/2)/pi^(p/2)) + (-1/2 * mahalanobis(my, center = mu, cov = Sig,
                                                                    inverted = FALSE)^bval) - log((2^(p/(2 * bval)) * gamma(p/(2 * bval)))) - 1/2 *
    log(det(Sig))
  if (!log.d) {
    return(exp(val))
  } else {
    return(val)
  }
}
cov_pe <- function(scale = NULL, beta = NULL, p = NULL) {
  #covariance of (elliptical) power exponential distribution; used to simulate SPE distributions
  Cov <- 2^(1/beta) * gamma((p + 2)/(2 * beta))/(p * gamma(p/(2 * beta))) * scale
  return(Cov)
}

mpower <- function(A, p = 1) {
  #matrix power
  E = eigen(A)
  M = E$vectors %*% diag(abs(E$values)^p) %*% t(E$vectors)
  return(M)
}

marginal_m_l <- function(my, mu, Sig, bval, p, log.d = FALSE) {
  #needed for calculating density of SPE
  val <- log(bval) + log(gamma(p/2)/pi^(p/2)) + (-1/2 * mahalanobis(my, center = mu, cov = Sig,
                                                                    inverted = FALSE)^bval) - log((2^(p/(2 * bval)) * gamma(p/(2 * bval)))) - 1/2 *
    log(det(Sig))
  if (!log.d) {
    return(exp(val))
  } else {
    return(val)
  }
}

dspe <- function(y, mu = rep(0, nrow(Sig)), Sig = diag(length(mu)), psi = rep(0, nrow(Sig)), bval = 1, p = nrow(Sig)) {
  #density of SPE
  marg <- marginal_m_l(y, mu, Sig, bval, p, log.d = TRUE)
  cond <- pnorm(as.numeric(t(psi) %*% mpower(Sig, p = -1/2) %*% (y - mu)), log.p = TRUE)
  return(list(dens = exp(log(2) + marg + cond), marg = exp(marg), cond = exp(cond)))
}
rspe <- function(n, location = rep(0, nrow(scale)), scale = diag(length(location)), beta = 1, psi = c(0, 0)) {
  #generate SPE distribution (random numbers)
  cc = 0
  i = 1
  p = length(location)
  cov_pe <- function(scale = NULL, beta = NULL, p = NULL) {
    Cov <- 2^(1/beta) * gamma((p + 2)/(2 * beta))/(p * gamma(p/(2 * beta))) * scale
    return(Cov)
  }
  psi = psi
  mu = location
  Sig = scale
  bval = beta
  cov <- cov_pe(scale = Sig, beta = beta, p = p) #get the corresponding variance for the density proposal we are using by getting the variance of the marginal distribution to which the truncation is being applied to.

  z = matrix(NA, (n + 2001), p)
  z[i, ] = location - rep(0.2, p)
  h = function(x) {
    mvtnorm::dmvnorm(x, mu, cov * 2) #multiply by 2
  }

  while (cc < (n + 2000)) {
    y = c(mvtnorm::rmvnorm(1, mu, cov * 2)) #multiply by 2
    alpha = min(1, (dspe(y, mu, Sig, psi, bval, p)$dens * h(z[i, ]))/(dspe(z[i, ], mu, Sig, psi, bval, p)$dens * h(y)))
    U = runif(1)
    if (U < alpha) {
      z[i + 1, ] <- y
      i <- i + 1
      cc <- cc + 1
    } else z[i, ] <- z[i, ]
  }
  return(z[-c(1:2001), ])
}

rpe <- function(n = NULL, beta = NULL, mean = NULL, scale = NULL){
  runifsphere<-function (n, p)
  {
    p <- as.integer(p)
    if (!is.integer(p)) stop("p must be an integer larger or equal than 2")
    if (p < 2) stop("p must be an integer larger or equal than 2")
    Mnormal <- matrix(rnorm(n * p, 0, 1), nrow = n)
    rownorms <- sqrt(rowSums(Mnormal^2))
    unifsphere <- sweep(Mnormal, 1, rownorms, "/")
    return(unifsphere)
  }
  rmvpex<-function (n, mean = rep(0, nrow(scale)), scale = diag(length(mean)), beta = 1)
  {
    p <- length(mean)
    if (!isSymmetric(scale, tol = sqrt(.Machine$double.eps))) stop("scale must be a symmetric matrix")
    if (p != nrow(scale)) stop("mean and scale have non-conforming size")
    ev <- eigen(scale, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("scale is numerically not positive definite")
    }
    scaleSqrt <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*%
      t(ev$vectors)
    radius <- (rgamma(n, shape = p/(2 * beta), scale = 2))^(1/(2 * beta))
    un <- runifsphere(n = n, p = p)
    mvpowerexp <- radius * un %*% scaleSqrt
    mvpowerexp <- sweep(mvpowerexp, 2, mean, "+")
    return(mvpowerexp)
  }
  return(rmvpex(n, mean = mean, scale = scale, beta = beta))
}
