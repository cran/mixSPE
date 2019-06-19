print.spemix <- function(x, ...) {
  #print fitted object
  cat("Call:\n")
  print.default(x$call)
  cat("-----------------------------------\n")
  cat(x$msc[which.max(x$msc[, 2]), 1], "groups and scale structure", attr(which.max(x$msc[,
                                                                                          2]), "names"), "was selected by the BIC.", "\n", "BIC was", x$msc[which.max(x$msc[,
                                                                                                                                                                            2]), 2], "\n", x$msc[which.max(x$msc[, 2]), 3], "parameters were estimated.\n",
      "Log-likelihood calculated at MLEs was", x$msc[which.max(x$msc[, 2]), 4], "\n")
  cat(x$msc[which.max(x$msc[, 5]), 1], "groups and scale structure", attr(which.max(x$msc[,
                                                                                          5]), "names"), "was selected by the ICL.", "\n", "ICL was", x$msc[which.max(x$msc[,
                                                                                                                                                                            5]), 5], "\n", x$msc[which.max(x$msc[, 5]), 3], "parameters were estimated.\n",
      "Log-likelihood calculated at MLEs was", x$msc[which.max(x$msc[, 5]), 4], "\n")
  cat("-----------------------------------\n")
  cat("Time taken:", x$time[3], "seconds.\n")
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


mspe <- function(verbose = FALSE, dat = NULL, seedno = 1, G = 1:4, start = "kmeans", kmeansinit=10, eps=0.005, maxit=2000, anneal = NULL, label=NULL, psistart="zero", modelnames = c("EIIE", "VIIE", "EEIE", "VVIE", "EEEE", "EEVE", "VVEE", "VVVE", "EIIV", "VIIV", "EEIV", "VVIV", "EEEV", "EEVV", "VVEV", "VVVV")) {
  #main workhorse function; fits mixture of SPE to data.
  #dat is data, seedno fixes the seed---important for initialization, esp. random,
  #start controls what kind of start: random or kmeans
  #eps is used in the Aitken's stopping criterion
  #maxit is the maximum number of GEM iterations allowed
  #anneal controls the annealing
  #label gives the already labelled observations
  #psistart=est means using an initial estimated skewness vector rather than zeros.
  ptm <- proc.time()
  if (is.null(dat))
    stop("data option is required")
  if(is.data.frame(dat)) dat <- data.matrix(dat)
  set.seed(seedno)
  models <- modelnames
  N <- nrow(dat) #length(membership) #Total number of observations
  p <- ncol(dat) #Dimensionality of x
  gr <- G
  results <- list()
  #################Functions used######################
  fnz <- function(git, psi, Sig, invSig, m_g, z, pi_g, nu, label) {
    #component label indicator
    z_g <- matrix(NA, nrow = N, ncol = git)
    fz <- z
    zold <- z
    xpart <- matrix(NA, N, git)
    kpart <- matrix(NA, N, git)
    pn <- matrix(NA, N, git)
    for (g in 1:git) {
      uplim_part <- t(psi[[g]]) %*% mpower(Sig[[g]], p = -1/2)
      uplim <- m_g[[g]] %*% t(uplim_part)
      pn[, g] <- drop(pnorm(uplim, log.p = TRUE))
      xpart[, g] <- try(-1/2 * log(det(Sig[[g]])) - 1/2 * (mahalanobis(as.matrix(dat), center = mu[[g]], cov = invSig[[g]], inverted = TRUE)) ^ betag[g], silent = TRUE)
      kpart[, g] <- try(log(p) + log(gamma(p/2)) - p/2 * log(pi) - log(gamma(1 + p/(2 * betag[g]))) -  (1 + p/(2 * betag[g])) * log(2), silent = TRUE)
      if (is.numeric(kpart[, g]) && is.numeric(xpart) && is.numeric(pn))
        z_g[, g] <- exp(log(pi_g[g]) + log(2) + kpart[, g] + xpart[, g] + pn[, g])
      if (!is.numeric(z_g[, g]))
        z_g[, g] <- Inf
    }
    if(any(!is.na(label))) z_g[!is.na(label),] <- unmap(label[!is.na(label)],git) #for classification
    # fz <- try(z_g/rowSums(z_g), silent = TRUE)
    fz <- try(z_g ^ nu / rowSums(z_g ^ nu), silent = TRUE) #nu for annealing...
    if (!is.numeric(fz))
      fz <- Inf
    return(list(z = fz, forz = z_g, zold = zold))
  }
  stopcheck <- function(it, incparll) {
    #stopping criterion
    if ((it > (maxit - 1)) || (incparll[it] == Inf) || (incparll[it] == incparll[it - 1])) {
      stop <- 1
    }
    if (incparll[it] != Inf && it < maxit) {
      if (incparll[it] != incparll[it - 1]) {
        ait1less <- (incparll[it] - incparll[(it - 1)])/(incparll[(it - 1)] - incparll[(it -
                                                                                          2)])
        llinf <- incparll[(it - 1)] + (1/(1 - ait1less)) * (incparll[it] - incparll[(it -
                                                                                       1)])
        if (it > 2) {
          stop <- 0 #in case the condition is not satisfied
          if (abs(llinf - incparll[(it - 1)]) < eps) {
            stop <- 1
          }
        }
      }
    }
    mon <- na.omit(incparll[1:it])
    if (!all(mon == cummax(mon))) {
      stop <- 1
    }
    return(stop)
  }
  ms_x <- function(it, z, git, pi_g, mu, invSig, m_g, Sig, psi, betag) {
    #Newton raphson for mu
    condition <- 1
    pi_gnew <- colMeans(as.matrix(z))
    munew <- list()
    for (g in 1:git) {
      mutrial = rep(NA, p)
      Sg <- mahalanobis(dat, center = mu[[g]], cov = invSig[[g]], inverted = TRUE)
      invsqsig <- mpower(Sig[[g]], p = -1/2)
      uplim_part <- t(psi[[g]]) %*% invsqsig
      uplim <- m_g[[g]] %*% t(uplim_part)
      dn <- drop(dnorm(uplim))
      pn <- drop(pnorm(uplim))
      dpratio <- dn/pn
      app_ind <- NULL
      app_ind <- which(uplim < -37) #indices for the approximation where dnorm/pnorm results in NaN or Inf...
      if (!is.null(app_ind)) dpratio[app_ind] <- unlist(lapply(-uplim[app_ind], FUN = function(x) 1 / (1/x - 1/x ^ 3 + 3/x ^ 5 - 3*5/x ^ 7 + 3*5*7/x ^ 9)))
      #Note that the negatives are taken as the division of a density and area cannot be negative.
      tupup <- t(uplim_part) %*% uplim_part
      non0 <- which(z[,g] != 0)
      fun_deriv <- function(mun) {
        betag[g] * colSums(z[non0, g] * Sg[non0] ^ (betag[g] - 1) * m_g[[g]][non0, ] %*% invSig[[g]]) - sum(z[non0, g] * dpratio[non0]) * t(uplim_part)
      }
      fun_hess <- function(mun) {
        Reduce('+', mapply(FUN = '*', lapply(X = non0, FUN = function(i) {invSig[[g]] %*% (dat[i, ] - mun) %*% t(invSig[[g]] %*% (dat[i, ] - mun))}), as.numeric(- 2 * z[non0, g] * betag[g] * (betag[g] - 1) * Sg[non0] ^ (betag[g] - 2)), SIMPLIFY = FALSE)) - sum(z[non0, g] * betag[g] * Sg[non0] ^ (betag[g] - 1)) * invSig[[g]] - sum(z[non0, g] * drop(uplim)[non0] * dpratio[non0]) * tupup - sum(z[non0, g] * (dpratio[non0]) ^ 2) * tupup
      }

      if (!(rcond(fun_hess(mu[[g]])) < sqrt(.Machine$double.eps))) {
        mutrial <- try(mu[[g]] - as.vector(solve(fun_hess(mu[[g]])) %*% fun_deriv(mu[[g]])), silent = TRUE)
      } else if (!(rcond(diag(diag(fun_hess(mu[[g]])))) < sqrt(.Machine$double.eps))) {
        mutrial <- try(mu[[g]] - as.vector(solve(diag(diag(fun_hess(mu[[g]])))) %*% fun_deriv(mu[[g]])), silent = TRUE)   #if less than tolerance, ignore non-diagonal terms.
      } else {
        mutrial <- try(mu[[g]] - as.vector(solve(diag(diag(fun_hess(mu[[g]])) + 10000)) %*% fun_deriv(mu[[g]])), silent = TRUE)
      }
      if (all(is.finite(mutrial))) {
        munew[[g]] <- mutrial
      }
      if (!all(is.finite(mutrial))) {
        munew[[g]] <- rep(NA, p)
        condition <- 0
      }
    }
    result <- list(pi_g = pi_gnew, mu = munew, pi_gold = pi_g, muold = mu, con = condition)
    return(result)
  }
  ms_betagV <- function(betag, git, z, m_g, invSig) {
    #Newton Raphson for beta - variable between groups
    betaold <- betag
    condition <- 1
    betanew = rep(NA, git)
    betan <- rep(NA, git)
    for (g in 1:git) {
      betan <- betag[g]
      Sg <- mahalanobis(m_g[[g]], center = rep(0, p), cov = invSig[[g]], inverted = TRUE)
      n_g <- colSums(z)[g]

      fun_be <- function(betan) { #derivative
        p/2 * n_g/betan^2 * digamma(1 + p/(2 * betan)) + p/2 * log(2)/betan^2 * n_g +
          sum(-z[, g]/2 * Sg^betan * as.numeric(log(Sg)))
      }
      fun_secbe <- function(betan) { #second derivative
        -p * n_g/betan^3 * digamma(1 + p/(2 * betan)) + p/2 * n_g/betan^2 * trigamma(1 +
                                                                                       p/(2 * betan)) * (-p/(2 * betan^2)) + (-p)/2 * log(2)/betan^3 * n_g + sum(-z[,
                                                                                                                                                                    g]/2 * Sg^betan * as.numeric(log(Sg))^2)
      }
      betanew[g] <- try(betan - fun_be(betan)/fun_secbe(betan), silent = TRUE)
      if (is.finite(betanew[g])) {
        if (betanew[g] < 0.05) {
          betanew[g] <- 0.05
        } else if (betanew[g] > 200) {
          betanew[g] <- 200
        }
      }
      if (!is.finite(betanew[g])) {
        betanew[g] <- NA
        condition <- 0
      }
    }
    result <- list(betag = betanew, betaold = betaold, con = condition)
    return(result)
  }
  ms_betagE <- function(betag, git, z, m_g, invSig) {
        #Newton Raphson for beta - equal between groups
    betanew <- NA
    betaold <- betag
    condition <- 1
    fn_be <- function(betan) {
      val <- numeric(git)
      valse <- numeric(git)
      for (g in 1:git) {
        Sg <- mahalanobis(m_g[[g]], center = rep(0, p), cov = invSig[[g]], inverted = TRUE)
        val[g] <- 0.5 * sum(z[, g] * Sg^betan * as.numeric(log(Sg)))
        valse[g] <- 0.5 * sum(z[, g] * Sg^betan * (as.numeric(log(Sg)))^2)
      }
      return(list(befi = sum(val), besec = sum(valse)))
    }
    fun_be <- function(betan) {
      N * p/(2 * betan^2) * (digamma(1 + p/(2 * betan)) + log(2)) - fn_be(betan)$befi
    }
    fun_secbe <- function(betan) {
      -p * N/betan^3 * digamma(1 + p/(2 * betan)) + p/2 * N/betan^2 * trigamma(1 +
                                                                                 p/(2 * betan)) * (-p/(2 * betan^2)) + (-p) * log(2)/betan^3 * N - fn_be(betan)$besec
    }
    ifelse(p > 1, betan <- betaold[1], betan <- betaold)
    betanew <- try(betan - fun_be(betan)/fun_secbe(betan), silent = TRUE)
    if (is.finite(betanew)) {
      if (betanew < 0.05) {
        betanew <- 0.05
      } else if (betanew > 200) {
        betanew <- 200
      }
    }
    if (!is.finite(betanew)) {
      betanew <- NA
      condition <- 0
    }
    result <- list(betag = rep(betanew, git), betaold = betaold, con = condition)
    return(result)
  }

  ms_psi <- function(dat, z, git, mu, m_g, Sig, psi) {
    #MM-update for skewness
    psinew <- NA
    condition <- 1
    for (g in 1:git) {
      mimi <- Reduce('+', lapply(X = 1:nrow(dat), FUN = function(i) {z[i, g] * m_g[[g]][i, ] %*% t(m_g[[g]][i, ])}))
      if (!(rcond(mimi) < sqrt(.Machine$double.eps))) {
        invmimi <- solve(mimi)
      invsqsig <- mpower(Sig[[g]], p = -1/2)
      uplim_part <- t(psi[[g]]) %*% invsqsig
      uplim <- m_g[[g]] %*% t(uplim_part)
      dn <- drop(dnorm(uplim))
      pn <- drop(pnorm(uplim))
      dpratio <- dn/pn
      app_ind <- NULL
      app_ind <- which(uplim < -37) #indices for the approximation from Irene's thesis. Needed for stability.
      if (!is.null(app_ind)) dpratio[app_ind] <- unlist(lapply(-uplim[app_ind], FUN = function(x) 1 / (1/x - 1/(x ^ 3) + 3/(x ^ 5) - 3*5/(x ^ 7) + 3*5*7/(x ^ 9))))
      val <- invmimi %*% (colSums(z[, g] * dpratio * m_g[[g]]) + drop(mimi %*% t(uplim_part)))
      # val <- invmimi %*% colSums(drop(uplim + dn/pn) * m_g)
      psi[[g]] <- drop(mpower(Sig[[g]], p = 1/2) %*% val)
      } else {
        condition <- 0
      }
    }
    return(list(psi = psi, con = condition))
  }
  testgrad.D = function(D = NULL, d = NULL, dat = NULL, Ak = NULL, b = NULL, g = NULL) {
    Dak = sweep(D, 2, 1/Ak, FUN = "*")
    invS = D %*% diag(1/Ak) %*% t(D)
    dat_x <- sweep(dat, 2, mu[[g]], "-") * (z[, g])^(1/(2 * betag[g]))
    q = mahalanobis(x = dat_x, center = rep(0, d), cov = invS, inverted = TRUE)
    qb_1 = as.numeric(q)^(b - 1)
    qb_1[which(q <= 0)] <- 0

    Rk = cov.wt(x = dat_x, wt = qb_1, center = rep(0, d), method = "ML")$cov * sum(qb_1)
    val = (2 * b * Rk) %*% Dak
    return(val)
  }
  testgrad.DE = function(D = NULL, d = NULL, dat = NULL, Ak = NULL, b = NULL) {
    Rk <- list()
    Dak <- list()
    vl <- list()
    for (g in 1:git) {
      Dak[[g]] = sweep(D[[g]], 2, 1/Ak[[g]], FUN = "*")
      invS = D[[g]] %*% diag(1/Ak[[g]]) %*% t(D[[g]])
      dat_x <- sweep(dat, 2, mu[[g]], "-") * (z[, g])^(1/(2 * b[g]))
      q = mahalanobis(x = dat_x, center = rep(0, d), cov = invS, inverted = TRUE)
      qb_1 = as.numeric(q)^(b[g] - 1)
      qb_1[which(q <= 0)] <- 0

      Rk[[g]] = cov.wt(x = dat_x, wt = qb_1, center = rep(0, d), method = "ML")$cov *
        sum(qb_1)
      vl[[g]] = (b[g] * Rk[[g]] %*% Dak[[g]])
    }
    val = 2 * Reduce("+", vl)
    return(val)
  }
  testvalb = function(dat = NULL, Ak = NULL, D = NULL, b = NULL, g = NULL) {
    d = ncol(dat)
    Dak = sweep(D, 2, 1/Ak, FUN = "*")
    invS = Dak %*% t(D)
    dat_x <- sweep(dat, 2, mu[[g]], "-") * (z[, g])^(1/(2 * betag[g]))
    q = mahalanobis(x = dat_x, center = rep(0, d), cov = invS, inverted = TRUE)
    val = sum(q^b)
    return(val)
  }
  testvalbE = function(dat = NULL, Ak = NULL, D = NULL, b = NULL) {
    d = ncol(dat)
    vl <- list()
    for (g in 1:git) {
      Dak = sweep(D[[1]], 2, 1/Ak[[g]], FUN = "*")
      invS = Dak %*% t(D[[1]])
      dat_x <- sweep(dat, 2, mu[[g]], "-") * (z[, g])^(1/(2 * b[g]))
      q = mahalanobis(x = dat_x, center = rep(0, d), cov = invS, inverted = TRUE)
      q[which(q < 0)] <- 0
      vl[[g]] = sum(q^b[g])
    }
    val = Reduce("+", vl)
    return(val)
  }
  projZ.X = function(Z = NULL, X = NULL) {
    XZ = t(X) %*% Z
    return(Z - X %*% (XZ + t(XZ))/2)
  }
  qr.Rplus = function(X) {
    z = qr(X)
    q = qr.Q(z, TRUE)
    r = qr.R(z, TRUE)
    q = q %*% diag(sign(diag(r)))
    r = diag(sign(diag(r))) %*% r
    return(list(q = q, r = r))
  }
  retractionZ.X = function(X, Z) {
    return(qr.Rplus(Z + X)$q)
  }
  newDb = function(D = NULL, Rdata = NULL, Ak = NULL, tmax = 100, b = NULL, g = NULL) {
    d = ncol(Rdata)
    eta = -1 * projZ.X(Z = testgrad.D(D = D, d = d, dat = Rdata, Ak = Ak, b = b, g = g),
                       X = D)
    tval0 = testvalb(dat = Rdata, Ak = Ak, D = D, b = b, g = g)
    tk = armijo(abs = c(1, 0.95, 1e-08), eta = eta, D = D, d = d, dat = Rdata, Ak = Ak,
                tval = tval0, mmax = 25, b = b, g = g)
    ifelse(tk == 0, Xk1 <- D, Xk1 <- retractionZ.X(Z = eta * tk, X = D))
    return(Xk1)
  }
  newDbE = function(D = NULL, Rdata = NULL, Ak = NULL, tmax = 100, b = NULL) {
    d = ncol(Rdata)
    eta = -1 * projZ.X(Z = testgrad.DE(D = D, d = d, dat = Rdata, Ak = Ak, b = b),
                       X = D[[1]])
    tval0 = testvalbE(dat = Rdata, Ak = Ak, D = D, b = b)
    tk = armijoE(abs = c(1, 0.95, 1e-08), eta = eta, D = D, d = d, dat = Rdata, Ak = Ak,
                 tval = tval0, mmax = 25, b = b)
    if (tk == 0)
      Xk1 = D[[1]]
    else Xk1 = retractionZ.X(Z = eta * tk, X = D[[1]])
    return(Xk1)
  }
  newaE <- function(s0 = NULL, dat = NULL, beta = NULL) {
    compsig <- list()
    condition <- 1
    for (g in 1:git) {
      w = dat[[g]] %*% diag(1/sqrt(s0[[g]]))
      ww = apply(w^2, 1, sum)^(beta[g] - 1)
      ww[which(!is.finite(ww))] <- 0
      if (any(ww > 0)) {
        d = ncol(dat[[g]])
        lam = cov.wt(w, wt = ww, center = rep(0, d), method = "ML")$cov * sum(ww) *
          beta[g]
        compsig[[g]] = diag(diag(s0[[g]]^(max(beta)/2)) %*% lam %*% diag(s0[[g]]^(max(beta)/2)))
      }
      Ap <- (1/N * Reduce("+", compsig))^(1/max(beta))
    }
    if (!any(ww > 0)) {
      Ap <- diag(p)
      condition <- 0
    }
    return(list(Ap = Ap, con = condition)) # matrix A raised to a power
  }
  newaEbless1 <- function(s0 = NULL, beta = NULL, dat = NULL) {
    compsig <- list()
    maha <- list()
    for (g in 1:git) {
      maha[[g]] <- mahalanobis(dat[[g]], center = rep(0, p), cov = diag(1/s0[[g]]),
                               inverted = TRUE)
    }
    for (g in 1:git) {
      compsig[[g]] <- diag(beta[g] * matrix(rowSums(sapply(1:N, function(i) {
        z[i, g] * (maha[[g]][i])^(beta[g] - 1) * (dat[[g]][i, ]) %*% t(dat[[g]][i,
                                                                                ])
      })), p, p))
    }
    A <- 1/N * Reduce("+", compsig)
    return(A)
  }
  newaV <- function(s0 = NULL, beta = NULL, dat = NULL, g = g, ng = ng) {
    maha <- mahalanobis(dat, center = rep(0, p), cov = diag(1/s0), inverted = TRUE)
    mahab <- maha^(beta - 1)
    mahab[which(maha == 0)] <- 0
    compA <- diag(beta/ng * matrix(rowSums(sapply(1:N, function(i) {
      z[i, g] * mahab[i] * (dat[i, ]) %*% t(dat[i, ])
    })), p, p))
    return(compA)
  }
  newa <- function(s0 = NULL, dat = NULL, beta = NULL, ng = NULL) {
    condition <- 1
    w = dat %*% diag(1/sqrt(s0))
    ww = apply(w^2, 1, sum)^(beta - 1)
    ww[which(!is.finite(ww))] <- 0
    if (any(ww > 0)) {
      d = ncol(dat)
      lam = cov.wt(w, wt = ww, center = rep(0, d), method = "ML")$cov * sum(ww) *
        beta/ng
      sig = diag(diag(s0^(beta/2)) %*% lam %*% diag(s0^(beta/2)))^(1/beta)
    }
    if (!any(ww > 0)) {
      sig <- diag(p)
      condition <- 0
    }
    return(list(sig = sig, con = condition))
  }
  armijo = function(abs = c(1, 0.1, 0.1), eta = NULL, D = NULL, d = NULL, dat = NULL,
                    Ak = NULL, tval = NULL, mmax = 100, b = NULL, g = NULL) {
    # abs is a numeric vector with alpha > 0, beta, sigma in (0,1)
    eta2 = sum(eta^2)
    m = 1
    obj = -1
    tkseq = seq(0, log(1e-08)/log(abs[2]), length.out = mmax)
    while (obj < 0 & m < mmax) {
      tk = exp(tkseq[m] * log(abs[2])) * abs[1]
      Reta = retractionZ.X(Z = eta * tk, X = D)
      Rval = testvalb(D = Reta, dat = dat, Ak = Ak, b = b, g = g) #- tval
      obj = tval - Rval #+ abs[3]* tk * eta2
      m = m + 1
    }
    if (obj < 0)
      tk = 0
    return(tk)
  }
  armijoE = function(abs = c(1, 0.1, 0.1), eta = NULL, D = NULL, d = NULL, dat = NULL,
                     Ak = NULL, tval = NULL, mmax = 100, b = NULL) {
    # abs is a numeric vector with alpha > 0, beta, sigma in (0,1)
    eta2 = sum(eta^2)
    Reta = list()
    m = 1
    obj = -1
    tkseq = seq(0, log(1e-08)/log(abs[2]), length.out = mmax)
    while (obj < 0 & m < mmax) {
      tk = exp(tkseq[m] * log(abs[2])) * abs[1]
      Reta[[1]] = retractionZ.X(Z = eta * tk, X = D[[1]])
      Rval = testvalbE(D = Reta, dat = dat, Ak = Ak, b = b) #- tval
      obj = tval - Rval #+ abs[3]* tk * eta2
      m = m + 1
    }
    if (obj < 0)
      tk = 0
    return(tk)
  }
  msEII <- function(it, z, Sig, git, betag, m_g) {
    condition <- 1
    Sigold <- Sig
    if (it == 1) {
      Sig[[1]] = diag(1)
    }
    Sigret <- list()
    mg <- array(NA, dim = c(N, p, git))
    for (g in 1:git) {
      mg[, , g] <- z[, g]^(0.5/betag[g]) * m_g[[g]]
    }
    fun_EII <- function(sigma) {
      p * N - sum(sapply(1:git, function(g) {
        betag[g] * sigma^(-betag[g]) * sum(rowSums(mg[, , g] * mg[, , g])^betag[g])
      }))
    }
    res <- try(uniroot(fun_EII, c(1e-04, (diag(Sig[[1]])[1] + 100)))$root, silent = TRUE)
    if (is.finite(res)) {
      for (g in 1:git) {
        Sigret[[g]] <- res * diag(p)
      }
    }
    if (!is.finite(res)) {
      Sigret[[g]] <- diag(NA, p)
      condition <- 0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition)
    return(result)
  }
  msVII <- function(it, z, Sig, git, betag, m_g) {
    condition <- 1

    Sigold <- Sig
    Signew <- rep(NA, git)
    Sigret <- list()
    for (g in 1:git) {
      Signew[g] <- try((betag[g] * sum(z[, g] * rowSums(m_g[[g]] * m_g[[g]])^betag[g])/(p *
                                                                                          colSums(z)[g]))^(1/betag[g]), silent = TRUE )
      if (is.finite(Signew[g]))
        Sigret[[g]] <- Signew[g] * diag(p)
    }
    if (!all(is.finite(Signew))) {
      for (g in 1:git) Sigret[[g]] <- diag(NA, p)
      condition <- 0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition)
    return(result)
  }
  msEEX <- function(it, z, Sig, git, betag, m_g, modelname) {
    Sigold <- Sig
    condition <- 1
    Sigret <- list()
    Anew <- list()
    Dnew <- list()
    Ddata <- list()
    DdataEbless1 <- list()
    if (it == 1) {
      for (g in 1:git) {
        Anew[[g]] <- rep(1, p)
        if (!substr(modelname, 3, 3) == "I")
          Dnew[[g]] <- diag(p)
      }
    }
    for (g in 1:git) {
      Anew[[g]] <- eigen(Sig[[g]])$values
      if (substr(modelname, 3, 3) == "E")
        Dnew[[g]] <- eigen(Sig[[1]])$vectors
      if (substr(modelname, 3, 3) == "V")
        Dnew[[g]] <- eigen(Sig[[g]])$vectors
    }
    if (substr(modelname, 3, 3) == "I") {
      for (g in 1:git) Ddata[[g]] <- (z[, g])^(1/(2 * betag[g])) * m_g[[g]]
      for (g in 1:git) DdataEbless1[[g]] <- m_g[[g]]
    }
    if (!substr(modelname, 3, 3) == "I") {
      for (g in 1:git) Ddata[[g]] <- (z[, g])^(1/(2 * betag[g])) * m_g[[g]] %*% Dnew[[g]]
      for (g in 1:git) DdataEbless1[[g]] <- m_g[[g]] %*% Dnew[[g]]
    }
    if (!substr(modelname, 3, 3) == "I") {
      for (g in 1:git) DdataEbless1[[g]] <- m_g[[g]] %*% Dnew[[g]]
    }
    if (all(betag < 1))
      Anew[[1]] <- try(newaEbless1(s0 = Anew, dat = DdataEbless1, beta = betag),
                       silent = TRUE)
    else if (any(betag >= 1)) {
      fnres <- newaE(s0 = Anew, dat = Ddata, beta = betag)
      Anew[[1]] <- fnres$Ap
      condition <- fnres$con
    }
    if (all(is.finite(Anew[[1]])) & all(Anew[[1]] > 0)) {
      if (git > 1) {
        for (g in 2:git) {
          Anew[[g]] <- Anew[[1]]
        }
      }
      if (all(Anew[[1]] > 0) && !substr(modelname, 3, 3) == "I") {
        if (substr(modelname, 3, 3) == "E") {
          Dnew[[1]] <- newDbE(D = Dnew, Ak = Anew, Rdata = dat, tmax = 100, b = betag)
          if (git > 1) {
            for (g in 2:git) {
              Dnew[[g]] <- Dnew[[1]]
            }
          }
        }
        if (substr(modelname, 3, 3) == "V") {
          Dnew[[g]] <- newDb(D = Dnew[[g]], Ak = Anew[[g]], Rdata = dat, tmax = 100,
                             b = betag[g], g = g)
        }
        for (g in 1:git) Sigret[[g]] <- Dnew[[g]] %*% diag(Anew[[g]]) %*% t(Dnew[[g]])
      }
      if (all(Anew[[1]] > 0) && substr(modelname, 3, 3) == "I") {
        for (g in 1:git) Sigret[[g]] <- diag(Anew[[g]])
      }
    } else if (!all(is.finite(Anew[[1]])) || !all(Anew[[1]] > 0)) {
      for (g in 1:git) Sigret[[g]] <- diag(NA, p)
      condition <- 0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition, D = Dnew, A = Anew)
    return(result)
  }
  msVVX <- function(it, z, Sig, git, betag, m_g, modelname) {
    Sigold <- Sig

    condition <- 1
    Sigret <- list()
    Anew <- list()
    Afin <- list()
    Dnew <- list()
    Dfin <- list()
    Ddata <- list()
    if (it == 1) {
      for (g in 1:git) {
        Anew[[g]] <- rep(1, p)
        if (substr(modelname, 3, 3) == "E")
          Dnew[[g]] <- diag(p)
      }
    }

    for (g in 1:git) {
      Anew[[g]] <- eigen(Sig[[g]])$values
      if (substr(modelname, 3, 3) == "E")
        Dnew[[g]] <- eigen(Sig[[1]])$vectors
    }
    if (substr(modelname, 3, 3) == "I") {
      for (g in 1:git) Ddata[[g]] <- (z[, g])^(1/(2 * betag[g])) * m_g[[g]]
    }
    if (substr(modelname, 3, 3) == "E") {
      for (g in 1:git) Ddata[[g]] <- (z[, g])^(1/(2 * betag[g])) * m_g[[g]] %*% Dnew[[g]]
    }
    for (g in 1:git) {
      if (betag[g] >= 1) {
        fnres <- newa(s0 = Anew[[g]], dat = Ddata[[g]], beta = betag[g], ng = sum(z[,
                                                                                    g]))
        if (fnres$con == 1)
          Anew[[g]] <- fnres$sig
        if (fnres$con == 0)
          Anew[[g]] <- diag(NA, p)
      }
      if (betag[g] < 1) {
        if (substr(modelname, 3, 3) == "E") {
          Anew[[g]] <- newaV(s0 = Anew[[g]], dat = m_g[[g]] %*% Dnew[[g]], beta = betag[g],
                             g = g, ng = sum(z[, g]))
        }
        if (substr(modelname, 3, 3) == "I") {
          Anew[[g]] <- newaV(s0 = Anew[[g]], dat = m_g[[g]], beta = betag[g], g = g,
                             ng = sum(z[, g]))
        }
      }
    }
    if (substr(modelname, 3, 3) == "E") {
      if (all(is.finite(unlist(Anew))) && all(unlist(Anew) > 0)) {
        Dnew[[1]] <- newDbE(D = Dnew, Ak = Anew, Rdata = dat, tmax = 100, b = betag)
        if (git > 1) {
          for (g in 2:git) {
            Dnew[[g]] <- Dnew[[1]]
          }
        }
        for (g in 1:git) Sigret[[g]] <- Dnew[[g]] %*% diag(Anew[[g]]) %*% t(Dnew[[g]])
      }
    }
    if (substr(modelname, 3, 3) == "I") {
      if (all(is.finite(unlist(Anew))) && all(unlist(Anew) > 0)) {
        for (g in 1:git) Sigret[[g]] <- diag(Anew[[g]])
      }
    }
    if (!all(is.finite(unlist(Anew))) || !all(unlist(Anew) > 0)) {
      for (g in 1:git) Sigret[[g]] <- diag(NA, p)
      condition <- 0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition, D = Dnew, A = Anew)
    return(result)
  }
  msVVV <- function(it, z, Sig, invSig, git, betag, m_g) {
    Sigold <- Sig

    condition <- 1
    Sigret <- list()
    Anew <- list()
    Dnew <- list()
    for (g in 1:git) {
      if (betag[g] < 1) {
        maha <- mahalanobis(m_g[[g]], center = rep(0, p), cov = invSig[[g]], inverted = TRUE)
        Sigret[[g]] <- 1/colSums(z)[g] * betag[g] * matrix(rowSums(sapply(1:N, function(i) {
          z[i, g] * maha[i]^(betag[g] - 1) * (m_g[[g]][i, ]) %*% t(m_g[[g]][i, ])
        })), p, p)
      }
      if (betag[g] >= 1) {
        if (it == 1) {
          Anew[[g]] <- rep(1, p)
          Dnew[[g]] <- diag(p)
        }
        Anew[[g]] <- eigen(Sig[[g]])$values
        Dnew[[g]] <- eigen(Sig[[g]])$vectors

        Ddata <- list()
        Ddata[[g]] <- (z[, g])^(1/(2 * betag[g])) * m_g[[g]] %*% Dnew[[g]]
        fnres <- newa(s0 = Anew[[g]], dat = Ddata[[g]], beta = betag[g], ng = sum(z[,
                                                                                    g]))
        if (fnres$con == 1)
          Anew[[g]] <- fnres$sig
        if (fnres$con == 0)
          Anew[[g]] <- diag(NA, p)
        if (all(is.finite(Anew[[g]])) && all(Anew[[g]] > 0)) {
          Dnew[[g]] <- newDb(D = Dnew[[g]], Ak = Anew[[g]], Rdata = dat, tmax = 100,
                             b = betag[g], g = g)
          Sigret[[g]] <- Dnew[[g]] %*% diag(Anew[[g]]) %*% t(Dnew[[g]])
        }
      }
    }
    if (!all(is.finite(unlist(Anew))) || !all(unlist(Anew) > 0)) {
      for (g in 1:git) Sigret[[g]] <- diag(NA, p)
      condition <- 0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition)
    return(result)
  }

  ncovpar <- function(modelname = NULL, p = NULL, G = NULL) {
    #Number of parameters fit in the specific model
    if (is.null(p))
      stop("p is null")
    if (is.null(G))
      stop("G is null")
    if (is.null(modelname))
      stop("modelname is null")
    if (substr(modelname, 1, 3) == "EII")
      npar = 1
    else if (substr(modelname, 1, 3) == "VII")
      npar = G
    else if (substr(modelname, 1, 3) == "EEI")
      npar = p
    else if (substr(modelname, 1, 3) == "VVI")
      npar = p * G
    else if (substr(modelname, 1, 3) == "EEE")
      npar = p * (p + 1)/2
    else if (substr(modelname, 1, 3) == "EEV")
      npar = G * p * (p + 1)/2 - (G - 1) * p
    else if (substr(modelname, 1, 3) == "VVE")
      npar = p * (p + 1)/2 + (G - 1) * p
    else if (substr(modelname, 1, 3) == "VVV")
      npar = G * p * (p + 1)/2
    else stop("modelname is not correctly defined")
    return(npar)
  }
  model.type <- function(modelname = NULL, it, z, Sig, invSig, git, betag, m_g) {
    #M-update dependent on model type
    if (is.null(modelname))
      stop("modelname is null")
    if (substr(modelname, 1, 3) == "EII")
      val = msEII(it, z, Sig, git, betag, m_g)
    else if (substr(modelname, 1, 3) == "VII")
      val = msVII(it, z, Sig, git, betag, m_g)
    else if (substr(modelname, 1, 3) == "EEI")
      val = msEEX(it = it, z = z, Sig = Sig, git = git, betag = betag, m_g = m_g, modelname = modelname)
    else if (substr(modelname, 1, 3) == "VVI")
      val = msVVX(it = it, z = z, Sig = Sig, git = git, betag = betag, m_g = m_g, modelname = modelname)
    else if (substr(modelname, 1, 3) == "EEE")
      val = msEEX(it = it, z = z, Sig = Sig, git = git, betag = betag, m_g = m_g, modelname = modelname)
    else if (substr(modelname, 1, 3) == "EEV")
      val = msEEX(it = it, z = z, Sig = Sig, git = git, betag = betag, m_g = m_g, modelname = modelname)
    else if (substr(modelname, 1, 3) == "VVE")
      val = msVVX(it = it, z = z, Sig = Sig, git = git, betag = betag, m_g = m_g, modelname = modelname)
    else if (substr(modelname, 1, 3) == "VVV")
      val = msVVV(it = it, z = z, Sig = Sig, invSig = invSig, git = git, betag = betag, m_g = m_g)
    else stop("modelname or covtype is not correctly defined")
    if (!is.list(val))
      val = list(sigma = val)
    return(val)
  }
  npar.model <- function(modelname = NULL, p = NULL, G = NULL) {
    #Total number of parameters fit
    val = NULL
    val[1] = G - 1 #for pi
    val[2] = G * p #for mean for data
    val[3] = ncovpar(modelname = modelname, p = p, G = G) #Covariance for data
    if (substr(modelname, 4, 4) == "V")
      val[4] = G #for unconstrained beta_g
    if (substr(modelname, 4, 4) == "E")
      val[4] = 1 #for constrained beta_g
    val[5] = G * p
    val = sum(val)
    return(val)
  }

  ###############################
  modelbic <- matrix(NA, nrow = length(gr) * length(models), ncol = 5)
  colnames(modelbic) <- c("G", "BIC", "df", "Loglik", "ICL")
  rownames(modelbic) <- rep(models, length(gr))
  mbi <- 1
  git <- gr[1]
  gitupper <- gr[length(gr)]
  zlist <- list()
  while (git < (gitupper + 1)) {
    if (start == "kmeans") {
      inz <- kmeans(dat, git, nstart = kmeansinit)
      initialz <- data.matrix(unmap(inz$cluster, git)) #Kmeans initialization
    } else if (start == "random") {
      initialz <- data.matrix(unmap(sample(1:git, N, replace = TRUE), git))
    }
    if(any(!is.na(label))) initialz[!is.na(label),] <- unmap(label[!is.na(label)],git) # for classification when some labels are known
    zlist[[git]] <- map(initialz)
    mit <- 1
    while (mit < length(models) + 1) {
      modelname <- models[mit]
      z <- initialz
      zold <- z
      #################Parameter initialization##################
      pi_g <- colMeans(initialz)
      mu <- list()
      Sig <- list()
      invSig <- list()
      m_g <- list()
      for (g in 1:git) {
        tmp <- cov.wt(as.matrix(dat), wt = (initialz[, g]), center = TRUE,  method = "ML")
        mu[[g]] <- tmp$center
        Sig[[g]] <- tmp$cov
        if (!(rcond(Sig[[g]]) < sqrt(.Machine$double.eps))) {
          invSig[[g]] <- solve(Sig[[g]])
        } else {
          Sig[[g]] <- invSig[[g]] <- diag(p)
        }
        m_g[[g]] <- sweep(dat, 2, mu[[g]], "-") #
      }
      betag <- rep(1, git) #initialize beta
      betaold <- NULL
      psi <- rep(list(rep(0, p)), git) #initialize skewness

      if (psistart=="est"){
        for (g in 1:git){
          psi[[g]]<-(colMeans(dat*z[,g])-apply(dat*z[,g],2,median))/sqrt(diag(Sig[[g]]))
        }
      }

      if (!is.null(anneal)) {
        ff <- fnz(git = git, psi = psi, Sig = Sig, invSig = invSig, m_g = m_g, z = z, pi_g = pi_g, nu = anneal[1], label=label)
      } else {
        ff <- fnz(git = git, psi = psi, Sig = Sig, invSig = invSig, m_g = m_g, z = z, pi_g = pi_g, nu = 1,label=label)
      }
      z <- ff$z
      Sigold <- list()
      incparll <- NULL
      posdef <- 1
      #######################GEM loop starts#####################
      it <- 1
      stop <- 0
      deg <- 1
      while (stop < 1) {
        if(it %% 1 == 0) {if (substr(models[mit], 4, 4) == "E") {
          bg <- ms_betagE(betag = betag, git = git, z = z, m_g = m_g, invSig = invSig)
          if (bg$con != 0)
            betag <- bg$betag
        } else {
          bg <- ms_betagV(betag = betag, git = git, z = z, m_g = m_g, invSig = invSig)
          if (bg$con != 0)
            betag <- bg$betag
        }
        if (any(!is.finite(betag)) || bg$con == 0) {
          deg <- 0
          stop <- 1
          incparll[it] <- Inf
          break
        }
        mst_x <- ms_x(it = it, z = z, git = git, pi_g = pi_g, mu = mu, invSig = invSig, m_g = m_g, Sig = Sig, psi = psi, betag = betag)
        pi_g <- mst_x$pi_g
        mu <- mst_x$mu
        if (mst_x$con == 0) {
          deg <- 0
          stop <- 1
          incparll[it] <- Inf
          break
        }
        }
        for (g in 1:git) m_g[[g]] <- sweep(dat, 2, mu[[g]], "-") #mean centered data
        psi_l <- ms_psi(dat, z, git, mu, m_g, Sig, psi)
        if (psi_l$con == 1) {
          psi <- psi_l$psi
        } else {
          deg <- 0
          stop <- 1
          incparll[it] <- Inf
          break
        }
        if(it %% 1 == 0) {mstep_Sig <- try(model.type(modelname = modelname, it = it, z = z, Sig = Sig, invSig = invSig, git = git, betag = betag, m_g = m_g), silent = TRUE)
        if (mstep_Sig$con != 0) Sig <- mstep_Sig$Sig
        for (g in 1:git) {
          if (any(eigen(Sig[[g]])$values < 0) || mstep_Sig$con == 0) {
            posdef <- 0
            incparll[it] <- Inf
            stop <- 1
            break
          }
        }
        } #If I update Sig every 2 iterations only, this speeds up the code a lot and still seems to work well in terms of the overall density fit; mixed results with the psi and Sig estimates.
        for (g in 1:git) {
          if (!(rcond(Sig[[g]]) < sqrt(.Machine$double.eps)))
            invSig[[g]] <- try(solve(Sig[[g]]), silent = TRUE)
          if (!is.numeric(invSig[[g]]) || rcond(Sig[[g]]) < sqrt(.Machine$double.eps)) {
            stop <- 1
            posdef <- 0
            incparll[it] <- Inf
            break
          }
        }
        #inverse of Sig---to be used in next iteration
        if (!is.null(anneal) & git > 1) {
          if (it < length(anneal)) {
            ff <- fnz(git = git, psi = psi, Sig = Sig, invSig = invSig, m_g = m_g, z = z, pi_g = pi_g, nu = anneal[it + 1], label = label)
          } else
            ff <- fnz(git = git, psi = psi, Sig = Sig, invSig = invSig, m_g = m_g, z = z, pi_g = pi_g, nu = 1, label=label)
        } else ff <- fnz(git = git, psi = psi, Sig = Sig, invSig = invSig, m_g = m_g, z = z, pi_g = pi_g, nu = 1, label = label)
        z <- ff$z

        if (is.null(anneal) | git == 1) {
          incparll[it] <- sum(log(rowSums(ff$forz)))
          if (!is.finite(incparll[it])) {
            stop <- 1
            posdef <- 0
            break
          }
        } else if (!is.null(anneal)) {
          if (it > (length(anneal))) {
            incparll[it] <- sum(log(rowSums(ff$forz)))
            if (!is.finite(na.omit(incparll[it]))) {
              stop <- 1
              posdef <- 0
              incparll[it] <- Inf
              break
            }
          }
        }
        #here posdef doesn't refer to positive definiteness not being met
        if (it > 2 && is.finite(incparll[it]) && is.null(anneal)) {
          stop <- stopcheck(it, incparll)
        } else if (!is.null(anneal)) {
          if (it > (length(anneal) + 2)) {
            stop <- stopcheck(it, incparll)
          }
        }
        if (git == 1 & it > 1) {
          if (abs(incparll[it] - incparll[(it - 1)]) < eps) stop <- 1
        }
        if (any(!is.finite(z))) {
          stop <- 1
          incparll[it] <- Inf
          break
        }
        if (it > maxit || posdef == 0) {
          stop <- 1
          break
        }
        it <- it + 1
        if (it %% 1000 == 0) {
          cat(unlist(psi),"\t")
          cat(it,"\t")
        }
      }#for EM - for it
      mon <- na.omit(incparll[1:it])
      if (deg == 1 && it > 1 && all(mon == cummax(mon)) && all(is.finite(betag)) && all(is.finite(z)) && posdef == 1 && mst_x$con == 1) {
        para <- npar.model(modelname = modelname, p = p, G = git)
        if (incparll[(it - 1)] == Inf || ncol(z) != git)
          bic <- -Inf
        else bic <- 2 * incparll[(it - 1)] - log(N) * para
        #Here, higher BIC is better
        forICL <- function(g) {
          sum(log(z[max.col(z) == g, g]))
        }
        icl <- bic + 2 * sum(sapply(unique(max.col(z)), forICL))
        results[[paste("G", git, modelname, sep = "_")]] <- list(bic = bic, icl = icl,
                                                                 z = z, para = para, ll = incparll, pi_g = pi_g, mu = mu, Sig = Sig, psi = psi,
                                                                 beta = betag)
      } else {
        results[[paste("G", git, modelname, sep = "_")]] <- list(bic = NA, icl = NA,
                                                                 z = z, para = NA, ll = incparll[2:(it - 1)], pi_g = NA, mu = NA, Sig = NA, psi = NA,
                                                                 beta = NA, var = NA)
      }
      if (deg == 1 && it > 1 && all(mon == cummax(mon)) && all(is.finite(betag)) && all(is.finite(z)) && posdef == 1 && mst_x$con == 1) {
        if (all(as.numeric(betag) < 200)) {
          modelbic[mbi, 1] <- git
          modelbic[mbi, 2] <- bic
          modelbic[mbi, 3] <- para
          modelbic[mbi, 4] <- incparll[(it - 1)]
          modelbic[mbi, 5] <- icl
        }
      } else {
        modelbic[mbi, 1] <- git
        modelbic[mbi, 4] <- NA
      }
      if (verbose == TRUE)
        cat(modelname, "done for G = ", git, "\n")
      mbi <- mbi + 1
      mit <- mit + 1
    }
    git <- git + 1
  } #for git
  bestbic <- which.max(modelbic[, 2])
  besticl <- which.max(modelbic[, 5])
  #msc stands for model selection criterion
  timetaken <- proc.time() - ptm
  val <- list(call = match.call(), time = timetaken, modelnames = models, msc = modelbic,
              bicclassification = map(results[[bestbic]]$z), iclclassification = map(results[[besticl]]$z),
              bicselection = results[[bestbic]], iclselection = results[[besticl]], zlist = zlist)
  class(val) <- "spemix"
  return(invisible(val))
}
