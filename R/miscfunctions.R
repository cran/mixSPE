map <- function(x) {
  return(apply(X = x, MARGIN = 1, FUN = which.max))
}
unmap <- function(mapz, git) {
  zunmap <- matrix(data = 0, nrow = length(mapz), ncol = git)
  alphabet <- 1:git
  for (u in 1:git) {
    zunmap[which(mapz == alphabet[u]), u] <- 1
  }
  return(zunmap)
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


# contourplot <- function(col=TRUE, datm = NULL, kk = NULL, a = seq(-6, 5, length.out = 100), b = seq(-6, 5, length.out = 100)){
#   #Rough code to create a contour plot for two-component but two-dimensional data; this can be easily extended.
#   #datm=data in matrix form, kk is the object that holds results from the skew power exponential (SPE)
#   #fit, a and b are the x and y axes.
#   z_fun <- function(shape, del, a, b, p) {
#     #calculate the density for a particular component.
#     cond_z <- marg_z <- z1 <- matrix(0, nrow = length(a), ncol = length(b))
#     colnames(z1) <- colnames(marg_z) <- colnames(cond_z) <- b
#     rownames(z1) <- rownames(marg_z) <- rownames(cond_z) <- a
#     for (i in 1:length(a)) {
#       for (j in 1:length(b)) {
#         y <- c(a[i], b[j]) #current coordinates
#         tmp <- dspe(y, mu, Sig, del, bval = shape, p)
#         z1[i, j] <- pi_g * tmp$dens #Put density here...
#         marg_z[i, j] <- tmp$marg
#         cond_z[i, j] <- tmp$cond
#       }
#     }
#     return(list(z1 = z1, margz = marg_z, cond_z = cond_z))
#   }
#   dat <- datm
#   mu = kk$bicselection$mu[[1]]
#   Sig = kk$bicselection$Sig[[1]]
#   pi_g = kk$bicselection$pi_g[1]
#   zvals1 <- z_fun(shape = kk$bicselection$beta[1], del = kk$bicselection$psi[[1]], a, b, p=2)
#   forclist <- list()
#   forclist[[1]] <- zvals1$z1 #component 1
#   mu = kk$bicselection$mu[[2]]
#   Sig = kk$bicselection$Sig[[2]]
#   pi_g = kk$bicselection$pi_g[2]
#   zvals2 <- z_fun(shape = kk$bicselection$beta[2], del = kk$bicselection$psi[[2]], a, b, p=2)
#   forclist[[2]] <- zvals2$z1 #component 2
#
#   forcall <- Reduce("+", forclist) #mixture
#
#   pch_v <- kk$bicclassification
#   pch_v[pch_v == 1] <- 19
#   pch_v[pch_v == 2] <- 17
#
#   if(col==TRUE) {filled.contour(a, b, forcall, levels = c(1e-4, 0.001, 0.01, 0.05, 0.1, 0.3), xlab=expression(x[1]),ylab=expression(x[2]), plot.axes = { points(dat, col = kk$bicclassification, pch = pch_v, cex = 1.5); axis(1,cex.axis=1.25 ); axis(2,cex.axis=1.25) }, plot.title={ title(main="Estimated",cex.main=1.25) }, color.palette = terrain.colors)} else{
#     plot(dat,col=kk$bicclassification, pch = pch_v, cex = 1.5)
#     contour(a, b, forcall, levels = c(1e-4, 0.001, 0.01, 0.05, 0.1, 0.3), xlab=expression(x[1]),ylab=expression(x[2]), add=TRUE)
#   }
# }
