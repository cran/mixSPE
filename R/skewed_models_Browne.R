malG <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL ) {
  r2 = crossprod( gam,  t(x) - mu )^2*1/lam
  colSums( r2 )
} #for non-XXI models: (x-mu)' gamma lam^-1 gamma' (x-mu)

malI <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL ) {
  r2 = ( t(x) - mu )^2*1/lam
  colSums( r2 )
} #for XXI models: (x-mu)' lam^-1 (x-mu)

logden.to.loglik <- function(logz=NULL, pi=NULL) {
  val = sum(log(  colSums( exp(logz + log(pi)) ) ))
  return(val)
}

logden.to.weights <- function(logden=NULL, pi=NULL) {
  G = length(pi)
  if (G > 1) {
    logdpi = logden + log(pi)
    maxz   = apply(logdpi, 2, max)
    logdpi = t(logdpi) - maxz
    den    = exp(logdpi)
    w      = den/rowSums(den)
  } else w = matrix(1,nrow=ncol(logden), ncol=G)
  return(w)
}

MAP <- function(w) {
  z = apply(w, 1, function(z) { which(z==max(z)) } )
  z = as.numeric(z)
  return( z)
}




EMn <- function(data=NULL, G=2, model=NULL, n=10, label =NULL, gpar0=NULL ) {
  if (is.null(gpar0)) gpar  = rgpar(data=data, G=G, wt=NULL, model=model)
  else gpar = gpar0

  logden.lab = create.logden.lab(n=nrow(data), G=G, label=label)

  getlogden   = create.logz(x=data, G=G, model=model)
  gpar.update = create.update(x=data, G=G, model=model)

  loglik = numeric(n)
  for (i in 1:n) {
    logden    = getlogden(gpar) + logden.lab
    loglik[i] = logden.to.loglik(logden, gpar$pi)
    wt        = logden.to.weights(logden, gpar$pi)
    gpar      = gpar.update(gpar, wt)

  }

  logden = getlogden(gpar) + logden.lab
  wt     = logden.to.weights(logden, gpar$pi)
  mapz   = MAP( wt  )
  numpar = npar.model(model, ncol(data), G)

  val = list(loglik= loglik, gpar=gpar, z=wt, map=mapz, label=label, numpar = numpar, maxLoglik = loglik[n])
  return(val)
}





newBetaE <- function(beta=NULL, p=NULL, wt=NULL, ng=NULL, del=NULL) {

  G    = length(beta)
  beta = mean(beta)
  beta0 = beta;

  #log( x0^(x0^b0) )
  #db.logdel  = log( del^( t(t(del)^beta) ) )   #logdel*t(t(del)^beta)
  #db.logdel2 = ( log( del^( t(t(del)^(beta/2) ) ) )  )^2
  db.logdel  =  log( del )*t(t(del)^beta)
  db.logdel2 =  log( del )^2*t(t(del)^beta)

  if (any(del == 0)) {
    xdel0 = which(del == 0)
    db.logdel[xdel0] = 0
    db.logdel2[xdel0] = 0
  }


  p2b  = 1+p/(2*beta)
  #f1   =  p/(2*beta^2) * ( digamma(p2b) + log(2) ) - sum( wt/2*logdel*t(t(del)^beta) )/sum(ng)
  f1   =  p/(2*beta^2) * ( digamma(p2b) + log(2) ) - sum( wt/2*db.logdel) /sum(ng)

  f2a  = -p/(beta^3) * digamma(p2b)  - p^2/(4*beta^4)*trigamma(p2b)
  #f2b  = - p*log(2)/beta^3 - sum( wt/2*(logdel)^2*t(t(del)^beta) )/sum(ng)
  f2b  = - p*log(2)/beta^3 - sum( wt/2*db.logdel2 ) /sum(ng)
  beta = beta - f1/(f2a + f2b)
  beta[ beta >= 20] = 20
  beta[ beta < 0] = beta0

  if (is.nan(beta)) {
    stop("here beta")
  }

  beta = rep(beta, G)
  return(beta)
}

newBetaV <- function(beta=NULL, p=NULL, wt=NULL, ng=NULL, del=NULL) {
  db.logdel  =  log( del )*t(t(del)^beta)
  db.logdel2 =  log( del )^2*t(t(del)^beta)

  if (any(del == 0)) {
    xdel0 = which(del == 0)
    db.logdel[xdel0] = 0
    db.logdel2[xdel0] = 0
  }

  p2b = 1+p/(2*beta)
  f1  =  p/(2*beta^2) * ( digamma(p2b) + log(2) ) - colSums( wt/2*db.logdel )/ng

  f2a  = -p/(beta^3) * digamma(p2b)  - p^2/(4*beta^4)*trigamma(p2b)
  f2b  = - p*log(2)/beta^3 - colSums( wt/2*db.logdel2 )/ng
  beta = beta - f1/(f2a + f2b)
  beta[ beta > 20] = 20

  return(beta)
}

newBetaD <- function(beta=NULL, p=NULL, wt=NULL, ng=NULL, del=NULL) {
  ### Do not estimate this quantity.
  beta
}







create.logden.lab <- function(n=NULL, G=NULL, label=NULL)	{
  mat = matrix(0, nrow=G, ncol=n)

  # known is a numeric with
  # 0 if unknown group membership
  # 1,2,3,.. for label of known group
  if (!is.null(label)) {
    neg.Inf0 = c(-Inf,0)
    kw     = label !=0
    for (k in 1:G) {
      labk    = numeric(G)
      labk[k] = 1
      sumk    = sum(label == k)
      mat[ , label == k] = neg.Inf0[labk +1]
    }
  }
  return(mat)
}





getGamLamVII <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  p  = ncol(mu);
  G  = length(beta)

  r2 = sapply(1:G, function(k) { colSums( (t(x) - mu[k,])^2  )  })
  lam = (  (beta/(p*ng) )*rowSums( t(wt)* t(r2)^(beta) ) )^( 1/beta )

  list(lam=matrix(lam, nrow=G, ncol=p), gam=NULL)
}










getGamLamEII <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  n    = nrow(x)
  G    = length(beta)
  r2   = sapply(1:G, function(k) { colSums( (t(x) - mu[k,])^2  )  })
  coef = rowMeans(t(wt)* t(r2)^(beta) )
  p   =  ncol(mu)
  lam0 = lam[1,1]

  maxBeta = max(beta)
  lam = ( sum( beta*coef*lam0^(maxBeta - beta) )/p )^(1/maxBeta)

  if ( any(is.nan(lam)) ) {
    stop("here - lam ")
  }

  list(lam=matrix(lam, nrow=G, ncol=p), gam=NULL)
}

getGamLamVVI <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  p  = ncol(mu);
  G  = length(beta)

  lam = sapply(1:G, function(k) {
    tr2 = (t(x) - mu[k,])^2
    if (beta[k] > 1) lam = lamBG1( tr2, lam[k,], beta[k], ng[k], wt=wt[,k])
    else lam = lamBL1(tr2, lam[k,], beta[k], ng[k], wt[,k])
    lam
  })

  list(lam=t(lam), gam=NULL)
}


lamBL1 <- function(tr2=NULL, lam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  del = colSums( tr2/lam )
  delb = del^(beta-1)
  delb[del == 0] = 0
  lam = beta/ng* colSums( (wt*delb* t(tr2  )  ) )
  lam
}

lamBG1 <- function(tr2, lam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  del = colSums( tr2/lam )
  delb = del^(beta-1)
  delb[del == 0] = 0
  eta = beta/ng* colSums( (wt*delb* t( tr2 * lam^(beta-1)  )  ) )
  lam = eta^(1/beta)
  lam
}




getGamLamEEI <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {

  p = ncol(mu);
  G = length(beta)

  if (any(beta >1) ) {
    maxBeta = max(beta)
    eta = sapply(1:G, function(k) {
      tr2 = (t(x) - mu[k,])^2
      del = colSums( tr2/lam[k,] )
      delb = del^(beta[k]-1)
      delb[del == 0] = 0
      #beta[k]*colMeans( (wt[,k]*del^(maxBeta-1)* t( tr2 * lam[k,]^(maxBeta-1)  )  ) )
      beta[k]*colMeans( (wt[,k]*delb* t( tr2 * lam[k,]^(maxBeta-1)  )  ) )
    })
    lam = rowSums(eta)^(1/maxBeta)
  } else {
    lam = sapply(1:G, function(k) {
      tr2 = (t(x) - mu[k,])^2
      del = colSums( tr2/lam[k,] )

      delb = del^(beta[k]-1)
      delb[del == 0] = 0

      colMeans( wt[,k]*delb*t(tr2) )*beta[k]
    })
    lam = rowSums(lam)
  }

  if ( any(is.nan(lam)) ) {
    stop("here - lam ")
  }

  list(lam=matrix(lam, G, p, byrow=TRUE), gam=NULL)
}




getGamLamVVV <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  G = nrow(mu)

  for (k in 1:G) {
    tr  = t(x) - mu[k,]
    mal = colSums( crossprod( gam[,,k],  tr )^2*1/lam[k,] ) ## mahalanobis distance

    malb1 = mal^(beta[k] - 1)
    malb1[mal == 0] = 0

    Gmat   = tr %*% ( t(tr) * (beta[k]*wt[,k]/ng[k]*malb1 ) )

    if (beta[k] > 1) {
      gam[,,k] = gamBG1(lam=lam[k,], gam=gam[,,k], beta=beta[k], Gmat=Gmat, tr=tr, wt=wt[,k], ng=ng[k])

      tr2 = crossprod( gam[,,k],  tr )^2
      del = colSums( tr2/lam[k,] )

      delb1 = del^(beta[k] - 1)
      delb1[del == 0] = 0

      eta = beta[k]/ng[k]* colSums( (wt[,k]*delb1* t( tr2 * lam[k,]^(beta[k]-1)  )  ) )
      lam[k,] = eta^(1/beta[k])

    } else {
      temp     =  eigen(Gmat)
      lam[k,]  = temp$values
      gam[,,k] = temp$vectors
    }

  }

  list(lam=lam, gam=gam)
}




gamBG1 <- function(lam=NULL, gam=NULL, beta=NULL, Gmat=NULL, tr=NULL,  wt=NULL, ng=NULL) {
  ## gam and lam update when beta greater than 1
  Grad  = 1/2*(Gmat) %*% gam %*% diag(1/lam)
  eta   = -1*projZ.X(Z=Grad, X=gam)
  tk    = armijo(abs=c(1, .95, 1e-8 ), eta=eta, gam=gam, tr=tr, lam=lam, beta=beta,wt=wt, mmax=25)
  gam = retractionZ.X(Z=eta*tk, X=gam)

  return(gam)
}


projZ.X = function(Z=NULL,X=NULL) {
  XZ = t(X) %*% Z
  return( Z - X %*% ( XZ + t(XZ) )/2 )
}

objGam = function(gam=NULL, tr = NULL, lam=NULL,  beta=NULL, wt=NULL) {
  mal = colSums( crossprod( gam,  tr )^2*1/lam )
  val = 1/2*mean(wt*mal^beta)

  return(val)
}

retractionZ.X = function(X, Z) { return(qr.Rplus(Z+X)$q) }

qr.Rplus = function(X) {
  z = qr(X)
  q = qr.Q(z,TRUE)
  r = qr.R(z,TRUE)
  q = q %*% diag(sign(diag(r)))
  r = diag(sign(diag(r))) %*% r
  return(list(q=q,r=r))
}

armijo = function(abs=c(1, .1, .1), eta=NULL, gam=NULL, tr=NULL, lam=NULL, beta=NULL, wt=NULL, mmax=100) {
  # abs is a numeric vector with alpha > 0, beta, sigma in (0,1)

  tval = objGam(gam=gam, tr=tr, lam=lam, beta=beta, wt=wt)

  eta2  = sum(eta^2)
  m   = 1
  obj = -1
  tkseq = seq(0, log(1e-8)/log(abs[2]), length.out=mmax)
  while (obj < 0 & m < mmax  ) {
    tk      = exp( tkseq[m]*log(abs[2])) * abs[1]
    gam.eta = retractionZ.X(Z=eta*tk, X=gam)
    Rval    = objGam(gam=gam.eta, tr = tr, lam=lam,  beta=beta, wt=wt) #- tval
    obj     = tval - Rval #+ abs[3]* tk * eta2
    m = m+1
  }
  if (obj < 0 ) tk = 0
  return(tk)
}

getGamLamEEE <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  G = length(beta)
  p = ncol(x)
  N = nrow(x)

  Gmat = sapply(1:G, function(k) {
    tr  = t(x) - mu[k,]
    mal = colSums( crossprod( gam[,,k],  tr )^2*1/lam[k,] )

    malb = mal^(beta[k] - 1)
    malb[mal == 0] = 0

    beta[k]*(tr %*% ( t(tr) * (wt[,k]*malb ) )  )
  } )
  Gmat = rowSums(Gmat)/N
  dim(Gmat) = c(p,p)

  if (any(beta > 1) ) {

    Grad = Gmat %*% gam[,,1] %*% diag(1/lam[1,])
    gam  = gamEBG1(lam=lam, gam=gam[,,1], beta=beta, Grad=Grad, x=x, mu=mu, wt=wt, ng=ng)
    gam=array(gam, c(p,p,G) )

    #lam = lamEBG1(x=x, mu=mu, lam=lam, gam=gam, beta=beta, ng=ng, wt=wt)
    # lam = sapply(1:G, function(k) {
    #   tr2 = crossprod( gam[,,k], t(x) - mu[k,] )^2
    #   del = colSums( tr2/lam[k,] )
    #   colMeans( wt[,k]*del^(beta[k]-1)*t(tr2) )*beta[k]
    # })
    # lam = rowSums(lam)
    maxBeta = max(beta)
    eta = sapply(1:G, function(k) {
      tr2 = crossprod( gam[,,k], t(x) - mu[k,] )^2
      del = colSums( tr2/lam[k,] )
      beta[k]*colMeans( (wt[,k]*del^(beta[k]-1)* t( tr2 * lam[k,]^(maxBeta-1)  )  ) )
    })
    lam = rowSums(eta)^(1/maxBeta)

    lam = matrix(lam, G, p, byrow = TRUE)

  } else {
    temp = eigen( Gmat )
    lam  = matrix(temp$values, G, p, byrow = TRUE)
    gam  = array(temp$vectors, c(p,p,G) )
  }

  list(lam=lam, gam=gam )
}



gamEBG1 <- function(lam=NULL, gam=NULL, beta=NULL, Grad=NULL, x=NULL, mu=NULL, wt=NULL, ng=NULL) {
  gam0 = gam


  eta  = -1*projZ.X(Z=Grad, X=gam)
  tk   = armijoE(abs=c(1, .95, 0.01 ), eta=eta, gam=gam, x=x, mu=mu, lam=lam, beta=beta, wt=wt, mmax=25)
  gam = retractionZ.X(Z=eta*tk, X=gam0)

  return(gam)
}


armijoE <- function(abs=c(1, .1, .1), eta=NULL, gam=NULL, x=NULL, mu=NULL, lam=NULL, beta=NULL, wt=NULL, mmax=100, obj) {
  # common GAM
  # abs is a numeric vector with alpha > 0, beta, sigma in (0,1)
  tval = objGamE(gam=gam, x=x, mu=mu, lam=lam, beta=beta, wt=wt)

  eta2  = sum(eta^2)
  m   = 1
  obj = -1
  tkseq = seq(0, log(1e-8)/log(abs[2]), length.out=mmax)
  while (obj < 0 & m < mmax  ) {
    tk      = exp( tkseq[m]*log(abs[2])) * abs[1]

    gam.eta = retractionZ.X(Z=eta*tk, X=gam)
    Rval    = objGamE(gam=gam.eta, x=x, mu=mu, lam=lam, beta=beta, wt=wt) #- tval
    obj     = tval - Rval #+ abs[3]* tk * eta2
    m = m+1
  }

  if (obj < 0 ) tk = 0
  return(tk)
}


objGamE <- function( gam=NULL, x = NULL, mu=NULL, lam=NULL,  beta=NULL, wt=NULL ) {
  G   = length(beta)
  mal = sapply(1:G, function(k) { colSums( crossprod( gam,  t(x) - mu[k,])^2*1/lam[k,] )   })
  sum( rowMeans( t(wt)*t(mal)^beta) )
}


lamEBG1 <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  G = length(beta)
  p = ncol(x)
  N = nrow(x)

  gam = gam[,,1]
  maxBeta = max(beta)
  eta = sapply(1:G, function(k) {
    tr2 = crossprod( gam,  t(x) - mu[k,] )^2
    del = colSums( tr2 /lam[k,] )

    delb1 = del^(beta[k] - 1)
    delb1[del == 0] = 0

    colSums( (wt[,k]*delb1* t( tr2 * lam[k,]^(maxBeta-1)  )  ) )
  })
  lam = rowSums(eta*maxBeta*1/N)^(1/maxBeta)

  lam
}

getGamLamVVE <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {
  G = length(beta)
  p = ncol(x)
  N = nrow(x)

  Gmat = sapply(1:G, function(k) {
    tr  = t(x) - mu[k,]
    mal = colSums( crossprod( gam[,,k],  tr )^2*1/lam[k,] )

    malb1 = mal^(beta[k] - 1)
    malb1[mal == 0] = 0

    beta[k]*(tr %*% ( t(tr) * (wt[,k]*malb1 ) )  )
  } )
  Gmat = rowSums(Gmat)/N
  dim(Gmat) = c(p,p)

  Grad = Gmat %*% gam[,,1] %*% diag(1/lam[1,])
  gam  = gamEBG1(lam=lam, gam=gam[,,1], beta=beta, Grad=Grad, x=x, mu=mu, wt=wt, ng=ng)
  gam=array(gam, c(p,p,G) )


  lam = sapply(1:G, function(k) {
    tr2 = crossprod( gam[,,k], t(x) - mu[k,] )^2
    if (beta[k] > 1) lam = lamBG1( tr2, lam[k,], beta[k], ng[k], wt=wt[,k])
    else lam = lamBL1(tr2, lam[k,], beta[k], ng[k], wt[,k])
    lam
  })

  list(lam=t(lam), gam=gam )
}


getGamLamEEV <- function(x=NULL, mu=NULL, lam=NULL, gam=NULL, beta=NULL, ng=NULL, wt=NULL) {

  G = length(beta)
  p = ncol(x)
  N = nrow(x)

  gam = sapply(1:G, function(k) {
    tr  = t(x) - mu[k,]
    mal = colSums( crossprod( gam[,,k],  tr )^2*1/lam[k,] )

    malb1 = mal^(beta[k] - 1)
    malb1[mal == 0] = 0

    Gmat = tr %*% ( t(tr) * (beta[k]*wt[,k]/ng[k]*malb1 ) )

    gamBG1(lam=lam[k,], gam=gam[,,k], beta=beta[k], Gmat=Gmat, tr=tr, wt=wt[,k], ng=ng[k])
  } )
  dim(gam) = c(p,p,G)

  lam0 = lam
  if (any(beta > 1) ) {
    #lam = sapply(1:G, function(k) {
    # tr2 = crossprod( gam[,,k], t(x) - mu[k,] )^2
    # del = colSums( tr2/lam[k,] )
    # colMeans( wt[,k]*del^(beta[k]-1)*t(tr2) )*beta[k]
    #})
    #lam = rowSums(lam)
    maxBeta = max(beta)
    eta = sapply(1:G, function(k) {
      tr2 = crossprod( gam[,,k], t(x) - mu[k,] )^2
      del = colSums( tr2/lam[k,] )
      beta[k]*colMeans( (wt[,k]*del^(beta[k]-1)* t( tr2 * lam[k,]^(maxBeta-1)  )  ) )
    })
    lam = rowSums(eta)^(1/maxBeta)
  } else {
    lam = sapply(1:G, function(k) {
      tr2 = crossprod( gam[,,k], t(x) - mu[k,] )^2
      del = colSums( tr2/lam[k,] )

      delb1 = del^(beta[k]-1)
      delb1[del == 0] = 0

      colMeans( wt[,k]*delb1*t(tr2) )*beta[k]
    })
    lam = rowSums(lam)
  }

  list(lam=matrix(lam, G, p, byrow=TRUE), gam=gam )
}

getall <- function(loglik) {
  if (length(loglik) <3) stop("must have at least 3 likelihood values")
  n = length(loglik)
  lm1 = loglik[n]
  lm  = loglik[(n-1)]
  lm_1  = loglik[(n-2)]
  am = (lm1 - lm)/(lm - lm_1)
  lm1.Inf = lm + (lm1 - lm)/(1-am)
  val = lm1.Inf - lm
  if (is.nan(val)) val=0
  if (val < 0) val= 1
  return( val )
}




igpar <- function(data=NULL, G=NULL, m=10, n=10, model=NULL) {
  ### m is the number of random starts
  ### n is the number of EM iterations

  z = numeric(m+2)
  gparn <- list()
  for (i in 1:length(z) ) gparn[[i]] = list()

  getlogden   = create.logz(x=data, G=G, model=model)

  gparn[[1]] = igpar.kmeans(data=data, G=G,  n=n, model=model)
  z[1] = logden.to.loglik(getlogden(gparn[[1]]), gparn[[1]]$pi)

  gparn[[2]] = igpar.hclust(data=data, G=G,  n=n, model=model)
  z[2] = logden.to.loglik(getlogden(gparn[[2]]), gparn[[2]]$pi)

  for (i in 3:length(z)) {
    gparn[[i]] = igpar.wt(data=data, G=G, n=n, model=model)
    z[i] = logden.to.loglik(getlogden(gparn[[i]]), gparn[[i]]$pi)
  }

  remove = is.infinite(z)
  z      = z[!remove]
  gparn  = gparn[!remove]
  #print( cbind(order(z,decreasing=TRUE), z[order(z,decreasing=TRUE)] ) )
  gpar = gparn[[ order(z,decreasing=TRUE)[1] ]]

  return(gpar)
}





igpar.wt <- function(data=NULL, G=2, n=10, wt=NULL, model=NULL) {

  if (is.null(wt)) {
    wt = matrix(rexp(nrow(data)*G), nrow(data), G)
    wt  =wt/rowSums(wt)
  }

  gpar = rgpar(data, G, wt=wt, model=model)
  gpar.update = create.update(x=data, G=G, model=model)

  try({ for (i in 1:n) gpar = gpar.update(x=data, G=G, model=model )  }, silent=TRUE)

  return(gpar)
}


combinewk <- function(weights=NULL, label=NULL)	{
  # known is a numeric with
  # 0 if unknown group membership
  # 1,2,3,.. for label of known group
  if (is.null(label)) stop('label is null')
  kw     = label !=0
  for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
  return(weights)
}

igpar.lab <- function(data=NULL, G=2, n=10, lab=NULL, model=NULL) {
  z = combinewk(weights=matrix(0, nrow=nrow(data), ncol=G), label=lab)
  gpar = igpar.wt(data=data, G=G, n=n, wt=z, model=model)
  return(gpar)
}



igpar.kmeans <- function(data=NULL, G=NULL, n=1, nstart = 1, model=NULL) {
  lw   = kmeans(data, centers=G, nstart = nstart)$cluster
  gpar = igpar.lab(data=data, G=G, n=n, lab=lw, model=model)
  return(gpar)
}

igpar.hclust <-  function(data=NULL, G=NULL, n=50, model=NULL) {
  lw   = cutree(hclust(dist(data), "average"),k=G)
  gpar = igpar.lab(data=data, G=G, n=n, lab=lw, model=model)
  return(gpar)
}


ncovpar <- function(modelname=NULL, p=NULL, G=NULL) {
  if (is.null(p)) stop("p is null")
  if (is.null(G)) stop("G is null")
  if (is.null(modelname)) stop("modelname is null")

  if (modelname == "EII") npar = 1
  else if (modelname == "VII") npar = G
  else if (modelname == "EEI") npar = p
  else if (modelname == "VEI") npar = p + G -1
  else if (modelname == "EVI") npar = p*G - G +1
  else if (modelname == "VVI") npar = p*G
  else if (modelname == "EEE") npar = p*(p+1)/2
  else if (modelname == "EEV") npar = G*p*(p+1)/2 - (G-1)*p
  else if (modelname == "VEV") npar = G*p*(p+1)/2 - (G-1)*(p-1)
  else if (modelname == "VVV") npar = G*p*(p+1)/2
  else if (modelname == "EVE") npar = p*(p+1)/2 + (G-1)*(p-1)
  else if (modelname == "VVE") npar = p*(p+1)/2 + (G-1)*p
  else if (modelname == "VEE") npar = p*(p+1)/2 + (G-1)
  else if (modelname == "EVV") npar = G*p*(p+1)/2 - (G-1)
  else stop("modelname is not correctly defined")

  return(npar)
}


getInitialization <- function(initialization=NULL, G=NULL, x=NULL, model=NULL) {
  # initialization is a list then is assumed to be a list with the
  #         same format as rgpar
  # initialization can be a z matrix.
  # initialization can be a vector of labels
  # a single number indicating the number of random starts
  #          in additon to kmeans

  if ( is.null( initialization ) )  gpar  = rgpar(data=x, G=G, wt=NULL, model=model)
  else if ( is.list(initialization) ) gpar = initialization
  else if ( initialization[1] == 0) gpar = igpar.kmeans(data=x, G=G, n=10, model=model)
  else if ( length(initialization) == 1 & initialization[1] > 0) gpar = igpar(data=x, G=G, m=initialization, n=10, model=model)
  else if ( is.matrix(initialization) ) gpar = igpar.wt(data=x, G=G, n=10, wt=initialization, model=model)
  else if ( is.vector(initialization) ) gpar = igpar.lab(data=x, G=G, n=10, lab=initialization, model=model)
  else stop("Wrong initialization")

  gpar
}





EM <- function(data=NULL, initialization=NULL, G=2, max.iter=100, epsilon=1e-2, label=NULL, model=NULL) {

  logden.lab = create.logden.lab(n=nrow(data), G=G, label=label)

  gpar = getInitialization(initialization, G=G, x=data, model=model )
  getlogden   = create.logz(x=data, G=G, model=model)
  gpar.update = create.update(x=data, G=G, model=model)
  ##

  loglik = numeric(max.iter)
  for (i in 1:3) {
    logden    = getlogden(gpar) + logden.lab
    loglik[i] = logden.to.loglik(logden, gpar$pi)
    wt        = logden.to.weights(logden, gpar$pi)
    gpar      = gpar.update(gpar, wt)
  }

  while ( ( getall(loglik[1:i]) > epsilon) & (i < (max.iter) ) )  {
    i = i+1
    logden    = getlogden(gpar) + logden.lab
    loglik[i] = logden.to.loglik(logden, gpar$pi)
    wt        = logden.to.weights(logden, gpar$pi)
    gpar      = gpar.update(gpar, wt)
  }

  logden = getlogden(gpar) + logden.lab
  wt     = logden.to.weights(logden, gpar$pi)
  mapz   = MAP( wt )
  numpar = npar.model(model, ncol(data), G)

  val = list(loglik= loglik[1:i], gpar=gpar, z=wt, map=mapz, label=label, numpar = numpar, maxLoglik = loglik[i])
  return(val)
}






EMGr <- function(data=NULL, initialization=NULL, iModel="EIIE", G=2, max.iter=500, epsilon=1e-2, label=NULL, modelSet="all", skewness=FALSE, keepResults=FALSE, seedno=1, scale=TRUE) {
  ## G can be a range of values.
  ## iModel is the model used for initialization.
  ## all models are started using the same z matrix resutling from initialization
  ## if modelSet="all" then all the model are sued default
  set.seed(seedno)
  if(is.data.frame(data)) data <- data.matrix(data)
  if(scale==TRUE) data <- scale(data)
  if (length(initialization) != 1 ) stop("Initialization can only be the number of random starts for this function")

  if (modelSet[1] == "all") {
    mne = c("EIIE", "VIIE", "EEIE", "VVIE", "EEEE", "EEVE", "VVEE", "VVVE")
    mnv = c("EIIV", "VIIV", "EEIV", "VVIV", "EEEV", "EEVV", "VVEV", "VVVV")
    mn = c(mne, mnv)
  } else {
    mn = modelSet
  }

  if (skewness) mn = paste(mn, "V", sep="")

  nam = list( mn, G )
  num.iter   = matrix(NA, nrow=length(mn), ncol=length(G), dimnames =nam )
  num.par    = matrix(NA, nrow=length(mn), ncol=length(G), dimnames =nam )
  loglik     = matrix(NA, nrow=length(mn), ncol=length(G), dimnames =nam )
  increasing = matrix(NA, nrow=length(mn), ncol=length(G), dimnames =nam )

  result = list()
  for (k in 1:length(G)) {
    gpar0     = getInitialization(initialization, G=G[k], x=data, model=iModel)
    getlogden = create.logz(x=data, G=G[k], model=iModel)
    z0        = logden.to.weights(getlogden(gpar0), gpar0$pi)

    result[[G[k]]] = list()
    result[[G[k]]] = lapply(1:length(mn), function(j) {
      test = try({ EM(data=data, G=G[k], model= mn[j], label=label, initialization=z0, max.iter=max.iter, epsilon=epsilon) }, silent=TRUE)
      if (length(test) == 1) {
        test = list(time=NA, loglik=rep(1,2), numpar=NA, maxLoglik=NA  )
      }
      test
    } )

    increasing[,k] = sapply(result[[G[k]]], function(tem) { all( round(diff(tem$loglik), 12) >= 0) })
    num.iter[,k] = sapply(result[[G[k]]], function(tem) { length(tem$loglik) })
    num.par[,k]  = sapply(result[[G[k]]], function(tem) { tem$numpar })
    loglik[,k]   = sapply(result[[G[k]]], function(tem) { tem$maxLoglik })
    loglik[which(!increasing,arr.ind = TRUE)] <- NA
  }
  BIC    = -2*loglik + num.par*log(nrow(data))
  maxBIC = as.numeric(arrayInd(which.min(BIC), dim(BIC)))

  if (keepResults) val = list(allModels=result)
  else val = list()
  val = append(val, list( bestmod=result[[ G[maxBIC[2]] ]][[ maxBIC[1] ]],loglik=loglik, num.iter=num.iter,  num.par=num.par, BIC=BIC, maxBIC=maxBIC ) )
  class(val) <- "spemix"
  return(invisible(val))
}




rgpar <- function(data=NULL, G=NULL, model=NULL, wt=NULL) {
  p = ncol(data)
  n = nrow(data)

  if (is.null(wt)) {
    wt = matrix( rexp(n*G), nrow=n, ncol=G )
    wt = t(apply(wt, 1, function(z) {z/sum(z)} ))
  }

  val = list()
  # Gxp
  val$mu  = t(sapply(1:G, function(k) { apply(data,2,weighted.mean, wt[,k] )} ) )
  eigval  = sapply(1:G, function(k) { apply(data, 2, var) })
  eigvec  = sapply(1:G, function(k) { diag(p) } )

  if (substr(model,1,2) == "EI") val$lam = matrix( mean(eigval), G, ncol=p)
  else if (substr(model,1,2) == "VI") val$lam = matrix(apply( eigval,2, mean),G,p)
  else if (substr(model,1,2) == "VV") val$lam = t(eigval)
  else if (substr(model,1,2) == "EE") val$lam =  matrix(apply(eigval,1,mean), G, ncol=p, byrow=TRUE)
  else stop("The given model, ", model, " is not valid at position 2")

  if (substr(model,3,3) == "I") val$gam = NULL
  else if (substr(model,3,3) == "V") val$gam = array(eigvec, c(p,p,G))
  else if (substr(model,3,3) == "E") val$gam = array(eigen(cov.wt(data, method="ML")$cov)$vectors, c(p,p,G))
  else stop("The given model, ", model, " is not valid at position 3")

  if (substr(model,4,4) == "V" ) val$beta  = rep(1/2, G)
  else if (substr(model,4,4) == "E" ) val$beta = rep(1/2, G)
  else if (substr(model,4,4) == "D" ) val$beta = rep(1, G)
  else stop("The given model, ", model, " is not valid at position 4")

  #if (substr(model,5,5) == "V" ) val$eta  = matrix( runif(p), G, ncol=p)
  #else if (substr(model,5,5) == "N" ) val$eta = matrix( runif(p), G, ncol=p)
  #else stop("The given model, ", model, " is not valid at position 5")

  if (substr(model,5,5) != "") val$eta  = matrix( 0, G, ncol=p)


  val$pi = apply(wt,2,mean)
  val$model = model
  return(val)
}



create.logz <- function(x=NULL, G=NULL, model=NULL) {
  n = nrow(x)
  p = ncol(x)

  if (substr(model,3,3) == "I") malfn <- malI
  else malfn <- malG

  if (substr(model,5,5) == "") {
    tskewpFn <- function(eta=NULL, mu=NULL) { 0 }
  } else {
    tskewpFn <- function(eta=NULL, mu=NULL) {
      skew  = sapply(1:G, function(k) { eta[k,] %*% (t(x) - mu[k,]) })
      tskewp     = pnorm(t(skew), log.p=TRUE) + log(2)
    }
  }


  function(gpar) {
    ## n x G
    mal  = sapply(1:G, function(k) {
      malfn(x=x, mu=gpar$mu[k,], lam=gpar$lam[k,], gam=gpar$gam[,,k])
    }) #mahalanobis distance

    #skew  = sapply(1:G, function(k) {
    #  gpar$eta[k,] %*% (t(x) - gpar$mu[k,])
    #})
    #tskewp     = pnorm(t(skew), log.p=TRUE)
    tskewp = tskewpFn(gpar$eta, gpar$mu)

    tmalb     = -1/2*t(mal)^gpar$beta
    logd.beta = -lgamma(1 + p/(2*gpar$beta) ) - (1 + p/(2*gpar$beta) )*log(2)
    logC      = log(p) - p/2*log(pi) + lgamma(p/2)
    tlogden   = rowSums( -1/2*log(gpar$lam) ) + (tmalb + tskewp + logd.beta) + logC

    return(tlogden) #transpose of log density
  }
}




newmuSkewed <- function(mu=NULL, x=NULL, lam=NULL, gam=NULL, wt=NULL, ng=NULL, beta=NULL, getDel=NULL, getInvS=NULL, eta=NULL) {
  ## Netwon Raphson - Direct
  invS = getInvS( lam = lam, gam = gam)

  tr   = t(x) - mu
  z0   = eta %*% tr
  del  = getDel( tr, gam, lam )

  mu0= mu
  if (beta > 1) {
    skew.wt = exp( dnorm(z0, log=TRUE) - pnorm(z0, log.p = TRUE)   )

    mu.grsk = -1* sum( wt*(skew.wt) )*eta
    mu.grml =  colSums( ( t(tr)*beta*wt*del^(beta-1) ) %*% invS )

    mu.grad = (mu.grml + mu.grsk)/ng
    mu.cov = tr %*% ( t(tr)*wt*del^(beta - 2) )
    mu.H1  = sum(wt*del^(beta-1)) * invS
    mu.H2  = 2*(beta-1) *  ( invS %*% mu.cov  %*% invS )

    # 2nd derivative for pnorm
    #mu.H3  = sum( wt*(z0*skew.wt + skew.wt^2)  )*outer(eta,eta)
    # MM update
    mu.H3 = 1*sum( wt  )*outer(eta,eta)
    #mu.H  = -1* ( beta*(mu.H1+mu.H2) + mu.H3 )/ng
    #mu     = mu - solve(mu.H, mu.grad)

    mu.H  =   ( beta*(mu.H1+mu.H2) + mu.H3  )/ng
    mu    = mu + solve(mu.H, mu.grad)
  }  else {
    if (any(del == 0) ) {
      xInf = which(del == 0)
      mu = x[xInf,]
    } else {
      skew.wt = exp( dnorm(z0, log=TRUE) - pnorm(z0, log.p = TRUE) )

      #colSums( ( x*beta*wt*del^(beta-1) ) )/
      H1  = invS*sum( beta*wt*del^(beta-1) ) + sum( wt  )*outer(eta,eta)
      b1   = -1* sum( wt*(skew.wt + 1*z0) )*eta
      b2   = invS %*% colSums( x*beta*wt*del^(beta-1) ) + outer(eta,eta) %*% colSums( x*wt )
      mu   = solve(H1, b1 + b2)
      #mu = colSums( ( x*beta*wt*del^(beta-1) ) )/sum( beta*wt*del^(beta-1) )
    }
  }

  if ( any(is.nan(mu)) ) {
    stop("here - mu ")
  }
 # print(mu)
  return(mu)
}




newmu0Skewed <- function(mu=NULL, x=NULL, lam=NULL, gam=NULL, wt=NULL, ng=NULL, beta=NULL, getDel=NULL, getInvS=NULL, eta=NULL) {
  ## Netwon Raphson - Direct
  p = ncol(x)

  invS = getInvS( lam = lam, gam = gam)
  tr   = t(x) - mu
  del = getDel( tr, gam, lam )

  mu0 =mu
  if (beta > 1) {
    mu.grad =  colSums( ( t(tr)*beta*wt*del^(beta-1) ) %*% invS )

    mu.cov = tr %*% ( t(tr)*wt*del^(beta - 2) )
    mu.H1  = sum(wt*del^(beta-1)) * invS
    mu.H2  = 2*(beta-1) *  ( invS %*% mu.cov  %*% invS )

    mu     = mu - solve(-(beta/ng*(mu.H1+mu.H2) ), mu.grad/ng)
  } else {
    if (any(del == 0) ) {
      xInf = which(del == 0)
      mu = x[xInf,]
    } else {
      #mu =  colSums( ( x*beta*wt*del^(beta-1) ) )/sum( beta*wt*del^(beta-1) )
      wtdb = wt*del^(beta-1)
      wtdb = wtdb/max(wtdb)
      mu   = colSums( x*wtdb )/sum( wtdb )
      # print("mu = formula")
      # stop('here')
    }
  }
  if ( any(is.nan(mu)) ) {
    stop("here - mu ")
  }
  return(mu)
}





newEtaV <- function(mu=NULL, x=NULL, wt=NULL, ng=NULL, eta=NULL ) {
  G = nrow(mu)
  t(sapply( 1:G, function(k) {
    newEta(mu=mu[k,], x=x, wt=wt[,k], ng=ng[k], eta=eta[k,] )
  }))
}


newEta <- function(mu=NULL, x=NULL, wt=NULL, ng=NULL, eta=NULL) {
  ## Netwon Raphson - Direct
  tr   = t(x) - mu
  z0   = as.numeric(eta %*% tr)

  skew.wt = exp( dnorm(z0, log=TRUE) - pnorm(z0, log.p = TRUE)   )
  grad.eta = colSums( wt*(skew.wt+z0) *t(tr) )/ng
  H.eta = ( tr %*% (t(tr)*wt ) )/ng

  eta     = solve(H.eta, grad.eta)
  return(eta)
}


newEtaD <- function(mu=NULL, x=NULL, wt=NULL, ng=NULL, eta=NULL ) {
  ## Non-estimated or don't estimate this quantity
  eta
}


newEtaNULL <- function(mu=NULL, x=NULL, wt=NULL, ng=NULL, eta=NULL ) {
  ## Non-estimated or don't estimate this quantity
  ## For the elliptical or non-skewed model
  NULL
}




create.update  <- function(x=NULL, G=NULL, model=NULL) {
  n = nrow(x)
  p = ncol(x)

  if ( substr(model,3,3) == "I" ) {
    malfn   <- malI
    getInvS <- function(lam=NULL, gam=NULL) { diag(1/lam) } #invSigma
    getdel  <- function(tr=NULL, gam=NULL, lam=NULL ) { colSums( (tr)^2*1/lam ) }
  } else {
    malfn   <- malG
    getInvS <- function(lam=NULL, gam=NULL) { gam %*% diag(1/lam) %*% t(gam) } #invSigma
    getdel  <- function(tr=NULL, gam=NULL, lam=NULL ) { colSums( crossprod( gam,  tr )^2*1/lam ) }
  }

  if (substr(model,4,4) == "V") getBeta = newBetaV
  else if (substr(model,4,4) == "D") getBeta = newBetaD #fixed
  else if (substr(model,4,4) == "E") getBeta = newBetaE
  else stop("The model is ", substr(model,4,4) )

  if (substr(model,1,3) == "VII") getGamLam = getGamLamVII
  else if (substr(model,1,3) == "EII") getGamLam = getGamLamEII
  else if (substr(model,1,3) == "VVI") getGamLam = getGamLamVVI
  else if (substr(model,1,3) == "EEI") getGamLam = getGamLamEEI
  else if (substr(model,1,3) == "VVV") getGamLam = getGamLamVVV
  else if (substr(model,1,3) == "EEE") getGamLam = getGamLamEEE
  else if (substr(model,1,3) == "VVE") getGamLam = getGamLamVVE
  else if (substr(model,1,3) == "EEV") getGamLam = getGamLamEEV
  else stop("The model is ", substr(model,1,3) )


  if (substr(model,5,5) == "D") getEta = newEtaD
  else if (substr(model,5,5) == "") getEta = newEtaNULL
  else getEta = newEtaV

  if (substr(model,5,5) == "") newMu = newmu0Skewed
  else newMu = newmuSkewed

  function(gpar=NULL, wt=NULL) {
    ng = colSums(wt)
    gpar$pi = ng/n

    gpar$mu = t(sapply(1:G, function(k) {
      newMu(mu=gpar$mu[k,], x=x, lam = gpar$lam[k,], gam = gpar$gam[,,k], wt=wt[,k], ng=ng[k], beta=gpar$beta[k], getDel=getdel, getInvS=getInvS, eta=gpar$eta[k,] )
    }))

    gpar$eta = getEta(gpar$mu, x, wt, ng, gpar$eta)

    tem = getGamLam(x=x, mu=gpar$mu, lam=gpar$lam, gam=gpar$gam, beta=gpar$beta, ng=ng, wt=wt)
    gpar$lam = tem$lam
    gpar$gam = tem$gam

    mal = sapply(1:G, function(k) {  malfn(x, gpar$mu[k,], gpar$lam[k,], gpar$gam[,,k] )  }   )
    gpar$beta = getBeta(beta=gpar$beta, p=p, wt=wt, ng=ng, del=mal)

    return(gpar)
  }
}




npar.model <- function(modelname=NULL, p=NULL, G=NULL) {
  val = numeric(3)
  val[1] = G-1  ## pi
  val[2] = G*p  ## mu
  val[3] = ncovpar(modelname= substr(modelname, 1, 3) , p=p, G=G)

  if (substr(modelname, 4, 4) == "V") val[4] = G
  else val[4] = 1

  if (substr(modelname,5,5) != "") val[5] = G * p
  else  val[5] = 0

  val = sum(val)
  return(val)
}



