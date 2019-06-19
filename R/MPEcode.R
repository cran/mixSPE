print.pemix <- function(x, ...){
  cat("Call:\n")
  print.default(x$call)
  cat("-----------------------------------\n")
  cat(x$msc[which.max(x$msc[,2]),1], "groups and scale structure",
      attr(which.max(x$msc[,2]),"names"), "was selected by the BIC.","\n",
      "BIC was",x$msc[which.max(x$msc[,2]),2],"\n",x$msc[which.max(x$msc[,2]),3],
      "parameters were estimated.\n","Log-likelihood calculated at MLEs was",
      x$msc[which.max(x$msc[,2]),4],"\n")
  cat(x$msc[which.max(x$msc[,5]),1], "groups and scale structure",
      attr(which.max(x$msc[,5]),"names"), "was selected by the ICL.","\n",
      "ICL was",x$msc[which.max(x$msc[,5]),5],"\n",x$msc[which.max(x$msc[,5]),3],
      "parameters were estimated.\n","Log-likelihood calculated at MLEs was",
      x$msc[which.max(x$msc[,5]),4],"\n")
  cat("-----------------------------------\n")
  cat("Time taken:", x$time[3], "seconds.\n")
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
mpe <- function(verbose = FALSE, dat = NULL, seedno = 1, G = 1:4,
                         start = "kmeans", kmeansinit=10,
                eps=0.005, maxit=5000, label = NULL, modelnames = c("EIIE", "VIIE", "EEIE", "VVIE",
                                                        "EEEE", "EEVE", "VVEE", "VVVE",
                                                        "EIIV", "VIIV", "EEIV", "VVIV",
                                                        "EEEV", "EEVV", "VVEV", "VVVV")){
  #main workhorse function; fits mixture of PE to data.
  #dat is data, seedno fixes the seed---important for initialization, esp. random,
  #start controls what kind of start: random or kmeans
  #eps is used in the Aitken's stopping criterion
  #maxit is the maximum number of GEM iterations allowed
  #anneal controls the annealing
  #label gives the already labelled observations
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
  model.type <- function(modelname=NULL) {
    if (is.null(modelname)) stop("modelname is null")
    if (substr(modelname,1,3) == "EII") val = msEII()
    else if (substr(modelname,1,3) == "VII") val = msVII()
    else if (substr(modelname,1,3) == "EEI") val = msEEX()
    else if (substr(modelname,1,3) == "VVI") val = msVVX()
    else if (substr(modelname,1,3) == "EEE") val = msEEX()
    else if (substr(modelname,1,3) == "EEV") val = msEEX()
    else if (substr(modelname,1,3) == "VVE") val = msVVX()
    else if (substr(modelname,1,3) == "VVV") val = msVVV()
    else stop("modelname or covtype is not correctly defined")
    if (!is.list(val)) val = list(sigma=val)
    return(val)
  }
  npar.model <- function(modelname=NULL, p=NULL, G=NULL) {
    val = NULL
    val[1] = G-1 #for pi
    val[2] = G*p #for mean for data
    val[3] = ncovpar(modelname= modelname, p=p, G=G) #Covariance for data
    if(substr(modelname,4,4)=="V") val[4] = G  #for unconstrained beta_g
    if(substr(modelname,4,4)=="E") val[4] = 1  #for constrained beta_g
    val = sum(val)
    return(val)
  }
  ncovpar <- function(modelname=NULL, p=NULL, G=NULL) {
    if (is.null(p)) stop("p is null")
    if (is.null(G)) stop("G is null")
    if (is.null(modelname)) stop("modelname is null")
    if (substr(modelname,1,3) == "EII") npar = 1
    else if (substr(modelname,1,3) == "VII") npar = G
    else if (substr(modelname,1,3) == "EEI") npar = p
    else if (substr(modelname,1,3) == "VVI") npar = p*G
    else if (substr(modelname,1,3) == "EEE") npar = p*(p+1)/2
    else if (substr(modelname,1,3) == "EEV") npar = G*p*(p+1)/2 - (G-1)*p
    else if (substr(modelname,1,3) == "VVE") npar = p*(p+1)/2 + (G-1)*p
    else if (substr(modelname,1,3) == "VVV") npar = G*p*(p+1)/2
    else stop("modelname is not correctly defined")
    return(npar)
  }
  map <- function(x){
    return(apply(X = x, MARGIN = 1, FUN = which.max))
  }
  unmap <- function(mapz, git){
    zunmap <- matrix(data = 0, nrow = length(mapz), ncol = git)
    alphabet <- 1:git
    for(u in 1:git){
      zunmap[which(mapz == alphabet[u]), u] <- 1
    }
    return(zunmap)
  }
  fnz <- function(label) {
    if(it>1){ zold<-z }
    z_g <- matrix(NA, nrow = N, ncol = git)
    fz <- z
    xpart <- matrix(NA, N, git)
    kpart <- matrix(NA, N, git)
    for (g in 1:git) {
      xpart[, g] <- try(det(Sig[[g]])^(-0.5) * exp(-0.5 * (mahalanobis(as.matrix(dat), center = mu[[g]],
                                                                       cov = invSig[[g]], inverted = TRUE))^betag[g]), silent = TRUE)
      kpart[, g] <- try(p * (gamma(p/2))/(pi^(p/2) * (gamma(1 + p/(2 * betag[g]))) * 2^(1 + p/(2 * betag[g]))),
                        silent = TRUE)
      if (is.numeric(kpart[, g]) && is.numeric(xpart))
        z_g[, g] <- pi_g[g] * kpart[, g] * xpart[, g]
      if (!is.numeric(z_g[, g]))
        z_g[, g] <- Inf
    }
    if(any(!is.na(label))) z_g[!is.na(label),] <- unmap(label[!is.na(label)],git) #for classification
    fz <- try(z_g/rowSums(z_g), silent = TRUE)
    if (!is.numeric(fz)) fz <- Inf
    return(list(z = fz, forz = z_g, zold = zold))
  }
  stopcheck <- function() {
    if ((it > (maxit - 1)) || (incparll[it] == Inf) || (incparll[it] == incparll[it - 1])) {stop <- 1}
    if (incparll[it] != Inf && it < maxit) {
      if (incparll[it] != incparll[it - 1]) {
        ait1less <- (incparll[it] - incparll[(it - 1)])/(incparll[(it - 1)] - incparll[(it - 2)])
        llinf <- incparll[(it - 1)] + (1/(1 - ait1less)) * (incparll[it] - incparll[(it - 1)])
        if (it > 3) {
          stop <- 0 #in case the condition is not satisfied
          if (abs(llinf - incparll[(it - 1)]) < eps) {stop <- 1}
        }
      }
    }
    t<-na.omit(incparll[1:it])
    if(!all(t == cummax(t))) {stop <- 1}
    return(stop)
  }
  ms_x <- function() {
    condition<-1
    muprev<-list()
    if(it>1){
      pi_gold<-pi_g
      muprev<-mu
    }
    pi_g <- colMeans(as.matrix(z))
    munew <- list()
    mun <- rep(NA,p)
    for (g in 1:git) {
      mutrial=rep(NA,p)
      if (it == 1) {
        temp = cov.wt(as.matrix(dat), wt = (z[, g]), center = TRUE, method = "ML")
        muold <- temp$center
      }
      if (it > 1) {
        muold <- mu[[g]]
      }
      Sg<-mahalanobis(m_g[[g]],center=rep(0,p),cov=invSig[[g]],inverted=TRUE)

      fun_deriv<-function(mun){betag[g]*rowSums(sapply(1:N,function(i){as.numeric(z[i,g]*Sg[i]^(betag[g]-1))*
                                                                         invSig[[g]]%*%(m_g[[g]][i,])}))}
      fun_hess<-function(mun){matrix(rowSums(sapply(1:N,
                                                    function(i){as.numeric(-z[i,g]*betag[g]*(betag[g]-1)*
                                                                             Sg[i]^(betag[g]-2))*(-1)*invSig[[g]]%*%(dat[i,]-mun)%*%
                                                                  t((-2)*invSig[[g]]%*%(dat[i,]-mun))+as.numeric(-z[i,g]*
                                                                                                                    betag[g]*Sg[i]^(betag[g]-1))*invSig[[g]] })) ,p,p)
      }
      mutrial <- try(muold - as.vector(solve(fun_hess(muold))%*%fun_deriv(muold)),silent=TRUE)
      if(all(is.finite(mutrial)))  {munew[[g]]<-mutrial}
      if(!all(is.finite(mutrial))) {
        munew[[g]]<-rep(NA,p)
        condition<-0
      }
    }
    result <- list(pi_g = pi_g, mu = munew, pi_gold = pi_g, muold = muprev, con=condition)
    return(result)
  }
  ms_betagV <- function() {
    if(it>1){betaold<-betag}
    condition<-1
    betanew = rep(NA,git)
    betan <- rep(NA,git)
    for (g in 1:git) {
      betan<-betag[g]
      Sg<-mahalanobis(m_g[[g]],center=rep(0,p),cov=invSig[[g]],inverted=TRUE)

      fun_be <- function(betan) {   #derivative
        p/2 * colSums(z)[g]/betan^2 * digamma(1 + p/(2 * betan)) +
          p/2 * log(2)/betan^2 * colSums(z)[g] + sum(-z[,g]/2*Sg^betan*as.numeric(log(Sg)))
      }
      fun_secbe <- function(betan) {    #second derivative
        -p * colSums(z)[g]/betan^3 * digamma(1 + p/(2 * betan)) +
          p/2 * colSums(z)[g]/betan^2 * trigamma(1 + p/(2 * betan)) * (-p/(2 * betan^2)) +
          (-p)/2 * log(2)/betan^3 * colSums(z)[g] + sum(-z[,g]/2*Sg^betan*as.numeric(log(Sg))^2)
      }
      betanew[g] <- try(betan - fun_be(betan)/fun_secbe(betan),silent=TRUE)
      if(is.finite(betanew[g])){
        if (betanew[g]<0.05){betanew[g]<-0.05} else if(betanew[g]>200){betanew[g]<-200}
      }
      if(!is.finite(betanew[g])){
        betanew[g] <- NA
        condition<-0
      }
    }
    result <- list(betag = betanew, betaold = betaold, con=condition)
    return(result)
  }
  ms_betagE <- function() {
    betanew <- NA
    if(it>1) betaold <- betag
    condition <- 1
    fn_be <- function(betan){
      val <- numeric(git)
      valse <- numeric(git)
      for(g in 1:git){
        Sg <- mahalanobis(m_g[[g]],center=rep(0,p),cov=invSig[[g]],inverted=TRUE)
        val[g] <- 0.5*sum(z[,g]*Sg^betan*as.numeric(log(Sg)) )
        valse[g]<- 0.5*sum(z[,g]*Sg^betan*(as.numeric(log(Sg)))^2 )
      }
      return(list(befi = sum(val),besec = sum(valse)))
    }
    fun_be<-function(betan){
      N*p/(2*betan^2)*(digamma(1+p/(2*betan)) + log(2)) - fn_be(betan)$befi
    }
    fun_secbe<-function(betan){
      -p * N/betan^3 * digamma(1 + p/(2 * betan)) +
        p/2 * N/betan^2 * trigamma(1 + p/(2 * betan)) * (-p/(2 * betan^2)) +
        (-p) * log(2)/betan^3 * N - fn_be(betan)$besec
    }
    ifelse(p>1, betan <- betaold[1], betan <- betaold)
    betanew <- try(betan - fun_be(betan)/fun_secbe(betan),silent=TRUE)
    if(is.finite(betanew)){
      if (betanew<0.05){betanew<-0.05} else if(betanew>200){betanew<-200}
    }
    if(!is.finite(betanew)){
      betanew <- NA
      condition<-0
    }
    result <- list(betag = rep(betanew,git), betaold = betaold, con = condition)
    return(result)
  }
  msEII <- function() {
    condition<-1
    if(it>1){Sigold<-Sig}
    if(it==1){Sig[[1]]=diag(1)}
    Sigret <- list()
    mg<-array(NA,dim=c(N,p,git))
    for(g in 1:git){mg[,,g]<-z[,g]^(0.5/betag[g])*m_g[[g]]}
    fun_EII <- function(sigma) {
      p * N - sum(sapply(1:git, function(g) {
        betag[g] * sigma^(-betag[g]) * sum(rowSums(mg[,,g]*mg[,,g])^betag[g])}))
    }
    res <- try(uniroot(fun_EII, c(0.0001, (diag(Sig[[1]])[1]+100)))$root, silent = TRUE)
    if(is.finite(res)) {for(g in 1:git){Sigret[[g]]<-res*diag(p)}}
    if(!is.finite(res)) {
      Sigret[[g]]<-diag(NA,p)
      condition<-0}
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition)
    return(result)
  }
  msVII <- function() {
    condition<-1
    if(it>1){Sigold<-Sig}
    Signew <- rep(NA,git)
    Sigret <- list()
    for (g in 1:git) {
      Signew[g] <- try((betag[g]*sum(z[,g]*rowSums(m_g[[g]]*m_g[[g]])^betag[g])/(p*colSums(z)[g]))
                       ^(1/betag[g]),silent=TRUE)
      if(is.finite(Signew[g])) Sigret[[g]]<-Signew[g]*diag(p)
    }
    if(!all(is.finite(Signew))) {
      for(g in 1:git) Sigret[[g]]<-diag(NA,p)
      condition<-0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition)
    return(result)
  }
  msEEX <- function() {
    if(it>1){Sigold<-Sig}
    condition <- 1
    Sigret <- list()
    Anew <- list()
    Dnew <- list()
    Ddata <- list()
    DdataEbless1 <- list()
    if(it==1){ for(g in 1:git) {
      Anew[[g]] <- rep(1,p)
      if(!substr(modelname,3,3)=="I") Dnew[[g]] <- diag(p)
    }  }
    if(it>1){  for(g in 1:git) {
      Anew[[g]] <- eigen(Sig[[g]])$values
      if(substr(modelname,3,3)=="E") Dnew[[g]] <- eigen(Sig[[1]])$vectors
      if(substr(modelname,3,3)=="V") Dnew[[g]] <- eigen(Sig[[g]])$vectors
    }  }
    if(substr(modelname,3,3)=="I") {
      for(g in 1:git) Ddata[[g]] <- (z[,g])^(1/(2*betag[g]))*m_g[[g]]
      for(g in 1:git) DdataEbless1[[g]] <- m_g[[g]]
    }
    if(!substr(modelname,3,3)=="I") {
      for(g in 1:git) Ddata[[g]] <- (z[,g])^(1/(2*betag[g]))*m_g[[g]]%*%Dnew[[g]]
      for(g in 1:git) DdataEbless1[[g]] <- m_g[[g]]%*%Dnew[[g]]
    }
    if(!substr(modelname,3,3)=="I") {for(g in 1:git) DdataEbless1[[g]] <- m_g[[g]]%*%Dnew[[g]]}
    if(all(betag<1))  Anew[[1]] <- try(newaEbless1(s0=Anew, dat=DdataEbless1, beta=betag),silent=TRUE)
    else if(any(betag>=1)) {
      fnres<-newaE(s0=Anew, dat= Ddata, beta=betag)
      Anew[[1]] <- fnres$Ap
      condition<-fnres$con
    }
    if(all(is.finite(Anew[[1]]))){
      if(git>1) { for(g in 2:git) {Anew[[g]] <- Anew[[1]] } }
      if(all(Anew[[1]]>0) && !substr(modelname,3,3)=="I") {
        if(substr(modelname,3,3)=="E"){
          Dnew[[1]] <- newDbE(D=Dnew, Ak=Anew, Rdata=dat,  tmax = 100, b=betag)
          if(git>1) { for(g in 2:git) {Dnew[[g]] <- Dnew[[1]] } } }
        if(substr(modelname,3,3)=="V"){
          Dnew[[g]] <- newDb(D = Dnew[[g]], Ak = Anew[[g]], Rdata = dat, tmax = 100, b = betag[g], g = g)
        }
        for(g in 1:git) Sigret[[g]] <- Dnew[[g]]%*%diag(Anew[[g]])%*% t(Dnew[[g]])
      }
      if(all(Anew[[1]]>0) && substr(modelname,3,3)=="I") {
        for(g in 1:git) Sigret[[g]] <- diag(Anew[[g]])
      }
    }
    else if(!all(is.finite(Anew[[1]])) || !all(Anew[[1]]>0)){
      for(g in 1:git) Sigret[[g]]<-diag(NA,p)
      condition<-0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition, D=Dnew, A=Anew)
    return(result)
  }
  msVVX <- function() {
    if(it>1){Sigold<-Sig}
    condition<-1
    Sigret <- list()
    Anew <- list()
    Afin <- list()
    Dnew <- list()
    Dfin <- list()
    Ddata <- list()
    if(it==1){ for(g in 1:git) {
      Anew[[g]] <- rep(1,p)
      if(substr(modelname,3,3)=="E") Dnew[[g]] <- diag(p)
    }   }
    if(it>1){
      for(g in 1:git){ Anew[[g]] <- eigen(Sig[[g]])$values
                       if(substr(modelname,3,3)=="E") Dnew[[g]] <- eigen(Sig[[1]])$vectors }
    }
    if(substr(modelname,3,3)=="I"){
      for(g in 1:git) Ddata[[g]] <- (z[,g])^(1/(2*betag[g]))*m_g[[g]]
    }
    if(substr(modelname,3,3)=="E"){
      for(g in 1:git) Ddata[[g]] <- (z[,g])^(1/(2*betag[g]))*m_g[[g]]%*%Dnew[[g]]
    }
    for(g in 1:git) {
      if(betag[g]>=1){
        fnres <- newa(s0=Anew[[g]], dat= Ddata[[g]], beta=betag[g],ng=sum(z[,g]))
        if(fnres$con==1) Anew[[g]] <- fnres$sig
        if(fnres$con==0) Anew[[g]] <- diag(NA,p)
      }
      if(betag[g]<1){
        if(substr(modelname,3,3)=="E"){Anew[[g]] <- newaV(s0=Anew[[g]], dat= m_g[[g]]%*% Dnew[[g]],
                                                          beta=betag[g], g=g, ng=sum(z[,g]))}
        if(substr(modelname,3,3)=="I"){Anew[[g]] <- newaV(s0=Anew[[g]], dat= m_g[[g]], beta=betag[g],
                                                          g=g, ng=sum(z[,g]))}
      }
    }
    if(substr(modelname,3,3)=="E"){
      if(all(is.finite(unlist(Anew))) && all(unlist(Anew)>0)) {
        Dnew[[1]] <- newDbE(D=Dnew, Ak=Anew, Rdata=dat,  tmax = 100, b=betag)
        if(git>1) { for(g in 2:git) {Dnew[[g]] <- Dnew[[1]] } }
        for(g in 1:git) Sigret[[g]] <- Dnew[[g]]%*%diag(Anew[[g]])%*% t(Dnew[[g]])
      } }
    if(substr(modelname,3,3)=="I"){
      if(all(is.finite(unlist(Anew))) && all(unlist(Anew)>0)) {
        for(g in 1:git) Sigret[[g]] <- diag(Anew[[g]])
      } }
    if(!all(is.finite(unlist(Anew))) || !all(unlist(Anew)>0)) {
      for(g in 1:git) Sigret[[g]]<-diag(NA,p)
      condition<-0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con = condition, D=Dnew, A=Anew)
    return(result)
  }
  msVVV <- function() {
    if(it>1){Sigold<-Sig}
    condition<-1
    Sigret <- list()
    Anew <- list()
    Dnew <- list()
    for (g in 1:git){
      if(betag[g]<1){
        maha<-mahalanobis(m_g[[g]],center=rep(0,p),cov=invSig[[g]],inverted=TRUE)
        Sigret[[g]] <- 1/colSums(z)[g]*betag[g]*matrix(rowSums(sapply(1:N,function(i){z[i,g]*
                                                                                        maha[i]^(betag[g]-1)*(m_g[[g]][i,])%*%t(m_g[[g]][i,])})),p,p)
      }
      if(betag[g]>=1){
        if(it==1){
          Anew[[g]] <- rep(1,p)
          Dnew[[g]] <- diag(p)
        }
        if(it>1){
          Anew[[g]] <- eigen(Sig[[g]])$values
          Dnew[[g]] <- eigen(Sig[[g]])$vectors
        }
        Ddata <- list()
        Ddata[[g]] <- (z[,g])^(1/(2*betag[g]))*m_g[[g]]%*%Dnew[[g]]
        fnres <- newa(s0=Anew[[g]], dat= Ddata[[g]], beta=betag[g],ng=sum(z[,g]))
        if(fnres$con==1) Anew[[g]] <- fnres$sig
        if(fnres$con==0) Anew[[g]] <- diag(NA,p)
        if(all(is.finite(Anew[[g]])) && all(Anew[[g]]>0)) {
          Dnew[[g]] <- newDb(D=Dnew[[g]], Ak=Anew[[g]], Rdata=dat,  tmax = 100, b=betag[g],g=g)
          Sigret[[g]] <- Dnew[[g]]%*%diag(Anew[[g]])%*% t(Dnew[[g]]) }
      }
    }
    if(!all(is.finite(unlist(Anew))) || !all(unlist(Anew)>0)) {
      for(g in 1:git) Sigret[[g]]<-diag(NA,p)
      condition<-0
    }
    result <- list(Sig = Sigret, Sigold = Sigold, con=condition)
    return(result)
  }
  varfn <- function(scmat) {
    Cov <- list()
    for (g in 1:git) {
      Cov[[g]] <- 2^(1/betag[g]) * gamma((p + 2)/(2 * betag[g]))/(p * gamma(p/(2 * betag[g]))) * scmat[[g]]
    }
    return(Cov)
  }
  valsave <- function(){
    return(list(z=ff$zold, betag=bg$betaold, pi_g=mst_x$pi_gold,mu=mst_x$muold, Sig=mstep_Sig$Sigold))
  }
  testgrad.D = function(D=NULL, d=NULL, dat=NULL, Ak=NULL, b=NULL,g=NULL) {
    Dak = sweep(D, 2, 1/Ak, FUN="*")
    invS = D %*% diag(1/Ak) %*% t(D)
    dat_x<-sweep(dat,2,mu[[g]],"-")*(z[,g])^(1/(2*betag[g]))
    q = mahalanobis(x=dat_x, center= rep(0, d), cov=invS, inverted = TRUE)
    qb_1 = as.numeric(q)^(b-1)
    qb_1[which(q<=0)] <- 0

    Rk = cov.wt(x=dat_x, wt= qb_1, center= rep(0, d), method="ML")$cov*sum(qb_1)
    val = (2*b*Rk) %*% Dak
    return( val )
  }
  testgrad.DE = function(D=NULL, d=NULL, dat=NULL, Ak=NULL, b=NULL) {
    Rk<-list()
    Dak <- list()
    vl <- list()
    for(g in 1:git){
      Dak[[g]] = sweep(D[[g]], 2, 1/Ak[[g]], FUN="*")
      invS = D[[g]] %*% diag(1/Ak[[g]]) %*% t(D[[g]])
      dat_x<-sweep(dat,2,mu[[g]],"-")*(z[,g])^(1/(2*b[g]))
      q = mahalanobis(x=dat_x, center= rep(0, d), cov=invS, inverted = TRUE)
      qb_1 = as.numeric(q)^(b[g]-1)
      qb_1[which(q<=0)] <- 0

      Rk[[g]] = cov.wt(x=dat_x, wt= qb_1, center= rep(0, d), method="ML")$cov*sum(qb_1)
      vl[[g]] = (b[g]*Rk[[g]] %*% Dak[[g]])
    }
    val = 2*Reduce("+",vl)
    return( val )
  }
  testvalb = function(dat =NULL, Ak=NULL, D=NULL, b=NULL,g=NULL) {
    d = ncol(dat)
    Dak = sweep(D, 2, 1/Ak, FUN="*")
    invS = Dak %*% t(D)
    dat_x<-sweep(dat,2,mu[[g]],"-")*(z[,g])^(1/(2*betag[g]))
    q = mahalanobis(x= dat_x, center= rep(0, d), cov=invS, inverted = TRUE)
    val = sum(q^b)
    return(val)
  }
  testvalbE = function(dat =NULL, Ak=NULL, D=NULL, b=NULL) {
    d = ncol(dat)
    vl<-list()
    for(g in 1:git){
      Dak = sweep(D[[1]], 2, 1/Ak[[g]], FUN="*")
      invS = Dak %*% t(D[[1]])
      dat_x<-sweep(dat,2,mu[[g]],"-")*(z[,g])^(1/(2*b[g]))
      q = mahalanobis(x= dat_x, center= rep(0, d), cov=invS, inverted = TRUE)
      q[which(q<0)]<-0
      vl[[g]] = sum(q^b[g])
    }
    val = Reduce("+",vl)
    return(val)
  }
  newDb = function(D=NULL, Rdata=NULL, Ak=NULL, tmax = 100, b=NULL,g=NULL) {
    d = ncol(Rdata)
    eta   = -1*projZ.X(Z=testgrad.D(D=D, d=d, dat = Rdata, Ak=Ak, b=b,g=g), X=D)
    tval0 = testvalb(dat = Rdata, Ak=Ak, D=D, b=b,g=g)
    tk    = armijo(abs=c(1, .95, 1e-8 ), eta=eta, D=D,d=d, dat = Rdata, Ak=Ak, tval=tval0, mmax= 25, b=b,g=g )
    ifelse(tk == 0,  Xk1 <- D, Xk1 <- retractionZ.X(Z=eta*tk, X=D)  )
    return( Xk1 )
  }
  newDbE = function(D=NULL, Rdata=NULL, Ak=NULL, tmax = 100, b=NULL) {
    d = ncol(Rdata)
    eta   = -1*projZ.X(Z=testgrad.DE(D=D, d=d, dat = Rdata, Ak=Ak, b=b), X=D[[1]])
    tval0 = testvalbE(dat = Rdata, Ak=Ak, D=D, b=b)
    tk    = armijoE(abs=c(1, .95, 1e-8 ), eta=eta, D=D,d=d, dat = Rdata, Ak=Ak, tval=tval0, mmax= 25, b=b)
    if (tk == 0)  Xk1 = D[[1]]
    else Xk1 = retractionZ.X(Z=eta*tk, X=D[[1]])
    return( Xk1 )
  }
  newaE <- function(s0=NULL, dat=NULL, beta=NULL) {
    compsig <- list()
    condition <- 1
    for(g in 1:git){
      w  = dat[[g]] %*% diag(1/sqrt( s0[[g]] )  )
      ww = apply(w^2, 1, sum)^( beta[g]-1)
      ww[which(!is.finite(ww))]<-0
      if(any(ww>0)){
        d = ncol(dat[[g]])
        lam = cov.wt(w, wt= ww, center=rep(0,d),method='ML')$cov*sum(ww)*beta[g]
        compsig[[g]] = diag( diag( s0[[g]]^(max(beta)/2) ) %*% lam %*%
                               diag(s0[[g]]^(max(beta)/2)) )
      }
      Ap <- (1/N*Reduce('+', compsig))^(1/max(beta))}
    if(!any(ww>0)) {Ap<-diag(p)
                    condition <- 0 }
    return(list(Ap=Ap,con=condition)) # matrix A raised to a power
  }
  newaEbless1 <- function(s0=NULL, beta=NULL, dat=NULL) {
    compsig <- list()
    maha <- list()
    for(g in 1:git){maha[[g]]<-mahalanobis(dat[[g]],center=rep(0,p),cov=diag(1/s0[[g]]),inverted=TRUE)}
    for(g in 1:git) {
      compsig[[g]]<-diag(beta[g]*matrix(rowSums(sapply(1:N,function(i){
        z[i,g]*(maha[[g]][i])^(beta[g]-1)*(dat[[g]][i,])%*%t(dat[[g]][i,])
      })),p,p))
    }
    A <- 1/N*Reduce('+', compsig)
    return(A)
  }
  newaV <- function(s0=NULL, beta=NULL, dat=NULL, g=g, ng=ng) {
    maha <- mahalanobis(dat,center=rep(0,p),cov=diag(1/s0),inverted=TRUE)
    mahab <- maha^(beta-1)
    mahab[which(maha==0)] <- 0
    compA<-diag(beta/ng*matrix(rowSums(sapply(1:N,function(i){
      z[i,g]*mahab[i]*(dat[i,])%*%t(dat[i,])    })),p,p))
    return(compA)
  }
  newa <- function(s0=NULL, dat=NULL, beta=NULL,ng=NULL) {
    condition<-1
    w  = dat %*% diag(1/sqrt( s0 )  )
    ww = apply(w^2, 1, sum)^( beta-1)
    ww[which(!is.finite(ww))]<-0
    if(any(ww>0)){
      d = ncol(dat)
      lam = cov.wt(w, wt= ww, center=rep(0,d),method='ML')$cov*sum(ww)*beta/ng
      sig = diag( diag( s0^(beta/2) ) %*% lam %*% diag(s0^(beta/2)) )^(1/beta)
    }
    if(!any(ww>0)) {sig<-diag(p)
                    condition <- 0 }
    return(list(sig=sig,con=condition))
  }
  projZ.X = function(Z=NULL,X=NULL) {
    XZ = t(X) %*% Z
    return( Z - X %*% ( XZ + t(XZ) )/2 )
  }
  armijo = function(abs=c(1, .1, .1), eta=NULL, D=NULL, d=NULL, dat=NULL, Ak=NULL,
                    tval=NULL, mmax=100, b=NULL,g=NULL) {
    # abs is a numeric vector with alpha > 0, beta, sigma in (0,1)
    eta2  = sum(eta^2)
    m   = 1
    obj = -1
    tkseq = seq(0, log(1e-8)/log(abs[2]), length.out=mmax)
    while (obj < 0 & m < mmax  ) {
      tk   = exp( tkseq[m]*log(abs[2])) * abs[1]
      Reta = retractionZ.X(Z=eta*tk, X=D)
      Rval = testvalb(D= Reta, dat=dat, Ak=Ak, b=b,g=g) #- tval
      obj  = tval - Rval #+ abs[3]* tk * eta2
      m = m+1
    }
    if (obj < 0 ) tk = 0
    return(tk)
  }
  armijoE = function(abs=c(1, .1, .1), eta=NULL, D=NULL, d=NULL, dat=NULL, Ak=NULL,
                     tval=NULL, mmax=100, b=NULL) {
    # abs is a numeric vector with alpha > 0, beta, sigma in (0,1)
    eta2  = sum(eta^2)
    Reta = list()
    m   = 1
    obj = -1
    tkseq = seq(0, log(1e-8)/log(abs[2]), length.out=mmax)
    while (obj < 0 & m < mmax  ) {
      tk   = exp( tkseq[m]*log(abs[2])) * abs[1]
      Reta[[1]] = retractionZ.X(Z=eta*tk, X=D[[1]])
      Rval = testvalbE(D= Reta, dat=dat, Ak=Ak, b=b) #- tval
      obj  = tval - Rval #+ abs[3]* tk * eta2
      m = m+1
    }
    if (obj < 0 ) tk = 0
    return(tk)
  }
  retractionZ.X = function(X, Z) {
    return(qr.Rplus(Z+X)$q)
  }
  qr.Rplus = function(X) {
    z = qr(X)
    q = qr.Q(z,TRUE)
    r = qr.R(z,TRUE)
    q = q %*% diag(sign(diag(r)))
    r = diag(sign(diag(r))) %*% r
    return(list(q=q,r=r))
  }
  ###############################
  modelbic <- matrix(NA, nrow = length(gr)*length(models), ncol = 5)
  colnames(modelbic) <- c("G", "BIC", "df", "Loglik", "ICL")
  rownames(modelbic) <- rep(models,length(gr))
  mbi <- 1
  git <- gr[1]
  gitupper <- gr[length(gr)]
  while (git < (gitupper + 1)) {
    matinit <- 0 #NULL
    if(start=="kmeans"){
      inz <- kmeans(dat,git,nstart=kmeansinit)
      initialz <- data.matrix(unmap(inz$cluster, git)) #Kmeans initialization
    } else if(start=="random"){
      initialz <- data.matrix(unmap(sample(1:git,N,replace=TRUE),git))
    }
    if(any(!is.na(label))) initialz[!is.na(label),] <- unmap(label[!is.na(label)],git) # for classification when some labels are known

    mit<-1
    while(mit<length(models)+1)
    {
      modelname<-models[mit]
      z <- initialz
      zold <- z
      #################Parameter initialization##################
      pi_g <- NULL
      pi_gold <- pi_g
      mu <- matrix(NA, nrow = git, ncol = p)
      muold <- mu
      Sig <- list()
      invSig <- list()
      for(g in 1:git){invSig[[g]] <- diag(p)}
      m_g <- list()
      for(g in 1:git){m_g[[g]] <- sweep(dat, 2, cov.wt(as.matrix(dat), wt = (initialz[, g]), center = TRUE,
                                                        method = "ML")$center,"-")}
      Sigold <- list()
      incparll <- NULL
      betag <- rep(0.9, git)
      betagold <- betag
      kurt <- NULL
      posdef<-1
      #######################GEM loop starts#####################
      it <- 1
      stop <- 0
      deg<-1
      while (stop < 1) {
        if (it > 1) {
          ff<- fnz(label=label)
          z <- ff$z
        }
        if(all(is.finite(z))){
          if (it > 1) {
            if(substr(models[mit],4,4)=="E"){
              bg <- ms_betagE()
              if(bg$con!=0) betag <- bg$betag}
            if(substr(models[mit],4,4)=="V"){
              bg <- ms_betagV()
              if(bg$con!=0) betag <- bg$betag}
            if(any(!is.finite(betag)) || bg$con==0) {
              stop<-1
              incparll[it] <- Inf
              break
            }
          }
          if(all(is.finite(betag)) && stop<1){
            mst_x <- ms_x()
            pi_g <- mst_x$pi_g
            mu <- mst_x$mu
            if(mst_x$con==0) deg<-0
            if(deg==1){
              for(g in 1:git) m_g[[g]] <- sweep(dat,2,mu[[g]],"-") #mean centered data
              mstep_Sig <- try(model.type(modelname),silent=TRUE)
              if(mstep_Sig$con!=0) Sig <- mstep_Sig$Sig
              for(g in 1:git){
                if(any(eigen(Sig[[g]])$values<0) || mstep_Sig$con==0) {
                  posdef<-0
                  incparll[it] <- Inf
                  stop<-1
                  break
                }
              }
              if(mst_x$con==0) stop<-1
              if(posdef==1 && mst_x$con==1){
                for(g in 1:git) {
                  if(!(rcond(Sig[[g]])<sqrt(.Machine$double.eps))) invSig[[g]] <- try(solve(Sig[[g]]),
                                                                                      silent=TRUE)
                  if(!is.numeric(invSig[[g]]) || rcond(Sig[[g]])<sqrt(.Machine$double.eps)){
                    stop <- 1
                    posdef <- 0
                    incparll[it] <- Inf
                    break
                  }
                }
                #inverse of Sig---to be used in next iteration
                if(it>1) valsave()
                if (it > 1) {incparll[it] <- sum(log(rowSums(fnz(label=label)$forz)))}
                if(it > 1 && !is.finite(incparll[it])) {stop<-1; posdef<-0}
                #here posdef doesn't refer to positive definiteness not being met
                if (it > 3 && is.finite(incparll[it])) stop <- stopcheck()
                if (git == 1 & it > 2) {
                  if (abs(incparll[it] - incparll[(it - 1)]) < eps)
                    stop <- 1
                }
              }
            }
          }
        }
        if (any(!is.finite(z))) {
          stop <- 1
          incparll[it] <- Inf
          break
        }
        if(it>maxit || posdef==0) {stop<-1; break}
        it <- it + 1
      } #for EM - for it
      t<-na.omit(incparll[1:it])
      if(!all(t == cummax(t))){
        z    <- valsave()$z
        betag<- valsave()$betag
        pi_g <- valsave()$pi_g
        mu   <- valsave()$mu
        Sig  <- valsave()$Sig
        incparll[it-1] <- NA
        it <- it-1
      }
      if(all(t == cummax(t)) && all(is.finite(betag)) && all(is.finite(z)) && posdef==1 && mst_x$con==1){
        for (g in 1:git) { kurt[g] <- (p^2 * gamma(p/2/betag[g]) * gamma((p + 4)/2/betag[g]))/
                             (gamma((p + 2)/2/betag[g]))^2 - p * (p + 2) }
        para <- npar.model(modelname=modelname, p=p, G=git)
        if (incparll[(it - 1)] == Inf || ncol(z) != git) bic <- -Inf else
          bic <- 2 * incparll[(it - 1)] - log(N) * para
        #Here, higher BIC is better
        forICL <- function(g) {
          sum(log(z[max.col(z)==g, g]))
        }
        icl <- bic + sum(sapply(unique(max.col(z)), forICL))
        results[[paste("G", git, modelname, sep = "_")]] <- list(bic = bic, icl = icl,
                                                                          z = z, para = para,
                                                                          ll = incparll[(it-1)],
                                                                          pi_g = pi_g, mu = mu,
                                                                          Sig = Sig, beta=betag,
                                                                          kurt=kurt, var=varfn(Sig))
      } else{
        results[[paste("G", git, modelname, sep = "_")]] <- list(bic = NA, icl = NA,  z = z,
                                                                          para = NA, ll = incparll[2:(it-1)],
                                                                          pi_g = NA, mu = NA, Sig = NA,
                                                                          beta=NA, kurt=NA, var=NA)
      }
      if(all(t == cummax(t)) && all(is.finite(betag)) && all(is.finite(z)) && posdef==1  && mst_x$con==1){
        if(all(as.numeric(betag)<200)){
          modelbic[mbi, 1] <- git
          modelbic[mbi, 2] <- bic
          modelbic[mbi, 3] <- para
          modelbic[mbi, 4] <- incparll[(it - 1)]
          modelbic[mbi, 5] <- icl
        }
      } else{
        modelbic[mbi, 1] <- git;     modelbic[mbi, 4] <- NA; }
      if(verbose==TRUE) cat(modelname, "done for G =", git, "\n")
      mbi <- mbi + 1
      mit <- mit + 1
    }
    git <- git + 1
  } #for git
  bestbic <- which.max(modelbic[,2])
  besticl <- which.max(modelbic[,5])
  #msc stands for model selection criterion
  timetaken <- proc.time() - ptm
  val <- list(call = match.call(), time = timetaken,
              modelnames = models,msc=modelbic,bicclassification=map(results[[bestbic]]$z),
              iclclassification=map(results[[besticl]]$z),
              bicselection=results[[bestbic]],iclselection=results[[besticl]])
  class(val) <- "pemix"
  return(invisible(val))
}
