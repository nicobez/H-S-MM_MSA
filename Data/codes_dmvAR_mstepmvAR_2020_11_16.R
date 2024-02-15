# # # # # # # # # # # # # # # # # # # # # # # #
# NOVEMBER 2020
#
# AR1 : m*(1-mu)=eta
# AR0 : mu, eta parameters are not used
# HMM : sojourn times are geometrical ==> rgeom(1,prob=1/m.soj.1)
# HSMM : sojourn time are shifted-Poisson (time ~ shift + Poisson(m-shift) ; variance = m-shift) 
#                      or shifted-NegBinom (time ~ shift + nbinom(size,prob) where size*(1-prob)/prob = m-shift ; variance = size*(1-prob)/prob²)
#
# state 1 = NON Fishing ; state 2 = Fishing
#
# m.V.1 : moyennes des variables dans l'état 1
# s2.V.1 : variances ; indépendance donc pas de covariances ==> vecteur et non matrice 
# mu.V.1 : coefficients des AR1
# m.soj : moyennes de temps de séjour dans les états ==> vecteur de dimension le nbr d'états
#         les unités sont en nombre de pas de temps
#
# Need f.np2mv() and f.mv2np()
# ==> codes_working_functions_2020_11_16.R
#
# The following codes work for either of the the two sojour time pdfs 
#    parametersP <- parameters[,,c(1,2,3),,]
#    parametersNB <- parameters[,,c(1,2,4),,]
#
# # # # # # # # # # # # # # # # # # # # # # # #


# Building models -------------------------------------------------------
raw.model.builder <- function(method="HMM",type="poisson",AR="AR0",theta){ 
  # theta = parameters[i.v,i.d,,,] set of relevant parameters
  tmp <- array(NA,dim=c(3,5,2))
  tmp <- theta[,,]
  theta <- tmp
  nstate <- 2 ; nvar <- 2
  # Defining initial state
  s0 = c(1,rep(0,nstate-1)) 
  #
  init0=list(m.V=list(m.V.1=theta[1:2,1,1],m.V.2=theta[1:2,1,2]),
             var.V=list(s2.V.1=diag(theta[1:2,2,1]),s2.V.2=diag(theta[1:2,2,2])),
             mu.V=list(mu.V.1=theta[1:2,3,1],mu.V.2=theta[1:2,3,2]),
             eta.V=list(eta.V.1=theta[1:2,4,1],eta.V.2=theta[1:2,4,2]),
             shift.soj=c(theta[3,5,]),
             m.soj=c(theta[3,1,]),
             var.soj=c(theta[3,2,]))
  if (AR=="AR1") {
    parms.emission = list(eta=init0$eta.V,mu=init0$mu.V,sigma=init0$var.V)
    dens.emission  = dmvnormAR.hsmm 
    rand.emission = rmvnormAR.hsmm
    mstep = mstep.mvAR
    x0 <- rep(0,nvar)
  } else {
    parms.emission = list(mu=init0$m.V,sigma=init0$var.V)
    dens.emission  = dmvnorm.hsmm 
    rand.emission = rmvnorm.hsmm
    mstep = mstep.mvnorm
    environment(mstep) <- .GlobalEnv # pour que la fonction cov.wt modifiée soit lue
  }
  transition <- matrix(0,nrow=nstate,ncol=nstate)
  if(method=="HMM"){
    diag(transition) <- sapply(1:nstate,function(i) max(1- 1/init0$m.soj[i],0))
    tmp <- sapply(1:nstate,function(i) transition[i,-i] <<- rdirichlet(1,rep(10,nstate-1))*(1-transition[i,i]))
    true_model <- hmmspec(s0,transition,parms.emission,dens.emission,rand.emission,mstep)
  } else {
    tmp <- sapply(1:nstate,function(i) transition[i,-i] <<- rdirichlet(1,rep(10,nstate-1)))
    if(type=="poisson") sojourn = list(shift  = init0$shift.soj,
                                       lambda = init0$m.soj-init0$shift.soj,
                                       type="poisson")
    if(type=="nbinom") {
      tmp1 <- f.mv2np(init0$m.soj[1]-init0$shift.soj[1],init0$var.soj[1]) 
      tmp2 <- f.mv2np(init0$m.soj[2]-init0$shift.soj[2],init0$var.soj[2])
      sojourn = list(shift = init0$shift.soj,
                     size = c(tmp1[1],tmp2[1]),
                     prob = c(tmp1[2],tmp2[2]),
                     type="nbinom")
    }
    true_model <- hsmmspec(s0,transition,parms.emission,sojourn,dens.emission,rand.emission,mstep)
  }
  true_model
}

vessel.model.builder <- function(method="HMM",type="poisson",AR="AR0", vessel=1, degrad=1, parameters){ 
  tmp <- parameters[vessel,degrad,,,]
  raw.model.builder(method=method,type=type, AR=AR, theta=tmp)
}

label.switching <- function(objet,AR="AR0"){
  # Rebuilt outputs if there has been switche between states
  # State 1 = NON Fishing ; state 2 = FISHING
  # There has been switch when mean(Vp.1) < mean(Vp.2)
  # Input : objet is output from h(s)mmfit
  res <- objet
  tmp <- objet$model
  if(AR=="AR0"){
    flag.switch <- tmp$parms.emission$mu[[1]][1] < tmp$parms.emission$mu[[2]][1]
    if(flag.switch){
      cat("Label switching ! \n")
      res$model$parms.emission$mu[[1]] <- tmp$parms.emission$mu[[2]]
      res$model$parms.emission$mu[[2]] <- tmp$parms.emission$mu[[1]]
      res$model$parms.emission$sigma[[1]] <- tmp$parms.emission$sigma[[2]]
      res$model$parms.emission$sigma[[2]] <- tmp$parms.emission$sigma[[1]]  
      res$model$transition[1,1] <- tmp$transition[2,2]  
      res$model$transition[1,2] <- tmp$transition[2,1]  
      res$model$transition[2,1] <- tmp$transition[1,2]  
      res$model$transition[2,2] <- tmp$transition[1,1]  
      res$yhat <- res$yhat+1 ; res$yhat[res$yhat==3] <- 1
      if(class(tmp)=="hsmmspec"){
        if(tmp$sojourn$type=="poisson"){
          res$model$sojourn$shift[1] <- tmp$sojourn$shift[2]
          res$model$sojourn$shift[2] <- tmp$sojourn$shift[1]
          res$model$sojourn$lambda[1] <- tmp$sojourn$lambda[2]
          res$model$sojourn$lambda[2] <- tmp$sojourn$lambda[1]
        }
        if(tmp$sojourn$type=="nbinom"){
          res$model$sojourn$shift[1] <- tmp$sojourn$shift[2]
          res$model$sojourn$shift[2] <- tmp$sojourn$shift[1]
          res$model$sojourn$size[1] <- tmp$sojourn$size[2]
          res$model$sojourn$size[2] <- tmp$sojourn$size[1] 
          res$model$sojourn$prob[1] <- tmp$sojourn$prob[2]
          res$model$sojourn$prob[2] <- tmp$sojourn$prob[1]  
        }
      }
    }
  } else {
    m1 <- tmp$parms.emission$eta[[1]][1]/(1-tmp$parms.emission$mu[[1]][1])
    m2 <- tmp$parms.emission$eta[[2]][1]/(1-tmp$parms.emission$mu[[2]][1])
    flag.switch <-  m1 < m2
    if(flag.switch){
      cat("Label switching ! \n")
      res$model$parms.emission$eta[[1]] <- tmp$parms.emission$eta[[2]]
      res$model$parms.emission$eta[[2]] <- tmp$parms.emission$eta[[1]]
      res$model$parms.emission$mu[[1]] <- tmp$parms.emission$mu[[2]]
      res$model$parms.emission$mu[[2]] <- tmp$parms.emission$mu[[1]]
      res$model$parms.emission$sigma[[1]] <- tmp$parms.emission$sigma[[2]]
      res$model$parms.emission$sigma[[2]] <- tmp$parms.emission$sigma[[1]]  
      res$model$transition[1,1] <- tmp$transition[2,2]  
      res$model$transition[1,2] <- tmp$transition[2,1]  
      res$model$transition[2,1] <- tmp$transition[1,2]  
      res$model$transition[2,2] <- tmp$transition[1,1]  
      res$yhat <- res$yhat+1 ; res$yhat[res$yhat==3] <- 1
      if(class(tmp)=="hsmmspec"){
        if(tmp$sojourn$type=="poisson"){
          res$model$sojourn$shift[1] <- tmp$sojourn$shift[2]
          res$model$sojourn$shift[2] <- tmp$sojourn$shift[1]
          res$model$sojourn$lambda[1] <- tmp$sojourn$lambda[2]
          res$model$sojourn$lambda[2] <- tmp$sojourn$lambda[1]
        }
        if(tmp$sojourn$type=="nbinom"){
          res$model$sojourn$shift[1] <- tmp$sojourn$shift[2]
          res$model$sojourn$shift[2] <- tmp$sojourn$shift[1]
          res$model$sojourn$size[1] <- tmp$sojourn$size[2]
          res$model$sojourn$size[2] <- tmp$sojourn$size[1] 
          res$model$sojourn$prob[1] <- tmp$sojourn$prob[2]
          res$model$sojourn$prob[2] <- tmp$sojourn$prob[1]  
        }
      }
    }
  }
  res
}



# Additional functions for AR --------------------------------------------------------

# Emission density function for a multivariate normal emission distribution
# assuming an AR process.
# Their postfix *.hsmm does not mean they are specific to HSMM.
# They can be used for HMM aswell.

dmvnormAR.hsmm <- function(x, j, model){
  N <- NROW(x);
  if(N==1){
    break("Error, Number of observations should be equal at least to 2")
  }
  if(NCOL(x)==1){
    dens.fun <- dnorm
    model$parms.emission$sigma <- lapply(model$parms.emission$sigma,
                                         function(sig) as.numeric(sqrt(sig)))
  }
  else{
    dens.fun <- dmvnorm
  }
  probs_m1 <- sapply(2:N,function(n){
    eta <- model$parms.emission$eta[[j]];
    mu <- model$parms.emission$mu[[j]];
    sigma <- model$parms.emission$sigma[[j]];
    ans <- dens.fun(x[n,],eta+mu*x[n-1,],sigma);# specification of AR
    ans[is.na(ans)] <- 1;
    return(ans)
  })
  probs <- c(1,probs_m1)
  return(probs)
}

#Performs re-estimation (M-Step) for a multivariate normal emission distribution
#assuming an AR process
mstep.mvAR <- function(x,wt){
  idx <- apply(is.na(x), 1, any)
  x <- x[!idx, , drop = FALSE]
  nobs <- nrow(x)
  wt <- wt[!idx, , drop = FALSE]
  nstate <- ncol(wt)
  new.mu <- lapply(1:nstate,function(i){
    gams <- wt[-1,i]
    vts <- matrix(x[-1,],ncol=ncol(x)); vtsm1 <- matrix(x[-nobs,],ncol=ncol(x))
    num <- sum(gams)*colSums(gams*vts*vtsm1)-colSums(gams*vtsm1)*colSums(gams*vts)
    denom <- sum(gams)*colSums(gams*vtsm1*vtsm1)-colSums(gams*vtsm1)^2
    num/denom
  })
  new.eta <- lapply(1:nstate,function(i){
    gams <- wt[-1,i]
    vts <- matrix(x[-1,],ncol=ncol(x)); vtsm1 <- matrix(x[-nobs,],ncol=ncol(x))
    mus.mat <- matrix(new.mu[[i]],ncol=ncol(x),nrow=nobs-1,byrow=T)
    num <- colSums(gams*(vts-mus.mat*vtsm1))
    denom <- sum(gams)
    num/denom
  })
  new.sigma <- lapply(1:nstate,function(i){
    gams <- wt[-1,i]
    vts <- matrix(x[-1,],ncol=ncol(x)); vtsm1 <- matrix(x[-nobs,],ncol=ncol(x))
    mus.mat <- matrix(new.mu[[i]],ncol=ncol(x),nrow=nobs-1,byrow=T)
    etas.mat <- matrix(new.eta[[i]],ncol=ncol(x),nrow=nobs-1,byrow=T)
    num <- colSums(gams*(vts-etas.mat-mus.mat*vtsm1)^2)
    denom <- sum(gams)
    num/denom
  })
  return(list(eta=new.eta,mu=new.mu,
              sigma = lapply(new.sigma,function(sig) diag(sig,ncol(x)))))
}



rmvnormAR.hsmm <- function(seq,model,x0=NULL){
  dim.x <- length(model$parms.emission$eta[[1]])
  etas <- model$parms.emission$eta;
  mus <- model$parms.emission$mu;
  if(dim.x==1){
    rand.fun <- rnorm 
    sigmas <- lapply(model$parms.emission$sigma,
                     function(sig) as.numeric(sqrt(sig)));
  }
  else{
    rand.fun <- rmvnorm
    sigmas <- model$parms.emission$sigma
  }
  if(is.null(x0)){
    x0 <- rand.fun(1,etas[[seq[1]]]/(1-mus[[seq[1]]]),sigmas[[seq[1]]]/(1-mus[[seq[1]]]^2))
  }
  res <- matrix(x0,nrow=length(seq),ncol=length(x0),byrow=T)
  for(i in 2:length(seq)){
    res[i,] <- rand.fun(1,
                       etas[[seq[i]]] + res[i-1,]*mus[[seq[i]]],
                       sigmas[[seq[i]]])
  }
  res
}

simulate.hmmspec.mod <- function (object, nsim, seed = NULL, rand.emission = NULL, ...) 
{
  if (!is.null(seed)) 
    set.seed(seed)
  if (is.null(rand.emission) & is.null(object$rand.emission)) 
    stop("rand.emission not specified")
  if (is.null(rand.emission)) 
    rand.emission = object$rand.emission
  if (length(nsim) == 1) {
    s1 = sim.mc(object$init, object$transition, nsim)
    if(all.equal(rmvnormAR.hsmm,rand.emission)){
      z <- list(...)
      x <- rmvnormAR.hsmm(s1,object,z$x0)
      ret = list(s = s1, x = x, N = NROW(x))
    }
    else{    
      x = sapply(s1, rand.emission, model = object)
      if (NCOL(x) > 1) 
        ret = list(s = s1, x = t(x), N = nsim)
      else ret = list(s = s1, x = x, N = nsim)
    }
    class(ret) <- "hsmm.data"
    ret
  }
  else mhsmm:::.sim.mhmm(object, nsim, rand.emission)
}

simulate.hsmmspec.mod <- function (object, nsim, seed = NULL, rand.emission = NULL, ...){
  # Modification of native mhsmm:::simulate.hsmmspec() to :
  #  - allow uploading rand.emission from the input object (line 6-7)
  #  - add shifted-nbinom case
  #  - account for AR1
  # Use of explicit call to mhsmm functions when required (mshmm:::)
  right.truncate = left.truncate = 0
  if (!is.null(seed)) 
    set.seed(seed)
  if (is.null(rand.emission) & is.null(object$rand.emission)) 
    stop("rand.emission not specified")
  if (is.null(rand.emission)) 
    rand.emission <- object$rand.emission
  if (length(nsim) == 1) {
    s0 = sim.mc(object$init, object$transition, nsim)
    if (object$sojourn$type == "poisson") {
      fn <- function(ii) mhsmm:::.rpois.hsmm.sojourn(1, object$sojourn$lambda[ii], 
                                                     object$sojourn$shift[ii])
    }
    else if (object$sojourn$type == "nbinom") {
      fn <- function(ii) object$sojourn$shift[ii]+rnbinom(1, size =  object$sojourn$size[ii], 
                                                          prob = object$sojourn$prob[ii])
    }
    else if (object$sojourn$type == "gamma") {
      fn <- function(ii) rgamma(1, shape = object$sojourn$shape[ii], 
                                scale = object$sojourn$scale[ii])
    }
    else if (object$sojourn$type == "logarithmic") {
      fn <- function(ii) mhsmm:::.rlog(1, object$sojourn$p[ii])
    }
    else if (object$sojourn$type == "nonparametric") {
      fn <- function(ii) sample(1:nrow(object$sojourn$d), 
                                1, prob = object$sojourn$d[, ii])
    }
    else stop("Unknown sojourn distribution")
    u = as.integer(round(sapply(s0, fn)))
    s1 = rep(s0, u)[(left.truncate + 1):(sum(u) - right.truncate)]
    if(all.equal(rmvnormAR.hsmm,rand.emission)){
      z <- list(...)
      x <- rmvnormAR.hsmm(s1,object,z$x0)
      ret = list(s = s1, x = x, N = NROW(x))
    }
    else{
      x = sapply(s1, function(i) rand.emission(i, object))
      if (NCOL(x) > 1) 
        ret = list(s = s1, x = t(x), N = NCOL(x))
      else ret = list(s = s1, x = x, N = NROW(x))
    }
    class(ret) <- "hsmm.data"
    ret
  }
  else {
    mhsmm:::.sim.mhsmm(nsim, object, object$sojourn$type, object$rand.emission)
  }
}


# Functions replacing generic functions of mhsmm package -------------------

hsmmspec <- function (init, transition, parms.emission, sojourn, dens.emission, 
                      rand.emission = NULL, mstep = NULL) 
{
  # Modification of mhsmm:::hsmmspec() to not-stop when sojourn time is nbinom
  # Use of explicit call to mhsmm functions when required (mshmm:::)
  if (is.null(dens.emission)) 
    stop("dens.emission not specified")
  if (length(init) != NROW(transition)) 
    stop("length(init)!=NROW(transition)")
  if (NROW(transition) != NCOL(transition)) 
    stop("NROW(transition)!=NCOL(transition)")
  if (is.null(sojourn$type)) 
    stop("Sojourn distribution type not specified.")
  if (all(sojourn$type != c("nonparametric", "gamma", "poisson","nbinom"))) 
    stop(paste("Invalid sojourn type specified (", sojourn$type, 
               ")"))
  ans = list(J = length(init), init = init, transition = transition, 
             parms.emission = parms.emission, sojourn = sojourn, 
             rand.emission = rand.emission, dens.emission = dens.emission, 
             mstep = mstep)
  class(ans) <- "hsmmspec"
  mhsmm:::.check.hsmmspec(ans)
  ans
}

simulate.hsmmspec <- function (object, nsim, seed = NULL, rand.emission = NULL, ...) 
{
  # Modification of native mhsmm:::simulate.hsmmspec() to :
  #  - allow uploading rand.emission from the input object (line 6-7)
  #  - add shifted-nbinom case
  # Use of explicit call to mhsmm functions when required (mshmm:::)
  right.truncate = left.truncate = 0
  if (!is.null(seed)) 
    set.seed(seed)
  if (is.null(rand.emission) & is.null(object$rand.emission)) 
    stop("rand.emission not specified")
  if (is.null(rand.emission)) 
    rand.emission <- object$rand.emission
  if (length(nsim) == 1) {
    s0 = sim.mc(object$init, object$transition, nsim)
    if (object$sojourn$type == "poisson") {
      fn <- function(ii) mhsmm:::.rpois.hsmm.sojourn(1, object$sojourn$lambda[ii], 
                                                     object$sojourn$shift[ii])
    }
    else if (object$sojourn$type == "gamma") {
      fn <- function(ii) rgamma(1, shape = object$sojourn$shape[ii], 
                                scale = object$sojourn$scale[ii])
    }
    else if (object$sojourn$type == "nbinom") {
      #fn <- function(ii) object$sojourn$shift[ii]+rnbinom(1, size =  object$sojourn$size[ii], 
      #                                                    prob = object$sojourn$prob[ii])
      # Modif 10/11/2023
	fn <- function(ii) mhsmm:::.rnbinom.hsmm.sojourn(1, size =  object$sojourn$size[ii], 
                                                          prob = object$sojourn$prob[ii],
							 shift = object$sojourn$shift[ii])
    }
    else if (object$sojourn$type == "logarithmic") {
      fn <- function(ii) mhsmm:::.rlog(1, object$sojourn$p[ii])
    }
    else if (object$sojourn$type == "nonparametric") {
      fn <- function(ii) sample(1:nrow(object$sojourn$d), 
                                1, prob = object$sojourn$d[, ii])
    }
    else stop("Unknown sojourn distribution")
    u = as.integer(round(sapply(s0, fn)))
    s1 = rep(s0, u)[(left.truncate + 1):(sum(u) - right.truncate)]
    x = sapply(s1, function(i) rand.emission(i, object))
    if (NCOL(x) > 1) 
      ret = list(s = s1, x = t(x), N = NCOL(x))
    else ret = list(s = s1, x = x, N = NROW(x))
    class(ret) <- "hsmm.data"
    ret
  }
  else {
    mhsmm:::.sim.mhsmm(nsim, object, object$sojourn$type, object$rand.emission)
  }
}

cov.wt <-function (x, wt = rep(1/nrow(x), nrow(x)), cor = FALSE, center = TRUE, 
                   method = c("unbiased", "ML")) 
{
  # Mofidication of native mshmm:::.cov.wt() to force it to be diagonal, i.e. independance between variables
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  else if (!is.matrix(x)) 
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x))) 
    stop("'x' must contain finite values only")
  n <- nrow(x)
  if (with.wt <- !missing(wt)) {
    if (length(wt) != n) 
      stop("length of 'wt' must equal the number of rows in 'x'")
    if (any(wt < 0) || (s <- sum(wt)) == 0) 
      stop("weights must be non-negative and not all zero")
    wt <- wt/s
  }
  if (is.logical(center)) {
    center <- if (center) 
      colSums(wt * x)
    else 0
  }
  else {
    if (length(center) != ncol(x)) 
      stop("length of 'center' must equal the number of columns in 'x'")
  }
  x <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)
  cov <- switch(match.arg(method), unbiased = crossprod(x)/(1 - 
                                                              sum(wt^2)), ML = crossprod(x))
  y <- list(cov = diag(diag(cov),nrow=nrow(cov)), center = center, n.obs = n)
  if (with.wt) 
    y$wt <- wt
  if (cor) {
    Is <- 1/sqrt(diag(cov))
    R <- cov
    R[] <- Is * cov * rep(Is, each = nrow(cov))
    y$cor <- R
  }
  y
}


