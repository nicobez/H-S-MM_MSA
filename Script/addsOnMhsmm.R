# # # # # # # # # # # # # # # # # # # # # # # #
#
# AR1 : m*(1-mu)=eta where mu=parameters[,,,3,] is the coef of auto-correlation
# AR0 : mu, eta parameters are not used
# HMM : sojourn times are geometrical ==> rgeom(1,prob=1/m.soj.1)
# HSMM : sojourn time are shifted-Poisson (time ~ shift + Poisson(m-shift) ; variance = m-shift) 
#                      or shifted-NegBinom (time ~ shift + nbinom(size,prob) where size*(1-prob)/prob = m-shift ; variance = size*(1-prob)/prob?)
#
# state 1 = NON Fishing ; state 2 = Fishing
#
# m.V.1 : means of state dependent variables in state 1
# s2.V.1 : variances of state dependent variables in state 1 ; independence ==> vector and not a matrix 
# mu.V.1 : coefficients for AR1
# m.soj : means of sojourn time in the different states ==> vector with length the number of states
#         units are in number of time steps
#
# The following codes work for either of the two sojour time pmfs 
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
    environment(mstep) <- .GlobalEnv # pour que la fonction cov.wt modifi?e soit lue
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

# -------------------------------------------------------
vessel.model.builder <- function(method="HMM",type="poisson",AR="AR0", vessel=1, degrad=1, parameters){ 
  tmp <- parameters[vessel,degrad,,,]
  raw.model.builder(method=method,type=type, AR=AR, theta=tmp)
}

# -------------------------------------------------------
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
# Emission density function for a multivariate normal emission distribution assuming an AR process.
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

# -------------------------------------------------------
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

# -------------------------------------------------------
simulate.hmmspec <- function (object, nsim, seed = NULL, rand.emission = NULL, ...) 
# modification of mhsmm:::simulate.hmmspec() to add the case AR1 (lines 10-15)
{
  if (!is.null(seed)) 
    set.seed(seed)
  if (is.null(rand.emission) & is.null(object$rand.emission)) 
    stop("rand.emission not specified")
  if (is.null(rand.emission)) 
    rand.emission = object$rand.emission
  if (length(nsim) == 1) {
    s1 = sim.mc(object$init, object$transition, nsim)
    if(isTRUE(all.equal(rmvnormAR.hsmm,rand.emission))){
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


#####################################################################################
# Functions replacing generic functions of mhsmm package 
# To call the native functions, one must use absolute calls like mhsmm:::hsmmspec()
# A simple call like hsmmspec() points on the local added modified version. 
#####################################################################################

# -------------------------------------------------------
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

# -------------------------------------------------------
simulate.hsmmspec <- function (object, nsim, seed = NULL, rand.emission = NULL, ...) 
{
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
    if(isTRUE(all.equal(rmvnormAR.hsmm,rand.emission))){
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

# -------------------------------------------------------
cov.wt <-function (x, wt = rep(1/nrow(x), nrow(x)), cor = FALSE, center = TRUE, 
                   method = "ML") 
{
  # Modification of native stats:::cov.wt() to force it to provide diagonal for multivariable cases, i.e. independance between variables
  # Remove the choice between "unbiased" and "ML" computations ; method='unbiased' gets a degenerated denominator if 1-sum(wt^2) = 0
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

# -------------------------------------------------------
hsmmfit <- function (x, model, mstep = NULL, M = NA, maxit = 100, tol= 1e-04,lock.transition = FALSE, 
          lock.d = FALSE, graphical = FALSE) 
{
  sojourn.distribution = model$sojourn$type
  #tol = 1e-04
  ksmooth.thresh = 1e-20
  shiftthresh = 1e-20
  J = nrow(model$transition)
  model$J = J
  if (is.null(mstep)) 
    if (is.null(model$mstep)) 
      stop("mstep not specified")
  else mstep = model$mstep
  mhsmm:::.check.hsmmspec(model)
  f = model$dens.emission
  if (mode(x) == "numeric" | mode(x) == "integer") {
    warning("x is a primitive vector.  Assuming single sequence.")
    NN = NROW(x)
  }
  else {
    NN = x$N
    x = x$x
  }
  if (is.na(M)) 
    M = max(NN)
  if (length(model$init) != J) 
    stop("length(model$init)!=J")
  if (NROW(x) != sum(NN)) 
    stop("NROW(x)!=sum(NN)")
  model <- mhsmm:::.build_d(model, M)
  new.model = model
  ll = rep(NA, maxit)
  rm(model)
  for (it in 1:maxit) {
    if (graphical) 
      plot.hsmm(list(model = new.model, J = J))
    p = sapply(1:J, function(state) f(x, state, new.model))
    if (any(is.na(p) | p == Inf)) 
      stop("NAs detected in b(x), check your supplied density function")
    if (any(apply(p, 1, max) == 0)) 
      stop("Some values have 0 pdf for all states!  Check your model parameters")
    estep_variables = .C("backward", transition = as.double(new.model$transition), 
                         init = as.double(new.model$init), p = as.double(p), 
                         d = as.double(new.model$d), D = as.double(new.model$D), 
                         timelength = as.integer(NN), J = as.integer(J), 
                         M = as.integer(rep(M, J)), L1 = double(NROW(x) * 
                                                                  J), N = double(NROW(x)), eta = double(M * J), 
                         F1 = double(J * NROW(x)), si = double(J * NROW(x)), 
                         gamma = double(J * NROW(x)), nsequences = as.integer(length(NN)), 
                         totallength = NROW(x), G = double(J * NROW(x)), 
                         PACKAGE = "mhsmm")
    if (any(is.nan(estep_variables$gamma))) {
      warning("NaNs detected in gamma.  Exiting...")
      return(estep_variables)
    }
    if (any(estep_variables$gamma < 0)) 
      estep_variables$gamma = zapsmall(estep_variables$gamma)
    if (any(estep_variables$eta < 0)) 
      estep_variables$eta = zapsmall(estep_variables$eta)
    if (any(estep_variables$N < 0)) 
      estep_variables$N = zapsmall(estep_variables$N)
    old.model = new.model
    state_wt <- matrix(estep_variables$gamma, ncol = J)
    if (any(colSums(state_wt) == 0)) 
      stop("Error: at least one state has an expected number of occurences equal to 0.\n This may be caused by bad starting parameters are insufficent sample size")
    new.model$parms.emission = mstep(x, state_wt)
    if (lock.d) {
      new.model$d = old.model$d
      new.model$D = old.model$D
    }
    else {
      if (sojourn.distribution == "nonparametric") {
        new.model$d = apply(matrix(estep_variables$eta, 
                                   ncol = J), 2, function(x) x/sum(x))
        new.model$sojourn$d <- new.model$d
      }
      else if (sojourn.distribution == "ksmoothed-nonparametric") {
        new.model$d = apply(matrix(estep_variables$eta + 
                                     1e-100, ncol = J), 2, function(x) x/sum(x))
        for (i in 1:J) {
          new.model$d[, i] = approx(density(which(new.model$d[, 
                                                              i] > ksmooth.thresh), weights = new.model$d[which(new.model$d[, 
                                                                                                                            i] > ksmooth.thresh), i], from = 1, n = M), 
                                    xout = 1:M)$y
          new.model$d[is.na(new.model$d[, i]), i] = 0
          new.model$d[, i] = (new.model$d[, i] + 1e-300)/sum(new.model$d[, 
                                                                         i])
        }
        new.model$sojourn$d <- new.model$d
      }
      else if (sojourn.distribution == "poisson") {
        new.model$d = apply(matrix(estep_variables$eta, 
                                   ncol = J), 2, function(x) x/sum(x))
        new.model$sojourn$lambda = numeric(J)
        new.model$sojourn$shift = numeric(J)
        for (i in 1:J) {
          eta = new.model$d[, i]
          maxshift = match(TRUE, eta > shiftthresh)
          Mtmp = tail(which(eta > shiftthresh), 1)
          new.model$sojourn$shift[i] = which.max(sapply(1:maxshift, 
                                                        function(shift) .dpois.hsmm.sojourn(x = maxshift:Mtmp, 
                                                                                            lambda = ((maxshift:Mtmp) - shift) %*% 
                                                                                              eta[maxshift:Mtmp], shift = shift, log = TRUE) %*% 
                                                          eta[maxshift:Mtmp]))
          new.model$sojourn$lambda[i] = ((new.model$sojourn$shift[i]:Mtmp) - 
                                           new.model$sojourn$shift[i]) %*% eta[new.model$sojourn$shift[i]:Mtmp]
          new.model$d[, i] = .dpois.hsmm.sojourn(1:M, 
                                                 new.model$sojourn$lambda[i], new.model$sojourn$shift[i])
        }
      }
      else if (sojourn.distribution == "nbinom") {
        new.model$d = matrix(nrow = M, ncol = J)
        new.model$sojourn$size = numeric(J)
        new.model$sojourn$shift = integer(J)
        new.model$sojourn$mu = numeric(J)
        new.model$sojourn$prob = numeric(J)
        eta = apply(matrix(estep_variables$eta, ncol = J), 2, function(x) x/sum(x))
        for (i in 1:J) {
          tmp = fitnbinom(eta[, i])
          new.model$sojourn$shift[i] = tmp[1]
          new.model$sojourn$size[i] = tmp[2]
          new.model$sojourn$mu[i] = tmp[3]
          new.model$sojourn$prob[i] = tmp[4]
          new.model$d[, i] = mhsmm:::.dnbinom.hsmm.sojourn(1:M, 
                                                   new.model$sojourn$size[i], new.model$sojourn$prob[i], 
                                                   new.model$sojourn$shift[i])
        }
      }
      else if (sojourn.distribution == "gamma") {
        new.model$d = matrix(estep_variables$eta, ncol = J)
        new.model$sojourn$shape = numeric(J)
        new.model$sojourn$scale = numeric(J)
        for (i in 1:J) {
          tmp = gammafit(1:M, wt = new.model$d[, i])
          new.model$sojourn$shape[i] = tmp$shape
          new.model$sojourn$scale[i] = tmp$scale
          new.model$d[, i] = dgamma(1:M, shape = tmp$shape, 
                                    scale = tmp$scale)
        }
      }
      else if (sojourn.distribution == "logarithmic") {
        new.model$d = apply(matrix(estep_variables$eta + 
                                     1e-100, ncol = J), 2, function(x) x/sum(x))
        new.model$sojourn$shape = numeric(J)
        for (i in 1:J) {
          new.model$sojourn$shape[i] = .logdistrfit(wt = new.model$d[, 
                                                                     i])
          new.model$d[, i] = .dlog(1:M, new.model$sojourn$shape[i])
        }
      }
      else if (sojourn.distribution == "lnorm") {
        eta = matrix(estep_variables$eta, ncol = J)
        new.model$d = matrix(nrow = M, ncol = J)
        new.model$sojourn$meanlog = numeric(J)
        new.model$sojourn$s.dlog = numeric(J)
        for (i in 1:J) {
          new.model$sojourn$meanlog[i] = weighted.mean(log(1:M), 
                                                       eta[, i])
          new.model$sojourn$s.dlog[i] = sqrt(cov.wt(data.frame(log(1:M)), 
                                                    eta[, i])$cov)
          new.model$d[, i] = dlnorm(1:M, new.model$sojourn$meanlog[i], 
                                    new.model$sojourn$s.dlog[i])
          new.model$d[, i] = new.model$d[, i]/sum(new.model$d[, 
                                                              i])
        }
      }
      else stop("Invalid sojourn distribution")
      new.model$D = apply(new.model$d, 2, function(x) rev(cumsum(rev(x))))
    }
    if (lock.transition) {
      new.model$init = old.model$init
      new.model$transition = old.model$transition
    }
    else {
      new.model$init = estep_variables$init
      new.model$init[new.model$init < 0] = 0
      new.model$transition = matrix(estep_variables$transition, 
                                    ncol = J)
      new.model$transition[new.model$transition < 0] = 0
    }
    ll[it] = sum(log(estep_variables$N))
    new.model$J = J
    if (it > 2) 
      if (abs(ll[it] - ll[it - 1]) < tol) 
        (break)()
  }
  class(new.model) <- "hsmmspec"
  ret = list(loglik = ll[!is.na(ll)], model = new.model, estep_variables = estep_variables, 
             M = M, J = J, NN = NN, f = f, mstep = mstep, yhat = apply(matrix(estep_variables$gamma, 
                                                                              ncol = J), 1, which.max))
  class(ret) <- "hsmm"
  ret
}

# ---------------------------------
fitnbinom <- function (eta) 
{
  # same as mhsmm:::.fitnbinom except that cov.wt() is now the local one
  shiftthresh = 1e-20
  maxshift = match(TRUE, eta > shiftthresh)
  Mtmp = tail(which(eta > shiftthresh), 1)
  fun1 <- function(shift) {
    m <- weighted.mean((maxshift:Mtmp) - shift, eta[maxshift:Mtmp])
    v <- as.numeric(stats:::cov.wt(data.frame((maxshift:Mtmp) - shift), 
                           wt = eta[maxshift:Mtmp])$cov)
    size <- if (v > m) 
      m^2/(v - m)
    else 100
    densfun <- function(par) sum(dnbinom((maxshift:Mtmp) - 
                                           shift, size = par[1], mu = par[2], log = TRUE) * 
                                   eta[maxshift:Mtmp])
    optim(c(size, m), densfun, control = list(fnscale = -1))$value
  }
  shift = which.max(sapply(1:maxshift, fun1))
  m <- weighted.mean((maxshift:Mtmp) - shift, eta[maxshift:Mtmp])
  v <- as.numeric(stats:::cov.wt(data.frame((maxshift:Mtmp) - shift), 
                         wt = eta[maxshift:Mtmp])$cov)
  size <- if (v > m) 
    m^2/(v - m)
  else 100
  densfun <- function(par) sum(dnbinom((maxshift:Mtmp) - shift, 
                                       size = par[1], mu = par[2], log = TRUE) * eta[maxshift:Mtmp])
  tmp = optim(c(size, m), densfun, control = list(fnscale = -1))$par
  c(shift = shift, size = tmp[1], mu = tmp[2], prob = tmp[1]/(sum(tmp)))
}

# ----------
f.np2mv <- function(n,p){
  # transforms the (n,p)=(size,prob) values of a NegBinom into the corresponding (mean,variance)
  m <- n*(1-p)/p
  v <- n*(1-p)/p^2
  c(m,v)
}

f.mv2np <- function(m,v){
  # transforms the (mean,variance) values of a NegBinom into the corresponding (n,p)=(size,prob) 
  n <- m^2/(v-m)
  p <- m/v
  c(n,p)
}






