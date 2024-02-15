##################
# Code for MSA experiments (estimation model = simulation model)
##################


##################
# Working variables 
##################
listVessel <- c(924675,"007")               
listDegrad <- rbind(c(1,2,4,8),c(1,3,6,12))
listMethod <- c("HMM","HSMM")
listAR <- c("AR0","AR1")
nSimu <- 3
nState <- 2
nVar <- 2 # nb of state dependent variables
x0 <- rep(0,nVar) # Initialization of the time series for AR1 

maxItHmm <- 200 # It is 1000 by default
maxItHsmm <- 20 # It is 100 by default
maxIt <- max(maxItHmm,maxItHsmm)

tolHmm <- 1e-04 # It is 1e-08 by default
tolHsmm <- 1e-03 # It is 1e-04 by defaut

##############
# Building experimental design in order to set the seed to the same values at each run
##################
expDesign <- expand.grid(listVessel,listMethod,listAR,1:4)
names(expDesign) <- c("vessel","method","AR", "degrad")
expDesign$tol <- ifelse(expDesign$method=="HMM",tolHmm,tolHsmm)
expDesign$maxIt <- ifelse(expDesign$method=="HMM",maxItHmm,maxItHsmm)


MSE <- function(listExp=1:1, listSimu = 1:1,nSimu=100){
  # nSimu controls the size of the outputs
  #################
  # Formating outputs
  ##################
  res <- list(
    #loglik=array(dim=c(nSimu,maxIt),
    #         dimnames=list(paste0("simu",seq(1,nSimu)),paste0("loglik",seq(1,maxIt)))),
    loglik=as.list(1:nSimu),
    time=array(dim=c(nSimu,2),
               dimnames=list(paste0("simu",seq(1,nSimu)),c("timeSimu","timeEstim"))),
    theta=array(dim=c(nSimu,3,5,nState),
                dimnames=list(paste0("simu",seq(1,nSimu)),
                              c("vp","vr","duree"),
                              c("m","var","mu","eta","shift"),
                              c("s1","s2"))),
    confusion.ML=array(dim=c(nSimu,2,2),
                       dimnames=list(paste0("simu",seq(1,nSimu)),
                                     c("yhat1","yhat2"),
                                     c("state1","state2"))),
    confusion.vit=array(dim=c(nSimu,2,2),
                        dimnames=list(paste0("simu",seq(1,nSimu)),
                                      c("yhat1.viterbi","yhat2.viterbi"),
                                      c("state1.viterbi","state2.viterbi")))
  )
  ####################
  for(iExp in listExp){ # iExp <- 8
    vessel <- as.character(expDesign[iExp,1]) ;  iVessel <- which(listVessel == vessel)
    method <- as.character(expDesign[iExp,2]) ;  iMethod <- which(listMethod == method)
    ar <- as.character(expDesign[iExp,3]) ;  iAR <- which(listAR == ar)
    iDegrad <- expDesign[iExp,4] ;  degrad <- listDegrad[iVessel,iDegrad]
    cat(paste("Experiment :", iExp,vessel,method,ar,"degrad",iDegrad,sep=" "), "\n")
    
    #############
    # Maximum duration of a sojourn time (useful for HSMM).
    # While the negBinom as an infinite support, we truncate it to the max of the its 99% quantiles in each state
    # so that, on average, only one simulation out of 100 is truncated.
    # Old version : M <- round(5*max(parameters[iVessel,iDegrad,3,1,]))
    #############
    tempF <- function(i){
      temp <- f.mv2np(parameters[iVessel,iDegrad,3,1,i]-parameters[iVessel,iDegrad,3,5,i],parameters[iVessel,iDegrad,3,2,i])
      qnbinom(0.99,temp[1],temp[2])
    }
    M <- max(sapply(1:nState, tempF))
    #rm(tempF,temp)
    
    #####################
    # Defining the model
    #####################
    modelInit <- raw.model.builder(method = method, type = "nbinom", AR = ar, theta=parameters[iVessel,iDegrad,,,])
    
    #####################
    # Defining the length of a simulation, i.e. "nsim" parameter.
    # For HMM, simulate.hmmspec(nsim) controls the overall length of the simulation.
    #   It is thus set to 250 times the average duration of a sequence state 1 state 2.
    # For HSMM, simulate.hsmmspec(nsim) controls the number of different visited states.
    #   It is thus set 500 to get 250 alternations of state 1 and 2.
    #####################
    if(method=="HMM") nMax <- 250*round(sum(parameters[iVessel,iDegrad,3,1,]))
    if(method=="HSMM") nMax <- 500 # (modif on 29/01/2019 ; it was set to 250 since 01/03/2019)
    
    #####################
    # Simulation-Estimation
    #####################
    for(iSimu in listSimu){ # iSimu <- 1
      seed <- set.seed((iExp-1)*100+iSimu) # from 1 to 3200
      if(iExp==4 & iSimu==72) seed <- set.seed((iExp-1)*100+iSimu+10000)
      if(iExp==12 & iSimu==73) seed <- set.seed((iExp-1)*100+iSimu+10000)
      ### Simulation
      tic()
      cat('simulation n° : ',iSimu,"\t")
      if(method == "HMM") simuModel <- simulate.hmmspec(modelInit,seed=seed,nsim=rep(nMax,1),x0=x0)
      if(method == "HSMM") simuModel <- simulate.hsmmspec(modelInit,seed=seed,nsim=rep(nMax,1),x0=x0)
      temp <- toc(quiet=T)
      cat(temp$callback_msg,"\t")
      res$time[iSimu,1] <- as.numeric(temp$toc-temp$tic)
      
      ### Save the simulated the speeds 
      simuObs <- list(x=simuModel$x,N=simuModel$N)
      #plot(simuModel)
      #addStates(simuModel$s)

      ### Estimation
      tic()
      cat('estimation n° : ',iSimu,"\t")
      if(method=="HMM") estim0 <- hmmfit(simuObs,modelInit,mstep=NULL,maxit=expDesign[iExp,]$maxIt,tol=expDesign[iExp,]$tol)
      if(method=="HSMM") estim0 <- hsmmfit(simuObs,modelInit,mstep=NULL,M=M,graphical=F,maxit=expDesign[iExp,]$maxIt,tol=expDesign[iExp,]$tol)
      temp <- toc(quiet=T)
      cat(temp$callback_msg,"\n")
      res$loglik[[iSimu]] <- estim0$loglik
      res$time[iSimu,2] <- as.numeric(temp$toc-temp$tic) 
      #if(method=="HMM" & length(estim0$loglik)==maxItHmm) cat("WARNING: Estimation nb", iSimu, "reached maxit ; check it with plot(estim0$loglik) \n") 
      #if(method=="HSMM" & length(estim0$loglik)==maxItHsmm) cat("WARNING: Estimation nb", iSimu, "reached maxit ; check it with plot(estim0$loglik) \n") 
      
      ### Label switching
      estim <- label.switching(estim0,AR=ar)
      
      ### Save the outputs of each estimation
      if(ar=="AR0"){
        res$theta[iSimu,1:2,1,1:nState] <- sapply(1:nState, function(j) estim$model$parms.emission$mu[[j]])
        res$theta[iSimu,1:2,2,1:nState] <- sapply(1:nState, function(j) diag(estim$model$parms.emission$sigma[[j]]))
      } else {
        res$theta[iSimu,1:2,3,1:nState] <- sapply(1:nState, function(i) estim$model$parms.emission$mu[[i]])
        res$theta[iSimu,1:2,4,1:nState] <- sapply(1:nState, function(i) estim$model$parms.emission$eta[[i]])
        res$theta[iSimu,1:2,2,1:nState] <- sapply(1:nState, function(i) diag(estim$model$parms.emission$sigma[[i]]))
        res$theta[iSimu,1:2,1,1:nState] <- res$theta[iSimu,1:2,4,1:nState]/(1-res$theta[iSimu,1:2,3,1:nState])
      }
      if(method =="HMM"){
        res$theta[iSimu,3,1,1:nState] <- 1/(1-diag(estim$model$transition))
      } else {
        res$theta[iSimu,3,1,1:nState] <- estim$model$sojourn$shift+estim$model$sojourn$mu
        res$theta[iSimu,3,2,1:nState] <- f.np2mv(estim$model$sojourn$size,estim$model$sojourn$prob)[c(3,4)]
        res$theta[iSimu,3,5,1:nState] <- estim$model$sojourn$shift
      }
      ### Confusion matrix without Viterbi
      res$confusion.ML[iSimu,1:2,1:2] <- table(estim$yhat,simuModel$s)
      ### Confusion matrix with Viterbi
      res$confusion.vit[iSimu,1:2,1:2] <- table(predict(estim,simuObs,method="viterbi")$s,simuModel$s)
      
      save(res,file=paste("Result/ResultatsR_MSA-nbinom",iVessel,iMethod,iAR,iDegrad,"exp",iExp,"simu",iSimu,"RData",sep="."))
    }
  }
}


#############################################
# Quality control of the estimations
# Given that tol and maxit are reduced, visual check of convergence in each simulation is recommanded

qualityControl <- function(nSimu=3){
  oldMai <- par('mai')
  par(mfrow=c(4,4))
  par(mai=oldMai/2)
  for(iExp in 1:32){ # iExp <- 1
    vessel <- as.character(expDesign[iExp,1]) ;  iVessel <- which(listVessel == vessel)
    method <- as.character(expDesign[iExp,2]) ;  iMethod <- which(listMethod == method)
    ar <- as.character(expDesign[iExp,3]) ;  iAR <- which(listAR == ar)
    iDegrad <- expDesign[iExp,4] ;  degrad <- listDegrad[iVessel,iDegrad]
    nIt <- NULL
    for(iSimu in 1:nSimu){
      load(paste("Result/ResultatsR_MSA-nbinom",iVessel,iMethod,iAR,iDegrad,"exp",iExp,"simu",iSimu,"RData",sep="."))
      nIt[iSimu] <- length(res$loglik[[iSimu]])
    }
    maxItRealised <- max(nIt)
    maxItPossible <- expDesign[iExp,]$maxIt
    plot(c(1,maxItRealised),c(-1,0),main=paste(vessel,method,ar,iDegrad,sep="."),type="n",
         xlab="",ylab="",las=1,xaxs="i",yaxs="i")
    for(iSimu in 1:nSimu){
      load(paste("Result/ResultatsR_MSA-nbinom",iVessel,iMethod,iAR,iDegrad,"exp",iExp,"simu",iSimu,"RData",sep="."))
      temp <- res$loglik[[iSimu]]
      tempMin <- min(temp)
      tempMax <- max(temp)
      tempCol <- "black"
      if(nIt[iSimu]==maxItPossible) tempCol <- "red"
      if(temp[1] > temp[nIt[iSimu]]) tempCol <- "blue"
      a <- 1/(tempMax-tempMin)
      b <- tempMax/(tempMin-tempMax)
      tempY <- a*temp+b
      lines(tempY,col=tempCol)#,type="b")
      #abline(v=length(temp),lty=3)
      text(maxItRealised/2,-0.8,iExp)
    }
    #axis(1,at=nIt,line=-0.5,labels=F,col.ticks =1)
    abline(v=nIt,col="grey",lty=3)
    #temp <- table(nIt)
    #par(new=T)
    #barplot(nIt)
  }
  par(mfrow=c(1,1))
  par(mai=oldMai)
}

############
iExp2VesselMethodArDegrad <- function(iExp=1){
  vessel <- as.character(expDesign[iExp,1]) ;  iVessel <- which(listVessel == vessel)
  method <- as.character(expDesign[iExp,2]) ;  iMethod <- which(listMethod == method)
  ar <- as.character(expDesign[iExp,3]) ;  iAR <- which(listAR == ar)
  iDegrad <- expDesign[iExp,4] ;  degrad <- listDegrad[iVessel,iDegrad]
  list(vessel=vessel,method=method,ar=ar,degrad=degrad)
  list(iVessel=iVessel,iMethod=iMethod,iAR=iAR,iDegrad=iDegrad)
}


#######################################
#######################################
#######################################

MSE(1:32,nSimu=3)

qualityControl()

#######################################
# This shows that for experiment n° 19, 23, 27, 28  
# the number of iterations should be increased
# Before running the full experiment for nSimu = 100 
iExpMaxIt <- c(19,23,27,28,31,32)
expDesign[iExpMaxIt,]$maxIt <- 60

MSE(iExpMaxIt,nSimu=3)

qualityControl()

MSE(listExp = 1:32,listSimu = 1:100, nSimu=100) 
# ==> plantage à iExp=12 et iSimu=73 pour cause de divergence de l'EM (décroissance de la logL) générant des calculs de variance des temps de séjour
# Inf dans cov.wt avec la méthode "unbiased" choisie par défaut 
# 06/02/2024 : modification de cov.wt(method="ML") pour avoir un calcul en 1/N et non en 1/(N-1)
# en pondéré N=sum(wt^2) qui pour les temps de séjour est parfois = 1 
# cov.wt(method="ML") résiste à ces cas où l'EM part dans la mauvaise direction

MSE(listExp = 12, listSimu=73)

qualityControl(nSimu=100)

MSE(listExp = 12, listSimu=1:100, nSimu=100)

MSE(listExp = 13:32,listSimu = 1:100, nSimu=100)

# For iExp==4, simu 72 is strange ; extending the nb of iteration is useless ==> removed
# For iExp==12, simu 73 gets decreasing logL (?) : extending the nb of iteration is useless ==> removed
# One can also modified the seed to

tempiEx <- 4
temp <- iExp2VesselMethodArDegrad(tempiEx)
nIt <- NULL
for(iSimu in 1:100){
  load(paste("Result/ResultatsR_MSA-nbinom",temp$iVessel,temp$iMethod,temp$iAR,temp$iDegrad,"exp",tempiEx,"simu",iSimu,"RData",sep="."))
  nIt[iSimu] <- length(res$loglik[[iSimu]])
}
which(nIt==expDesign$maxIt[tempiEx])


tempiEx <- 12
temp <- iExp2VesselMethodArDegrad(tempiEx)
nIt <- NULL
for(iSimu in 1:100){
  load(paste("Result/ResultatsR_MSA-nbinom",temp$iVessel,temp$iMethod,temp$iAR,temp$iDegrad,"exp",tempiEx,"simu",iSimu,"RData",sep="."))
  nIt[iSimu] <- length(res$loglik[[iSimu]])
}
which(nIt==max(nIt))


MSE(listExp = 4, listSimu=72, nSimu=100)
MSE(listExp = 12, listSimu=73, nSimu=100)

### targeted increase of maxIt (some simu reach maxIt - red lines in qualityControl)
### iExp = c(11,15,24) ==> 40
### iExp = c(19,23,27,31) ==> 80
expDesign[c(11,15,24),]$maxIt <- 40
expDesign[c(19,23),]$maxIt <- 80
expDesign[c(27,31),]$maxIt <- 120 # the default for HSMM was 100

MSE(listExp = c(11,15,24), listSimu=1:100, nSimu=100)
MSE(listExp = c(19,23), listSimu=1:100, nSimu=100)

qualityControl(nSimu=100)

MSE(listExp = c(27,31), listSimu=1:100, nSimu=100)

qualityControl(nSimu=100)

expDesign[31,]$maxIt <- 150 
MSE(listExp = 31, listSimu=1:100, nSimu=100)
qualityControl(nSimu=100)



