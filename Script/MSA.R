##################
# Code for MSA experiments (estimation model = simulation model)
##################


##################
# Working variables 
##################

nState <- 2
nVar <- 2 # nb of state dependent variables
x0 <- rep(0,nVar) # Initialization of the time series for AR1 

maxItHmm <- 200 # It is 1000 by default
maxItHsmm <- 150 # It is 100 by default
maxIt <- max(maxItHmm,maxItHsmm)

tolHmm <- 1e-04 # It is 1e-08 by default
tolHsmm <- 1e-03 # It is 1e-04 by defaut

##################
# Building experimental design in order to set the seed to the same values at each run
##################

expDesign <- expand.grid(listVessel,listMethod,listAR,1:4)
names(expDesign) <- c("vessel","method","AR", "degrad")
expDesign$tol <- ifelse(expDesign$method=="HMM",tolHmm,tolHsmm)
expDesign$maxIt <- ifelse(expDesign$method=="HMM",maxItHmm,maxItHsmm)

##################
# Building the function for running MSA
# listExp: vector of the experiment to run
# listSimu: vector of the simulation-estimation to run (between 1 and nSimu)
# outRes: boolean whether the output should be returned or not at the end. 
#         This is useful for one particular experiment and one particular simulation.
#         Else, only the last simu of the experiment is returned.
#         Default is F.
#
# Each experiment x simulation is saved into an external RData file
##################

MSE <- function(listExp=1:1, listSimu = 1:1, outRes = F){
  #################
  # Formating outputs for ONE simu of a given experiment
  ##################
  res <- list(
    loglik=NULL,
    time=list("simu"=NULL,"estim"=NULL),
    theta=array(dim=c(3,5,nState),
                dimnames=list(c("vp","vr","duree"),
                              c("m","var","mu","eta","shift"),
                              c("s1","s2"))),
    confusion.ML=array(dim=c(2,2),
                       dimnames=list(c("yhat1","yhat2"),
                                     c("state1","state2"))),
    confusion.vit=array(dim=c(2,2),
                        dimnames=list(c("yhat1.viterbi","yhat2.viterbi"),
                                      c("state1.viterbi","state2.viterbi")))
  )
  for(iExp in listExp){ # iExp <- 8
    vessel <- as.character(expDesign[iExp,1]) ;  iVessel <- which(listVessel == vessel)
    method <- as.character(expDesign[iExp,2]) ;  iMethod <- which(listMethod == method)
    ar <- as.character(expDesign[iExp,3]) ;  iAR <- which(listAR == ar)
    iDegrad <- expDesign[iExp,4] ;  degrad <- listDegrad[iVessel,iDegrad]
    cat(paste("Experiment :", iExp,vessel,method,ar,"degrad",iDegrad,sep=" "), "\n")
    #############
    # Definition of the maximum duration of a sojourn time (useful for HSMM).
    # While the negBinom as an infinite support, we truncate it to the max of the its 99% quantiles in each state
    # so that, on average, only one simulation out of 100 is truncated.
    # Old version : M <- round(5*max(parameters[iVessel,iDegrad,3,1,]))
    #############
    tempF <- function(i){
      temp <- f.mv2np(parameters[iVessel,iDegrad,3,1,i]-parameters[iVessel,iDegrad,3,5,i],parameters[iVessel,iDegrad,3,2,i])
      qnbinom(0.99,temp[1],temp[2])
    }
    M <- max(sapply(1:nState, tempF))
    #####################
    # Definition of the model
    #####################
    modelInit <- raw.model.builder(method = method, type = "nbinom", AR = ar, theta=parameters[iVessel,iDegrad,,,])
    #####################
    # Definition of the length of a simulation, i.e. "nsim" parameter of simulate.hmmspec() and simulate.hsmmspec().
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
      res$time$simu <- as.numeric(temp$toc-temp$tic)
      ### Save the simulated the speeds 
      simuObs <- list(x=simuModel$x,N=simuModel$N)
      ### Estimation
      tic()
      cat('estimation n° : ',iSimu,"\t")
      if(method=="HMM") estim0 <- hmmfit(simuObs,modelInit,mstep=NULL,maxit=expDesign[iExp,]$maxIt,tol=expDesign[iExp,]$tol)
      if(method=="HSMM") estim0 <- hsmmfit(simuObs,modelInit,mstep=NULL,M=M,graphical=F,maxit=expDesign[iExp,]$maxIt,tol=expDesign[iExp,]$tol)
      temp <- toc(quiet=T)
      cat(temp$callback_msg,"\n")
      res$loglik <- estim0$loglik
      res$time$estim <- as.numeric(temp$toc-temp$tic) 
      ### Label switching
      estim <- label.switching(estim0,AR=ar)
      ### Save the outputs of each estimation
      if(ar=="AR0"){
        res$theta[1:2,1,1:nState] <- sapply(1:nState, function(j) estim$model$parms.emission$mu[[j]])
        res$theta[1:2,2,1:nState] <- sapply(1:nState, function(j) diag(estim$model$parms.emission$sigma[[j]]))
      } else {
        res$theta[1:2,3,1:nState] <- sapply(1:nState, function(i) estim$model$parms.emission$mu[[i]])
        res$theta[1:2,4,1:nState] <- sapply(1:nState, function(i) estim$model$parms.emission$eta[[i]])
        res$theta[1:2,2,1:nState] <- sapply(1:nState, function(i) diag(estim$model$parms.emission$sigma[[i]]))
        res$theta[1:2,1,1:nState] <- res$theta[1:2,4,1:nState]/(1-res$theta[1:2,3,1:nState])
      }
      if(method =="HMM"){
        res$theta[3,1,1:nState] <- 1/(1-diag(estim$model$transition))
      } else {
        res$theta[3,1,1:nState] <- estim$model$sojourn$shift+estim$model$sojourn$mu
        res$theta[3,2,1:nState] <- f.np2mv(estim$model$sojourn$size,estim$model$sojourn$prob)[c(3,4)]
        res$theta[3,5,1:nState] <- estim$model$sojourn$shift
      }
      ### Confusion matrix without Viterbi
      res$confusion.ML[1:2,1:2] <- table(estim$yhat,simuModel$s)
      ### Confusion matrix with Viterbi
      res$confusion.vit[1:2,1:2] <- table(predict(estim,simuObs,method="viterbi")$s,simuModel$s)
      ### Saving the final resul into a dedicated RData
      save(res,file=paste("Result/ResultatsR_MSA-nbinom",iVessel,iMethod,iAR,iDegrad,"exp",iExp,"simu",iSimu,"RData",sep="."))
    }
  }
  if(outRes) res
}


##################
# Building the function for quality control of the estimations.
# setting: number of the setting for which the control is performed (setting = 1 or 2)
# nSimu: number of simulations represented
# save: save=T means the plot is saved and overright any prior file with the same name
#
# Given that tol and maxit are reduced, visual check of convergence in each simulation is recommended
# Cases for which the log-likelihood decreases along the EM are represented in BLUE
# Cases for which the maximum possible number of iterations is reached are represented in RED
##################
qualityControl <- function(setting = 1, nSimu=3, save=T){
  listExp <- setting + c(0,8,16,24,4,12,20,28,2,10,18,26,6,14,22,30)
  oldMai <- par('mai')
  m <- matrix(1:20,4,5,byrow=T)
  layout(m, widths=c(1,rep(4,4)),heights=rep(4,4)) #layout.show(20)
  par(mai=c(0.3,0.4,0.1,0.05))#rep(0.075,4))
  k <- 0
  for(iExp in listExp){ # iExp <- 1
    vessel <- as.character(expDesign[iExp,1]) ;  iVessel <- which(listVessel == vessel)
    method <- as.character(expDesign[iExp,2]) ;  iMethod <- which(listMethod == method)
    ar <- as.character(expDesign[iExp,3]) ;  iAR <- which(listAR == ar)
    iDegrad <- expDesign[iExp,4] ;  degrad <- listDegrad[iVessel,iDegrad]
    nIt <- NULL
    for(iSimu in 1:nSimu){ #iSimu <- 1
      load(paste("Result/ResultatsR_MSA-nbinom",iVessel,iMethod,iAR,iDegrad,"exp",iExp,"simu",iSimu,"RData",sep="."))
      nIt[iSimu] <- length(res$loglik)
    }
    maxItRealised <- max(nIt)
    maxItPossible <- expDesign[iExp,]$maxIt
    if(k%%5 == 0){
      plot(0,0,bty="n",type="n",xaxt="n",yaxt="n")
      text(-0.1,0,paste(method,ar,sep=","),lwd=2,cex=1.25,srt=90)
      k <- k+1
    }
    plot(c(1,maxItRealised),c(-1,0),type="n",xlab="",ylab="",las=1,xaxs="i",yaxs="i")
    for(iSimu in 1:nSimu){
      load(paste("Result/ResultatsR_MSA-nbinom",iVessel,iMethod,iAR,iDegrad,"exp",iExp,"simu",iSimu,"RData",sep="."))
      temp <- res$loglik
      tempMin <- min(temp) ; tempMax <- max(temp)
      tempCol <- vesselCol[setting]
      if(nIt[iSimu]==maxItPossible) tempCol <- "black"
      if(temp[1] > temp[nIt[iSimu]]) tempCol <- "black"
      a <- 1/(tempMax-tempMin)
      b <- tempMax/(tempMin-tempMax)
      tempY <- a*temp+b
      lines(tempY,col=tempCol)#,type="b")
      text(maxItRealised/2,-0.8,paste0("Experiment: ",iExp))
    }
    axis(3,at=nIt,line=0,labels=F,col.ticks = vesselCol[setting])
    k <- k+1
  }
  if(save) dev.print(device = png, file = paste0("Result/figureQualityControl",setting,".png"),width=600,height=600)
  par(mai=oldMai)
}

############
iExp2VesselMethodArDegrad <- function(iExp=1){
  vessel <- as.character(expDesign[iExp,1]) ;  iVessel <- which(listVessel == vessel)
  method <- as.character(expDesign[iExp,2]) ;  iMethod <- which(listMethod == method)
  ar <- as.character(expDesign[iExp,3]) ;  iAR <- which(listAR == ar)
  iDegrad <- expDesign[iExp,4] ;  degrad <- listDegrad[iVessel,iDegrad]
  list(iVessel=iVessel,iMethod=iMethod,iAR=iAR,iDegrad=iDegrad,vessel=vessel,method=method,ar=ar,degrad=degrad)
}


#######################################
#######################################
#######################################

# Check that it runs properly for all experiments WITH ONLY ONE SIMU
MSE(listExp=1:32,listSimu=1:1)
qualityControl(setting=1,nSimu=1,save=F)
qualityControl(setting=2,nSimu=1,save=F)

# Run it for nSimu = 100
# !!!!!! BE CAREFULL : it lasts 10 hours and 20 minutes !!!!!
# It is thus commented by default.
# Please, uncomment if you really want to run the MSA with 32 experiments and 100 simulations
# MSE(1:32,1:100)
qualityControl(setting=1,nSimu=100,save=T)
qualityControl(setting=2,nSimu=100,save=T)

# Extract the computing times
simuTime <- 0
estimTime <- 0
for(iExp in 1:32){ # iExp <- 1
  vessel <- as.character(expDesign[iExp,1]) ;  iVessel <- which(listVessel == vessel)
  method <- as.character(expDesign[iExp,2]) ;  iMethod <- which(listMethod == method)
  ar <- as.character(expDesign[iExp,3]) ;  iAR <- which(listAR == ar)
  iDegrad <- expDesign[iExp,4] ;  degrad <- listDegrad[iVessel,iDegrad]
  for(iSimu in 1:100){ #iSimu <- 1
    load(paste("Result/ResultatsR_MSA-nbinom",iVessel,iMethod,iAR,iDegrad,"exp",iExp,"simu",iSimu,"RData",sep="."))
    simuTime <- simuTime + res$time$simu
    estimTime <- estimTime + res$time$estim
  }
}
simuTime/3600
estimTime/3660
(simuTime+estimTime)/3600

