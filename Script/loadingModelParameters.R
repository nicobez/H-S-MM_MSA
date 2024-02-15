# Settings -------------------------------------------
listVessel <- c(924675,"007")               
listDegrad <- rbind(c(1,2,4,8),c(1,3,6,12))
listMethod <- c("HMM","HSMM")
listAR <- c("AR0","AR1")
listSize <-c(50,250)
nSimu <- 100
nVessel <- length(listVessel)

# Full experimental design i.e. for 2 vessels -----------------------------
expDesign <- expand.grid(listVessel,listMethod,listAR,1:4)
dimnames(expDesign)[[2]] <- c("vessel","method","AR", "degrad")
expDesignSize <- dim(expDesign)[1]

#####################
# Raw parameters
#####################

# Loading model parameters -------------------------------------------
# Model parameters are fitted once for all on the basis of raw data.
# To speed up the process, they are uploaded without re-estimating them. 
# The code for the estimation can be obtained on demand. 
load("Data/parameters.nbinom")
parameters <- parameters.nbinom
rm(parameters.nbinom)

