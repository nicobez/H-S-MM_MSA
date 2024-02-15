ro <- seq(0.1,0.9,0.1)
rro <- sqrt(1-ro**2)

L <- 10000

#temp <- sapply(0:L, function(i) ro**i)
#temp2 <- c(1,sapply(1:L, function(i) sqrt(1-ro**2)/ro**i))
#RO <- temp%*%t(temp2)
#RO[col(RO)>row(RO)] <- 0

nSimu=100

v <- seq(50,150,by=10)
s2vinV <- NULL
A <- NULL
for(i.A in 1:length(ro)){
  for(i.v in 1:length(v)){
    tempCut <- cut(1:(L+1),breaks=seq(1,(L+1),by=v[i.v]))
    tempRes <- NULL
    for(i in 1:nSimu){
      #simu <- RO%*%rnorm(L+1)
      simu <- rnorm(1)
      for(k in 1:L) simu[k+1] <- ro[i.A]*simu[k] + rro[i.A]*rnorm(1)
      tempRes[i] <- var(tapply(simu,tempCut,FUN="mean"))
    }
    s2vinV[i.v] <- mean(tempRes)
  }
  
  plot(log10(v),log10(s2vinV),asp=1)
  lsfit(log(v),log(s2vinV))$coef
  A[i.A] <- exp(lsfit(log(v),log(s2vinV))$coef[1])
}

plot(ro,A,type="b",xlim=c(0,1),ylim=c(0,10))
lines(seq(0,1,0.01),-1/log(seq(0,1,0.01)),col=2) # integral range given by \int exp(-lambda t)dt = -1/lambda=-1/log(ro) 
lines(seq(0,1,0.01),1/(1-seq(0,1,0.01)),col=3) # integral range given by \sum ro^k = 1/(1-ro)

f.integralRange <- approxfun(x=ro,y=A,rule=2)
lines(seq(0,1,0.01),f.integralRange(seq(0,1,0.01)))

