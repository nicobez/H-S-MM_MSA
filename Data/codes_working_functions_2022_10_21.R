
# Working Function f.* ----------------------------------------------------

f.lvrais_marginale <- function( s, Z,type="poisson"){
  # (marginal) Log-Likelyhood of shifted data for Poisson or NegBinom 
  # If the variance and the mean of the shifted data are not consistent with a NegBinom, the LL is -Inf
  N    <- length(Z)
  Z <- Z-s
  Zm <- mean(Z)
  Zvar <- var(Z)
  if(type=="poisson")
    res <- N  * (Zm) * ( 1-log(Zm)) - sum(log(gamma(Z+1)))
  if(type=="nbinom"){
    n <- Zm^2/(Zvar-Zm)
    p <- Zm/Zvar
    res <- -Inf
    if(Zvar > Zm) {
      N0 <- sum(Z==0)
      res <- N0 + sum(log(1/(Z[Z>0]*beta(n,Z[Z>0]))))+N*n*log(p)+log(1-p)*N*Zm
    }
  }
  res   
}

f.shift_mle <- function(Z,type="poisson"){
  # ML estimation of the shift of Poisson or NegBinom 
  res <- which.max(sapply(0:min(Z), function(s_){f.lvrais_marginale(s_,Z,type=type) }))
  (0:min(Z))[res]
}

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


f.hellinger.geom.vs.shifted.nbinom <- function(m,v,shift){
	x <- 0:floor(10*m)
	p <- 1/m ; q <- 1-p
	Gi <- q^x*p
	#Gi <- dgeom(x,prob=1/(1+m)) # WARNING rgeom(x,p) generates 0 ! and mean = 1/p
	tmp <- f.mv2np(m-shift,v)
	NBi <- c(rep(0,shift),dnbinom(x[1:(length(x)-shift)],size=tmp[1],prob=tmp[2]))
	sqrt(1-sum(sqrt(Gi*NBi)))
}


f.sse.nbinom <- function(v,x){
  # sum of square errors between empirical proportions and nbinom pdf of parameters v=c(size,prob)  
  #sum((dnbinom(x[,1],size=v[1],prob=v[2])-x[,2])^2)
  if (v[2] <0) cat(v,"\n")
  sum((pnbinom(x[,1],size=v[1],prob=v[2])-cumsum(x[,2]))^2)
  
}
f.sse.nbinom2 <- function(v,x){
  # sum of square errors between empirical proportions and nbinom pdf of parameters v=c(size,prob)  
  #sum((dnbinom(x[,1],size=v[1],prob=v[2])-x[,2])^2)
  sum((pnbinom(x[,1]-v[3],size=v[1],prob=v[2])-cumsum(x[,2]))^2)
}

f.est.nbinom <- function(d){
  # Inference of shifted.nbinom parameters by least squares.
  # For each possible shift value, SSE is computed and the least one is selected.
  # The optim requires initial values. To avoid problems for nbinom when the empiriacl mean is larger than the variance,
  # initial values are rather set to :
  #   size = n ~ m²/var instead of m²/(var-m) that can be < 0 when var < m
  #   prob = p = n/(m+n) that corresponds to the theoretical values that naturally lies in [0,1]
  # Input : d is the vector of the sojourn times
  d1 <- sort(unique(d))
  p1 <- sapply(d1, function(i){sum(d==i)})/length(d)
  tmp <- which.min(sapply(0:min(d1), function(s){optim(c((mean(d)-s)^2/f.var(d),1/(1+f.var(d)/(mean(d)-s))),f.sse.nbinom,x=cbind(d1-s,p1),lower=c(0,0),upper=c(500,1),method="L-BFGS-B")$value}))
  s <- (0:min(d1))[tmp]
  res <- optim(c((mean(d)-s)^2/f.var(d),1/(1+f.var(d)/(mean(d)-s))),
               f.sse.nbinom,x=cbind(d1-s,p1))$par
  list(size=res[1],prob=res[2],shift=s)
}

f.sse.poisson <- function(v,x){
  # sum of square errors between empirical proportions and poisson pdf of parameters v=lambda  
  #sum((dpois(x[,1],lambda=v[1])-x[,2])^2)
  sum((ppois(x[,1],lambda=v[1])-cumsum(x[,2]))^2)
  
}

f.est.poisson <- function(d){
  # Inference of shifted.nbinom parameters by least squares.
  # For each possible shift value, SSE is computed and the least one is selected.
  # The optim is initialised to the empirical mean.
  # Input : d is the vector of the sojourn times
  d1 <- sort(unique(d))
  p1 <- sapply(d1, function(i){sum(d==i)})/length(d)
  tmp <- which.min(sapply(0:min(d1), function(s){optim(mean(d)-s,f.sse.poisson,x=cbind(d1-s,p1),method="Brent",lower=0,upper=mean(d))$value}))
  s <- (0:min(d1))[tmp]
  res <- optim(mean(d)-s,f.sse.poisson,x=cbind(d1-s,p1),method="Brent",lower=0,upper=mean(d))$par
  list(lambda=res[1],shift=s)
}

f.var <- function(x,na.rm=T){
  # Computing proper empirical variance (the R function var() provides an estimated variance in 1/(n-1))
  mean(x^2,na.rm=na.rm)-mean(x,na.rm=na.rm)^2
}
f.cov <- function(x,y,na.rm=T){
  # Computing proper empirical covariance (the R function cov() provides an estimated covariance in 1/(n-1))
  mean(x*y,na.rm=na.rm)-mean(x,na.rm=na.rm)*mean(y,na.rm=na.rm)
}

f.accuracy <- function(M) {
  # Accuracy associated to a confusion matrix
  (M[1,1]+M[2,2])/sum(M)
}

f.odd <- function(x) {
  # Provides TRUE if x is odd
  x%%2 == 1
}

f.tau <- function(x){
  # Overlapping of speeds' Gaussian pdfs between states (one value per speed variable) 
  # x is an array with dim= 4*5*2 of which we use [1:2,1:2,]
  res <- c(NA,NA)
  for(j in 1:2){
    m1 <- x[j,1,1]
    s1 <- sqrt(x[j,2,1])
    m2 <- x[j,1,2]
    s2 <- sqrt(x[j,2,2])
    #x.inter <- (m1*s2+m2*s1)/(s1+s2)
    b <- 2*(m2*s1^2 - m1*s2^2)
    a <- s2^2-s1^2
    c <- m1^2*s2^2-m2^2*s1^2-2*s1^2*s2^2*(log(s2)-log(s1))
    delta <- b^2 -4*a*c
    x1 <- (-b+sqrt(delta))/(2*a)
    x2 <- (-b-sqrt(delta))/(2*a)
    if(m1==m2){
      temp1 <- sort(c(x1,x2))
      x1 <- temp1[1] ; x2 <- temp1[2] # make sure x1 < x2
      if(s2 > s1){
        tau <- pnorm(x2,m2,s2)-pnorm(x1,m2,s2)+(1-pnorm(x2,m1,s1)+pnorm(x1,m1,s1))
      } else {
        tau <- pnorm(x2,m1,s1)-pnorm(x1,m1,s1)+(1-pnorm(x2,m2,s2)+pnorm(x1,m2,s2))
      }
    }
    if(m1<m2){
      x.inter <- (x1 > m1 & x1 < m2)*x1 + (x2 > m1 & x2 < m2)*x2 
      tau <- 1-pnorm(x.inter,mean=m1,sd=s1)+pnorm(x.inter,mean=m2,sd=s2)
    } else {
      x.inter <- (x1 > m2 & x1 < m1)*x1 + (x2 > m2 & x2 < m1)*x2
      tau <- 1-pnorm(x.inter,mean=m2,sd=s2)+pnorm(x.inter,mean=m1,sd=s1)
    }
    res[j] <- tau
  } 
  res
}

f.curvedarrow <- function(from, to, lwd=2, lty=1, lcol="black", arr.col=lcol, 
                          arr.pos=0.5, curve=1, dr=0.01, endhead=FALSE, segment = c(0,1), ...)   {
  # copy from diagram package
  dpos  <- to-from
  angle <- atan(dpos[2]/dpos[1])*180/pi         # angle between both
  if (is.nan(angle)) return
  mid   <- 0.5*(to+from)                        # midpoint of ellipsoid arrow
  dst   <- dist(rbind(to, from))                # distance from-to
  ry    <- curve*dst                            # small radius of ellepsoid
  aFrom<-0                                      # angle to and from
  aTo<-pi
  if ( from[1] <= to[1]) {
    aFrom <- pi
    aTo <- 2*pi
  }
  
  if (segment [1] != 0)
    From <- segment[1] * aTo + (1-segment[1]) * aFrom
  else
    From <- aFrom
  
  if (segment [2] != 1)
    To <- segment[2] * aTo + (1-segment[2]) * aFrom
  else
    To <- aTo
  
  meanpi <- arr.pos * aTo + (1-arr.pos) * aFrom
  if (endhead) To <- meanpi
  
  
  plotellipse(rx=dst/2,  ry=ry, mid=mid, angle=angle, from = From, to = To,
              lwd=lwd, lty=lty, lcol=lcol)
  ell <- getellipse(rx=dst/2, ry=ry, mid=mid, angle=angle,
                    from=1.001*meanpi, to=0.999*meanpi, dr= 0.002)       #Changed from -0.002
  Arrows(ell[1,1], ell[1,2], ell[nrow(ell),1], ell[nrow(ell),2],
         code=1, lcol=lcol, arr.col=arr.col, ...)
  curvedarrow <- c(ell[nrow(ell),1], ell[nrow(ell),2])
}


# # # modification de la fonction radarchart{fmsb}
### notamment ajout de la variable "threshold"
f.myradar <- function (df, axistype = 0, seg = 4, pty = 16, pcol = 1:8, plty = 1:6, 
                       plwd = 1, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 3, 
                       cglwd = 1, cglcol = "navy", axislabcol = "blue", title = "", 
                       maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlabels = NULL, 
                       vlcex = NULL, caxislabels = NULL, calcex = NULL, paxislabels = NULL, 
                       palcex = NULL, threshold = NULL,...) 
{
  if (!is.data.frame(df)) {
    cat("The data must be given as dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  if(maxmin){dfmax <- unique(as.numeric(df[1,]))}
  if(is.null(threshold)) {plot(c(-1.2, 1.2), c(-1.2, 1.2), type = "n", frame.plot = FALSE, 
                               axes = FALSE, xlab = "", ylab = "", main = title, asp = 1, 
                               ...)}
  if(!is.null(threshold) & maxmin){plot(c(-2, 2), c(-2,2), type = "n", frame.plot = FALSE, 
                                        axes = FALSE, xlab = "", ylab = "", main = title, asp = 1, 
                                        ...)}
  theta <- seq(90, 450, length = n + 1) * pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xx * (i + CGap)/(seg + CGap), yy * (i + CGap)/(seg + 
                                                             CGap), lty = cglty, lwd = cglwd, border = cglcol)
    if (axistype == 1 | axistype == 3) 
      CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5) 
      CAXISLABELS <- sprintf("%3.2f", i/seg)
    if (!is.null(caxislabels) & (i < length(caxislabels))) 
      CAXISLABELS <- caxislabels[i + 1]
    if (axistype == 1 | axistype == 3 | axistype == 4 | 
        axistype == 5) {
      if (is.null(calcex)) 
        text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
             col = axislabcol)
      else text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS, 
                col = axislabcol, cex = calcex)
    }
  }
  if(!is.null(threshold)){
    if(centerzero){
      polygon(xx * 2, yy * 2, lty = 3, lwd = cglwd, border = cglcol)
    }
    if(!centerzero){
      polygon(xx * (2*seg+1)/(1+seg), yy * (2*seg+1)/(1+seg), lty = 3, lwd = cglwd, border = cglcol)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1, lwd = cglwd, lty = cglty, 
           length = 0, col = cglcol)
  }
  else {
    arrows(xx/(seg + CGap), yy/(seg + CGap), xx * 1, yy * 
             1, lwd = cglwd, lty = cglty, length = 0, col = cglcol)
  }
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels)) 
    PAXISLABELS <- paxislabels
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex)) 
      text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol)
    else text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol, 
              cex = palcex)
  }
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) 
    VLABELS <- vlabels
  if (is.null(vlcex)){ 
    if(is.null(threshold)) text(xx * 1.2, yy * 1.2, VLABELS)
    if(!is.null(threshold)) text(xx * 2.05, yy * 2.05, VLABELS)
  }
  else {
    if(is.null(threshold)) text(xx * 1.2, yy * 1.2, VLABELS, cex = vlcex)
    if(!is.null(threshold)) text(xx * 2.05, yy * 2.05, VLABELS, cex = vlcex)
  }
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg + CGap) + (df[i, ] - df[2, ])/(df[1, 
                                                         ] - df[2, ]) * seg/(seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n", i, 
                  df[i, ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1, 
                              1)
            }
            xxleft <- xx[left] * CGap/(seg + CGap) + 
              xx[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            yyleft <- yy[left] * CGap/(seg + CGap) + 
              yy[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            xxright <- xx[right] * CGap/(seg + CGap) + 
              xx[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + 
                                                                                            CGap)
            yyright <- yy[right] * CGap/(seg + CGap) + 
              yy[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + 
                                                                                            CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[j] <- xx[j] * (yyleft * xxright - yyright * 
                                 xxleft)/(yy[j] * (xxright - xxleft) - 
                                            xx[j] * (yyright - yyleft))
            yys[j] <- (yy[j]/xx[j]) * xxs[j]
          }
          else {
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j] * CGap/(seg + CGap) + xx[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, 
                                                 j]) * seg/(seg + CGap)
          yys[j] <- yy[j] * CGap/(seg + CGap) + yy[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, 
                                                 j]) * seg/(seg + CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], col = pfcols[i - 
                                                                                                      2])
      }
      else {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], density = pdensities[i - 
                                                                                                              2], angle = pangles[i - 2], col = pfcols[i - 
                                                                                                                                                         2])
      }
      points(xx * scale, yy * scale, pch = ptys[((i-1)*dim(df)[2]+1):(i*dim(df)[2])], 
             col = pcols[i - 2],cex=1.75)
    }
  }
}



