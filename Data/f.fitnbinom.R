f.fitnbinom <- function (eta) 
{
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