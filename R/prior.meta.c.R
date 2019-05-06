prior.meta <- function(prior.type="custom"){

string2 <-   
" delta.n[i] ~ dnorm(0, tau.n)
  delta.u[i] ~ dnorm(0, tau.u)
  delta.s[i] ~ dnorm(0, tau.s)
  }

# priors
  alpha.n ~  dnorm(0, 0.16)
  alpha.a ~ dnorm(0, 0.16)
  alpha.s ~  dnorm(0, 0.25)
  alpha.b ~  dnorm(0, 0.25)
  alpha.u ~  dnorm(0, 0.25)
  alpha.v ~  dnorm(0, 0.25)

  tau.n ~ dgamma(2, 2)
  sigma.n <- 1/sqrt(tau.n)
  
  tau.u ~ dgamma(2, 2)
  sigma.u <- 1/sqrt(tau.u)
  u1out <- phi(alpha.u/sqrt(1+sigma.u^2))
  v1out <- phi(alpha.v)
  CACE <- u1out-v1out
  
  s1out <- ilogit(alpha.s/sqrt(1 + (16^2*3/(15^2*pi^2))*sigma.s^2))
  tau.s ~ dgamma(2, 2)
  sigma.s <- 1/sqrt(tau.s)
  
  b1out <- ilogit(alpha.b)
  }"

return(string2)
}
