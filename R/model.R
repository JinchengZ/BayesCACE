model.full <- function(prior.type="diffuse", delta.u=T, delta.v=T, 
                       delta.n=T, delta.a=T, delta.s=T, delta.b=T, pho=T){
  
  if(prior.type=="diffuse" & delta.u & delta.v & 
     delta.n & delta.a & delta.s & delta.b & pho){
    modelstring<-"model{
  for (i in 1:I) {
    
    prob[i, 1] <- (pi_n[i]*(1-s1[i]) + pi_c[i]*(1-v1[i]))
    prob[i, 2] <- (pi_n[i]*s1[i] + pi_c[i]*v1[i])
    prob[i, 3] <- (pi_a[i]*(1-b1[i]))
    prob[i, 4] <- (pi_a[i]*b1[i])
    prob[i, 5] <- (pi_n[i]*(1-s1[i]))
    prob[i, 6] <- (pi_n[i]*s1[i])
    prob[i, 7] <- (pi_c[i]*(1-u1[i])+pi_a[i]*(1-b1[i]))
    prob[i, 8] <- (pi_c[i]*u1[i]+pi_a[i]*b1[i])
    
    R[i, 1:4] ~ dmulti(prob[i, 1:4], N0[i])
    R[i, 5:8] ~ dmulti(prob[i, 5:8], N1[i])
    
    probit(u1[i]) <- alpha.u + delta.u[i]
    delta.u[i] ~ dnorm(0, tau.u)
    probit(v1[i]) <- alpha.v + delta.v[i]
    delta.v[i] ~ dnorm(0, tau.v)
    
    n[i] <- alpha.n + delta.n[i]
    delta.n[i] ~ dnorm(0, tau.n)
    a[i] <- alpha.a + delta.a[i]
    delta.a[i] ~ dnorm(0, tau.a)
    pi_n[i] <- exp(n[i])/(1+exp(n[i])+exp(a[i]))
    pi_a[i] <- exp(a[i])/(1+exp(n[i])+exp(a[i]))
    pi_c[i] <- 1-pi_a[i]-pi_n[i]
    
    logit(s1[i]) <- alpha.s1 + delta.s[i]
    delta.s[i] ~ dnorm(0, tau.s)
    logit(b1[i]) <- alpha.b1 + delta.b[i]
    delta.b[i] ~ dnorm(0, tau.b)
  } 
    
    CACE <- phi(alpha.u/sqrt(1+sigma.u^2))-phi(alpha.v/sqrt(1+sigma.v^2))
    
    pin <- exp(alpha.n)/(1+exp(alpha.n)+exp(alpha.a))
    pia <- exp(alpha.a)/(1+exp(alpha.n)+exp(alpha.a))
    pic <- 1-pia-pin
    u1out <- phi(alpha.u/sqrt(1+sigma.u^2))
    v1out <- phi(alpha.v/sqrt(1+sigma.v^2))
    s1out <- ilogit(alpha.s1)
    b1out <- ilogit(alpha.b1)
    
    # priors
    alpha.n ~  dnorm(0, 0.16)
    tau.n ~ dgamma(2, 2)
    sigma.n <- 1/sqrt(tau.n)
    alpha.a ~ dnorm(0, 0.16)
    tau.a ~ dgamma(2, 2)
    sigma.a <- 1/sqrt(tau.a)
    
    alpha.s1 ~  dnorm(0, 0.25)
    tau.s ~ dgamma(2, 2)
    sigma.s <- 1/sqrt(tau.s)
    alpha.b1 ~  dnorm(0, 0.25)
    tau.b ~ dgamma(2, 2)
    sigma.b <- 1/sqrt(tau.b)
    
    alpha.u ~  dnorm(0, 0.25)
    tau.u ~ dgamma(2, 2)
    sigma.u <- 1/sqrt(tau.u)
    alpha.v ~  dnorm(0, 0.25)
    tau.v ~ dgamma(2, 2)
    sigma.v <- 1/sqrt(tau.v)
} "
  }

  if(!is.element(prior.type,c("diffuse"))){
    stop("specified prior type is wrong.")
  }
  
  return(modelstring)
}