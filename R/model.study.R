#' @export
#' 
model.study <- function(prior.type="default"){

if(prior.type == "default"){ 
string1 <- "model{
  prob[1] <- (pi.n*(1-s1) + pi.c*(1-v1))
  prob[2] <- (pi.n*s1 + pi.c*v1)
  prob[3] <- (pi.a*(1-b1))
  prob[4] <- (pi.a*b1)
  prob[5] <- (pi.n*(1-s1))
  prob[6] <- (pi.n*s1)
  prob[7] <- (pi.c*(1-u1)+pi.a*(1-b1))
  prob[8] <- (pi.c*u1+pi.a*b1)
  
  R[1:4] ~ dmulti(prob[1:4], N0)
  R[5:8] ~ dmulti(prob[5:8], N1)
  
  pi.n <- exp(n)/(1+exp(n)+exp(a))
  pi.a <- exp(a)/(1+exp(n)+exp(a))
  pi.c <- 1-pi.a-pi.n
  probit(u1) <- alpha.u
  probit(v1) <- alpha.v
  logit(s1) <- alpha.s
  logit(b1) <- alpha.b
  
  CACE <- u1-v1
"
string2 <-
  "# priors
  n ~ dnorm(0, 0.16)
  a ~ dnorm(0, 0.16)
  alpha.s ~ dnorm(0, 0.25)
  alpha.b ~ dnorm(0, 0.25)
  alpha.u ~ dnorm(0, 0.25)
  alpha.v ~ dnorm(0, 0.25)
  }
  "
modelstring <- paste(string1, string2, sep="\n")
}

  
else if (prior.type == "custom"){
  string2 <- prior.study(prior.type)
  modelstring <- paste(string1, string2, sep="\n")
}
  
if(!is.element(prior.type,c("default", "custom"))){
  stop("specified prior type should be either 'default' or 'custom'.")
}
  
  return(modelstring)
}
