prior.single <- function(prior.type="custom"){
  
  string2 <-   
  "# priors
  n ~ dnorm(0, 0.01)
  a ~ dnorm(0, 0.01)
  alpha.s ~ dnorm(0, 0.01)
  alpha.b ~ dnorm(0, 0.01)
  alpha.u ~ dnorm(0, 0.01)
  alpha.v ~ dnorm(0, 0.01)
  }"

return(string2)
}
