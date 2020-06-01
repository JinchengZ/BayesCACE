#' This function performs the Bayesian hierarchical model method for meta-analysis 
#' when the dataset has complete compliance information for all studies, 
#' as described in the paper Section 2.2.2, the Bayesian hierarchical model.
#' @importFrom stats update complete.cases
#' @import Rdpack
#' @import rjags
#' @import coda
#' @export 
#' @title Bayesian hierarchical models for CACE meta-analysis with complete compliance data
#' @param data a input dataset the same structure as the example data `epidural_c`, 
#' containing multiple rows referring to multiple studies in a meta-analysis. 
#' @param param a character string vector indicating the parameters to be tracked and estimated. 
#' By default the following parameters (see \code{details}) are included: \eqn{\theta^{CACE}} 
#' (\code{CACE}), \eqn{E(u_{i1})} (\code{u1out}), \eqn{E(v_{i1})} (\code{v1out}), \eqn{E(s_{i1})} (\code{s1out}), 
#' \eqn{E(b_{i1})} (\code{b1out}), \eqn{\pi_a} (\code{pia}), \eqn{\pi_n} (\code{pin}), and 
#' \eqn{\pi_c=1-\pi_a-\pi_n} (\code{pic}). 
#' Users can modify the string vector to only include parameters of interest besides \eqn{\theta^{CACE}}. 
#' @param prior.type the default priors are used by the default assignment `prior.type="default"`.
#' Like the function \code{\link{cace.study}}, weakly informative priors \eqn{\alpha_n, \alpha_a \sim 
#' N(0, 2.5^2)} and \eqn{\alpha_s, \alpha_b, \alpha_u, \alpha_v \sim N(0, 2^2)} are assigned to the 
#' means of these transformed parameters:
#' \eqn{\pi_{in}=\frac{\exp(n_i)}{1+\exp(n_i)+\exp(a_i)}}, \eqn{\pi_{ia}=\frac{\exp(a_i)}{1+\exp(n_i)+\exp(a_i)}}, 
#' where \eqn{n_i=\alpha_n+\delta_{in}}, \eqn{a_i=\alpha_a+\delta_{ia}}, \eqn{logit(s_{i1})=\alpha_s + \delta_{is}},
#' \eqn{logit(b_{i1})=\alpha_b + \delta_{ib}}, \eqn{probit(u_{i1})=\alpha_u + \delta_{iu}}, 
#' and \eqn{probit(v_{i1})=\alpha_v + \delta_{iv}}. 
#' The default prior of random effects are illstrated in `Details`.
#' Alternatively, this function allows users to specify their own prior distributions by saving a separate 
#' `R` file \code{prior.meta.R} under the same directory with the model file, and assigning the argument 
#' `prior.type = "custom"`. See example in `Details`. 
#' Users can modify the above customized file \code{prior.meta.R} to assign their preferred prior 
#' distributions. Note that same as the function \code{\link{cace.study}}, the function cannot
#' combine the default priors with partial user-defined prior distributions. Thus users need to 
#' be careful when choosing the customized priors: the pre-defined `R` file \code{prior.meta.R} must 
#' include distributions for all hyper-parameters. 
#' @param delta.n,delta.a,delta.u,delta.v,delta.s,delta.b,cor logical values indicating whether the 
#' corresponding random effect is included in the model. The default model sets all of these arguments 
#' to `TRUE`. Note that \eqn{\rho} (\code{cor}) can only be included when both \eqn{\delta_{in}} 
#' (\code{delta.n}) and \eqn{\delta_{ia}} (\code{delta.a}) are set to `TRUE`. Otherwise, a warning 
#' occurs and the model continues running by forcing `delta.n = TRUE` and `delta.a = TRUE`. 
#' @param study.specific a logical value indicating whether to calculate the study-specific 
#' \eqn{\theta^{CACE}_i}. If `TRUE`, the model will first check the logical status of arguments 
#' \code{delta.u} and \code{delta.v}. If both are `FALSE`, meaning that neither response rate \eqn{u_{i1}} 
#' or \eqn{v_{i1}} is modeled with a random effect, then the study-specific \eqn{\theta^{CACE}_i} is 
#' the same across studies. The function gives a warning and continues by making `study.specific = FALSE`. 
#' Otherwise, the study-specific \eqn{\theta^{CACE}_i} are estimated and saved as the parameter \code{cacei}.
#' 
#' @inheritParams cace.study
#' 
#' @details  
#' In a meta-analysis, the CACE can also be estimate using the joint likelihood from the Bayesian hierarchical
#' model. This method is systematically introduced in Zhou, et al. 2019. The log likelihood contribution of 
#' trial \eqn{i} is given by the equstion in \code{\link{cace.study}} `Details`, by adding a subscript \eqn{i} 
#' to each parameter. Then the log likelihood for all trials in the meta-analysis is \eqn{log{L}(\beta)=\sum_i 
#' log L_i(\beta_i)}. 
#' Because the studies are probably not exactly identical in their eligibility criteria, measurement techniques, 
#' study quality, etc., differences in methods and sample characteristics may introduce heterogeneity to the 
#' meta-analysis. One way to model the heterogeneity is to use a random-effects model.
#' Let \eqn{f(\beta_i | \beta_0, \Sigma_0)} be the distributions described above of all parameters 
#' \eqn{\beta_i=(\pi_{ia}, \pi_{in}, s_{i1}, b_{i1}, u_{i1}, v_{i1})}, where \eqn{\beta_0} is the vector 
#' of mean hyper-parameters \eqn{(\alpha_n, \alpha_a, \alpha_s, \alpha_b, \alpha_u, \alpha_v)}, and 
#' \eqn{\Sigma_0} is the diagonal covariance matrix containing \eqn{\Sigma_{ps}},  \eqn{\sigma^2_s}, 
#' \eqn{\sigma^2_b}, \eqn{\sigma^2_u} and \eqn{\sigma^2_v}.  
#' For the random effects, we have \eqn{\delta_{is} \sim N(0,\sigma^2_s)}, 
#' \eqn{\delta_{ib} \sim N(0,\sigma^2_b)}, 
#' \eqn{\delta_{iu} \sim N(0,\sigma^2_u)}, and 
#' \eqn{\delta_{iv} \sim N(0,\sigma^2_v)}, as response rates are assumed to be independent between latent classes. 
#' A \eqn{Gamma(2, 2)} hyper-prior distribution is assigned to the precision parameters \eqn{\sigma^{-2}_s}, 
#' \eqn{\sigma^{-2}_b}, \eqn{\sigma^{-2}_u} and \eqn{\sigma^{-2}_v}, which corresponds to a 95\% interval of 
#' \eqn{(0.6, 2.9)} for the corresponding standard deviations, allowing moderate heterogeneity in the 
#' response rates. In a reduced model with one of \eqn{\delta_{in}} or \eqn{\delta_{ia}} set to 0, the prior 
#' of the other precision parameter is also assumed to be \eqn{Gamma (2, 2)}, which gives moderate 
#' heterogeneity for latent compliance classes probabilities, whereas for the full model, \eqn{(\delta_{in}, 
#' \delta_{ia})^T \sim N(0, \Sigma_{ps})}, the prior for the variance-covariance matrix 
#' \eqn{\Sigma_{ps}} is \eqn{InvWishart(I, 3)}, where \eqn{I} is the identity matrix. 
#' Because the function `cace.meta.c()` is more complicated depending on the choice of random effects, 
#' as an illustration this is an example of the customized prior distributions file when assigning 
#' `delta.n = TRUE`, `delta.a = FALSE`, `delta.u = TRUE`, `delta.v = FALSE`, `delta.s = TRUE`, 
#' and `cor = FALSE` to function \code{cace.meta.c()}.
#' \preformatted{ prior.meta <- function(prior.type = "custom") {
#' string2 <-  "delta.n[i] ~ dnorm(0, tau.n)
#'   delta.u[i] ~ dnorm(0, tau.u)
#'   delta.s[i] ~ dnorm(0, tau.s) 
#'   # priors
#'   alpha.n ~  dnorm(0, 0.16)
#'   alpha.a ~ dnorm(0, 0.16)
#'   alpha.s ~  dnorm(0, 0.25)
#'   alpha.b ~  dnorm(0, 0.25)
#'   alpha.u ~  dnorm(0, 0.25)
#'   alpha.v ~  dnorm(0, 0.25)
#'   tau.n ~ dgamma(2, 2)
#'   sigma.n <- 1/sqrt(tau.n)
#'   tau.u ~ dgamma(2, 2)
#'   sigma.u <- 1/sqrt(tau.u)
#'   u1out <- phi(alpha.u/sqrt(1 + sigma.u^2))
#'   v1out <- phi(alpha.v)
#'   CACE <- u1out - v1out
#'   s1out <- ilogit(alpha.s/sqrt(1 + (16^2 * 3/(15^2 * pi^2)) * sigma.s^2))
#'   tau.s ~ dgamma(2, 2)
#'   sigma.s <- 1/sqrt(tau.s)
#'   b1out <- ilogit(alpha.b)"
#' return(string2) }
#' }
#' \eqn{\theta^{CACE}_i=u_{i1}-v_{i1}} for study \eqn{i}, so for the meta-analysis, the overall CACE is 
#' \eqn{\theta^{CACE}=E(\theta^{CACE}_i)=E(u_{i1})-E(v_{i1})}. When a random effect \eqn{\delta_{iu}} or 
#' \eqn{\delta_{iv}} is not assigned in the model, \eqn{E(u_{i1})=g^{-1}(\alpha_u)} and 
#' \eqn{E(v_{i1})=g^{-1}(\alpha_v)}. Otherwise, \eqn{E(u_{i1})} and \eqn{E(v_{i1})} can be estimated by 
#' integrating out the random effects, e.g., \eqn{E(u_{i1})=\int^{+\infty }_{-\infty }{g^{-1}(\alpha_u+t)}\sigma^{-1}_u \phi (\frac{t}{\sigma_u})dt}, 
#' where \eqn{\phi(\cdot)} is the standard Gaussian density. If the function \eqn{g(\cdot)} is the probit link, 
#' this expectation has a closed form: \eqn{E(u_{i1})= \Phi(\frac{\alpha_u}{\sqrt{1+{\sigma}^2_u}})}. 
#' If the link function \eqn{g(\cdot)} is logit, a well-established approximation 
#' \eqn{E(u_{i1}) \approx \text{logit}^{-1}(\frac{\alpha_u}{\sqrt{1+{C^2\sigma}^2_u}})} can be used, 
#' where \eqn{C=\frac{16\sqrt{3}}{15\pi}}. The above formulas also apply to \eqn{E(v_{i1})}, 
#' the expected response rate of a complier in the control group.
#' 
#' The default parameters to be monitored depend on which parameters are modeled as random effects. 
#' For example, \code{u1out} refers to \eqn{E(u_{i1})}, where for the probit link, \eqn{E(u_{i1})=
#' \Phi({\alpha_u})} if \eqn{\delta_u} is not specified in the model, and \eqn{E(u_{i1})=
#' \Phi(\frac{\alpha_u}{\sqrt{1+{\sigma}^2_u}})} when the random effect \eqn{\delta_u} is included. 
#' `DIC` is the penalized deviance, calculated as the sum of `D.bar` and `pD`, 
#' where `D.bar` is the posterior expectation of the deviance, reflecting the model 
#' fit, and `pD` reflects the effective number of parameters in the model. 
#' `D.bar` is usually lower when more parameters are included in the model, but 
#' complex models may lead to overfitting. Thus `DIC` balances the model's fit 
#' against the effective number of parameters. 
#' Generally a model with smaller DIC is preferred. However, it is difficult 
#' to conclude what constitutes an important improvement in DIC. Following 
#' Lunn D et al. 2012, we suggest that a reduction of less than 5 is not a substantial
#' improvement. 
#' When fitting models to a particular dataset, it is usually uncertain which random 
#' effect variables should be included in the model. The function `cace.meta.c()` 
#' allows users to specify candidate models with different random effects, and thus 
#' to conduct a forward/backward/stepwise model selection procedure to choose the best
#' fitting model. 
#' 
#' @examples
#' \dontrun{
#' data("epidural_c", package = "BayesCACE")
#' set.seed(123)
#' out.meta.c <- cace.meta.c(data = epidural_c, conv.diag = TRUE, 
#' mcmc.samples = TRUE, study.specific = TRUE)
#' # By calling the object smry from the output list out.meta.c, posterior estimates 
#' # (posterior mean, standard deviation, posterior median, 95\% credible interval, and 
#' # time-series standard error) are displayed.
#' out.meta.c$smry
#' out.meta.c$DIC
#' }
#' @seealso \code{\link[BayesCACE]{cace.study}}, \code{\link[BayesCACE]{cace.meta.ic}}
#' @references {
#' \insertRef{zhou2019bayesian}{BayesCACE}
#' \insertRef{zhou2020software}{BayesCACE}
#' \insertRef{lunn2012bugs}{BayesCACE}
#' \insertRef{zeger1988models}{BayesCACE}
#' } 

cace.meta.c <-
  function(data, 
           param = c("CACE", "u1out", "v1out", "s1out", "b1out", 
                   "pic", "pin", "pia"),
           prior.type = "default", 
           delta.n = TRUE, delta.a = TRUE, delta.u = TRUE, delta.v = TRUE, 
           delta.s = TRUE, delta.b = TRUE, cor = TRUE, 
           digits = 3, n.adapt = 1000, n.iter = 100000,
           n.burnin = floor(n.iter/2), n.chains = 3, n.thin = max(1,floor((n.iter-n.burnin)/100000)),
           conv.diag = FALSE, mcmc.samples = FALSE, study.specific = FALSE)    {
    ## check the input parameters
    options(warn=1)
    
    if(missing(data)) stop("need to specify data")
    if(!missing(data) ){
      study.id <- data$study.id[complete.cases(data)]
      n000<-data$n000[complete.cases(data)]
      n001<-data$n001[complete.cases(data)]
      n010<-data$n010[complete.cases(data)]
      n011<-data$n011[complete.cases(data)]
      n100<-data$n100[complete.cases(data)]
      n101<-data$n101[complete.cases(data)]
      n110<-data$n110[complete.cases(data)]
      n111<-data$n111[complete.cases(data)]
      cat("NA is not allowed in the input data set;\n")
      cat("the rows containing NA are removed.\n")
    }

    if(length(study.id)!=length(n000) | length(n000)!=length(n001) | length(n001)!=length(n010) | 
       length(n010)!=length(n011) | length(n011)!=length(n100) | length(n100)!=length(n101) |
       length(n101)!=length(n110) | length(n110)!=length(n111) )
      stop("study.id, n000, n001, n010, n011, n100, n101, n110, and n111 have different lengths. \n")
    
    if ((!(delta.n & delta.a)) & cor){
      warning("'cor' can be assigned as TRUE only if both delta.n and delta.a are TRUE.\n
              the model is continued by forcing delta.n=TRUE and delta.a=TRUE")
      delta.n=TRUE
      delta.a=TRUE
    }
    
    if ((!(delta.u|delta.v)) & study.specific){
      warning("no random effect is assigned to the response rate u1 or v1, \n
           study-specific CACE is the same across studies. \n
           the model is continued by making 'study.specific=FALSE'. \n
           to make a CACE forestplot, please run 'cace.study' to estimate study level CACEs. \n")
      study.specific=FALSE
    }
      
    
    Ind <- rep(1, 7)
    if (!delta.n) Ind[1]=0
    if (!delta.a) Ind[2]=0
    if (!delta.u) Ind[3]=0
    if (!delta.v) Ind[4]=0
    if (!delta.s) Ind[5]=0
    if (!delta.b) Ind[6]=0
    if (!cor) Ind[7]=0
    
    ## jags model
    modelstring<-model.meta.c(prior.type, Ind)

    ## jags data
    if(prior.type == "default"){
      Ntol <- n000+n001+n010+n011+n100+n101+n110+n111
      N0 <- n000+n001+n010+n011
      N1 <- n100+n101+n110+n111
      R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111)
      I <- length(Ntol)
      pi <- pi
      data.jags <- list(N0=N0, N1=N1, R=R, I=I, Ind=Ind, pi=pi)
    }

    
    ## jags initial value
    rng.seeds<-sample(1000000,n.chains)
    init.jags <- vector("list", n.chains)
    for(ii in 1:n.chains){
      init.jags[[ii]] <- list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
    }
    
    ## parameters to be paramed in jags
    if(!is.element("CACE",param)) param<-c("CACE",param)

    fullparam <- c("CACE", "u1out", "v1out", "s1out", "b1out", 
                   "pic", "pin", "pia", "alpha.n", "alpha.a", 
                   "alpha.u", "alpha.v", "alpha.s", "alpha.b", 
                   "sigma.n", "sigma.a", "rho", "Sigma.rho",
                   "sigma.u", "sigma.v", "sigma.s", "sigma.b")
    if(!any(is.element(param, fullparam))) stop("parameters must be specified from the following:  
                  CACE, u1out, v1out, s1out, b1out, 
                   pic, pin, pia, alpha.n, alpha.a, 
                alpha.u, alpha.v, alpha.s, alpha.b, 
                sigma.n, sigma.a, rho, Sigma.rho,
                sigma.u, sigma.v, sigma.s, sigma.b")
    if(study.specific) param<-c("cacei",param)
    
    ## run jags
    cat("Start running MCMC...\n")
    jags.m<-jags.model(file=textConnection(modelstring),data=data.jags,inits=init.jags,
                       n.chains=n.chains,n.adapt=n.adapt)
    update(jags.m,n.iter=n.burnin)
    jags.out<-coda.samples.dic(model=jags.m,variable.names=param,n.iter=n.iter,thin=n.thin)
    
    out<-NULL
    out$Ind <- Ind
    out$model<-"cace.meta.c"
    
    smry<-summary(jags.out$samples)
    smry<-cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")], 
                smry$statistics[,c("Naive SE","Time-series SE")])
    smry<-signif(smry,digits=digits)
    out$smry <- smry

    #dic
    dev<-jags.out$dic[[1]] # mean deviance
    pen<-jags.out$dic[[2]] # pD
    pen.dev<-dev+pen # DIC
    dic.stat<-rbind(dev,pen,pen.dev)
    rownames(dic.stat)<-c("D.bar","pD","DIC")
    colnames(dic.stat)<-""
    out$DIC<-dic.stat

    for (i in 1:length(fullparam)){
      if(is.element(fullparam[i],param)) {out[[fullparam[i] ]]<-smry[c(fullparam[i]), ]}
    } 
    
    if(conv.diag){
      cat("MCMC convergence diagnostic statistics are calculated and saved in conv.out\n")
      conv.out<-gelman.diag(jags.out$samples,multivariate=FALSE)
      out$conv.out<-conv.out$psrf
    }
    
    if(mcmc.samples){
      out$mcmc.samples<-jags.out$samples
    }

    class(out)<-"cace.Bayes"
    return(out)
    options(warn=0)   
}


