plot.noncomp <- 
  function(data) {
  if(missing(data)) stop("need to specify a data set with complete compliance information. \n")
  
  if (any((data$n000==0 & data$n001==0 & data$n010==0 & data$n011==0) |
      (data$n100==0 & data$n101==0 & data$n110==0 & data$n111==0)) )
    stop ("noncompliance plot can only be made for complete data. \n")

  data$n00s <- data$n000+data$n001
  data$n01s <- data$n010+data$n011
  data$n10s <- data$n100+data$n101
  data$n11s <- data$n110+data$n111
  data$n0ss <- data$n00s+data$n01s
  data$n1ss <- data$n10s+data$n11s                 
  data$p10 <- data$n01s/data$n0ss   # P(T=1|R=0)
  data$p01 <- data$n10s/data$n1ss   # P(T=0|R=1)

  tmp.b <- t(sapply(split(data, data$study.id), function(x) binom.test(x$n01s, x$n0ss)$conf.int))
  tmp.a <- t(sapply(split(data, data$study.id), function(x) binom.test(x$n10s, x$n1ss)$conf.int))
  data$lower.p10 <- tmp.b[,1]
  data$upper.p10 <- tmp.b[,2]
  data$lower.p01 <- tmp.a[,1]
  data$upper.p01 <- tmp.a[,2]
  
  expit <- function(x){y=exp(x)/(1+exp(x)) 
    return(y)}
  outp10 <- rma(xi=n01s, ni=n0ss, measure="PLO", method="ML", data=data)  #p10
  p10_out <- paste(formatC(round(expit(outp10$beta), 3), format='f', digits=3 ), " (",
                   formatC(round(expit(outp10$ci.lb), 3), format='f', digits=3 ),",", 
                   formatC(round(expit(outp10$ci.ub), 3), format='f', digits=3 ), ")", sep="")
  
  outp01 <- rma(xi=n10s, ni=n1ss, measure="PLO", method="ML", data=data)  #p01
  p01_out <- paste(formatC(round(expit(outp01$beta), 3), format='f', digits=3 ), " (",
                   formatC(round(expit(outp01$ci.lb), 3), format='f', digits=3 ),",", 
                   formatC(round(expit(outp01$ci.ub), 3), format='f', digits=3 ), ")", sep="")
  
  tabletext<-cbind(
    c("Study (Author, Year)", paste(data$study), "Overall"),
    c("P(T=0|R=1)", 
      paste(formatC(round(data$p01, 3), format='f', digits=3 ), " (",
            formatC(round(data$lower.p01, 3), format='f', digits=3 ),",", 
            formatC(round(data$upper.p01, 3), format='f', digits=3 ), ")", sep=""), 
      p01_out  
    ),
    c("P(T=1|R=0)", 
      paste(formatC(round(data$p10, 3), format='f', digits=3 ), " (",
            formatC(round(data$lower.p10, 3), format='f', digits=3 ),",", 
            formatC(round(data$upper.p10, 3), digits=3 ), ")", sep=""), 
      p10_out  
    )
  )
  
  tabletext<-cbind(c("Study (Author, Year)", paste(data$study.name), "Overall"))

  forestplot(tabletext,
     graph.pos = 2,
     hrzl_lines = gpar(lwd=1, col="#444444"),
     legend_args = fpLegend(pos = list("top", "inset"=.03, "align"="horizontal")),
     legend = c("P(T=0|R=1)", "P(T=1|R=0)"),
     fn.ci_norm = c(fpDrawCircleCI, fpDrawNormalCI),
     boxsize = .15, # We set the box size to better visualize the type
     line.margin = .36, # We need to add this to avoid crowding
     is.summary=c(TRUE, rep(FALSE, nrow(data)), TRUE),
     mean = cbind(c(NA, data$p01, round(expit(outp01$beta), 3)), 
                  c(NA, data$p10, round(expit(outp10$beta), 3))),
     lower = cbind(c(NA, data$lower.p01, round(data$lower.p01, 3)), 
                   c(NA, data$lower.p10, round(data$lower.p10, 3))),
     upper = cbind(c(NA, data$upper.p01, round(data$upper.p01, 3)), 
                   c(NA, data$upper.p10, round(data$upper.p10, 3))),
     clip =c(0, 1),
     col=fpColors(box=c("darkred", "blue"), summary=c("darkred", "blue")),
     grid = structure(c(0, 0.5), gp = gpar(lty = 2, col = "#CCCCFF")), 
     xlab="Nocompliance Rates")

  invisible()
  }

