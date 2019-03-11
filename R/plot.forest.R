plot.forest <- 
  function(data, obj) {
    if(missing(obj)) stop("need to specify obj, the object generated from one of the following functions:
         cace.meta.c, cace.meta.ic.\n")
    if (!inherits(obj, "cace.Bayes"))
      stop("Use only with 'cace.Bayes' objects, the output generated from one of the following functions:
           cace.meta.c, cace.meta.ic.\n")
    
    if(missing(data)) stop("need to specify data")
    
    if(!class(obj$smry) %in% c("matrix", "list") ){
      stop("'obj$smry' must be a matrix or list generated from one of the following functions:
           cace.meta.c, cace.meta.ic.")
    }
    if (class(obj$smry) == "matrix"){
      x <- obj$smry
      outcacei=x[substr(row.names(x), start=1, stop=4)=="cace", ]
    }
    else if (class(obj$smry) == "list"){
      outcacei <- obj$CACE
      if (nrow(outcacei) == 1) stop("forestplot cannot be made for a single study")
    }
    
    if(!nrow(data)==nrow(outcacei)) 
      stop("data and study-specific CACE have different lengths")
    
    data$miss.r0 <- ifelse((data$n000==0 & data$n001==0 & data$n010==0 & data$n011==0), 1, 0)
    data$miss.r1 <- ifelse((data$n100==0 & data$n101==0 & data$n110==0 & data$n111==0), 1, 0)
    data$miss <- ifelse((data$miss.r0==1|data$miss.r1==1), 1, 0)
    tmp1 <- data[order(data$miss.r0, data$miss.r1),]
    tmp2 <- cbind(tmp1, outcacei)
    tmp3 <- tmp2[order(tmp2$study.id),]

    study.name <- tmp3$study.name
    study.id <- tmp3$study.id
    cace_median <- round(as.vector(tmp3[, "50%"]), 3)
    cace_mean <- round(as.vector(tmp3[, "Mean"]), 3)
    cace_lower <- round(as.vector(tmp3[, "2.5%"]), 3)
    cace_upper <- round(as.vector(tmp3[, "97.5%"]), 3)

    values <- list(
        median  = c(NA, NA, cace_median, x["CACE", "50%"]), 
        mean  = c(NA, NA, cace_mean, x["CACE", "Mean"]), 
        lower = c(NA, NA, cace_lower, x["CACE", "2.5%"]), 
        upper = c(NA, NA, cace_upper, x["CACE", "97.5%"]) )
    
    tabletext<-cbind(
      c("", "Study", paste(study.name), "Overall"),
      
      c("", "CACE", paste(formatC(cace_mean, format='f', digits=3 ), " (",
                          formatC(cace_lower, format='f', digits=3 ),", ", 
                          formatC(cace_upper, format='f', digits=3 ), ")", sep="") , 
        paste(round(x["CACE", "Mean"], 3), " (", 
              round(x["CACE", "2.5%"], 3), ", ", 
              round(x["CACE", "97.5%"], 3), ")", sep="" ) )
    )
    
    if(obj$model %in% c("cace.single", "cace.meta.c")) {
      forestplot(tabletext, 
                 boxsize = .2, 
                 hrzl_lines = gpar(lwd=1, col="#444444"),
                 mean = values$mean,
                 lower = values$lower,
                 upper = values$upper,
                 new_page = TRUE, 
                 is.summary=c(TRUE, TRUE, rep(FALSE,length(outcacei)), TRUE),
                 clip=c(-1.0, 1.0),
                 col=fpColors(box="royalblue",line="darkblue", summary="royalblue")
      ) 
    }
    
    else if (obj$model=="cace.meta.ic") {
      
      missind <- c(NA, NA, tmp3$miss, 1)
      
      CACE.complete <- CACE.marginal <- values
      CACE.marginal$mean[missind==0] <- CACE.marginal$median[missind==0] <- 
        CACE.marginal$lower[missind==0] <- CACE.marginal$upper[missind==0] <- NA
      CACE.complete$mean[missind==1] <- CACE.complete$median[missind==1] <- 
        CACE.complete$lower[missind==1] <- CACE.complete$upper[missind==1] <- NA
      
      forestplot(tabletext, 
                 boxsize = .2, 
                 hrzl_lines = gpar(lwd=1, col="#444444"),
                 mean = cbind(CACE.complete$mean, CACE.marginal$mean),
                 lower = cbind(CACE.complete$lower, CACE.marginal$lower),
                 upper = cbind(CACE.complete$upper, CACE.marginal$upper),
                 new_page = TRUE, 
                 is.summary=c(TRUE, TRUE, rep(FALSE,nrow(outcacei)), TRUE, TRUE),
                 clip=c(-1.0, 1.0), 
                 lty.ci = c(1, 5),
                 col=fpColors(box="royalblue",line="darkblue", summary="royalblue")
      ) 
    }
    
    invisible()
}

