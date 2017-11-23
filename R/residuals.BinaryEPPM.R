residuals.BinaryEPPM <-
function (object, type = c("spearson", "deviance", "pearson", 
    "response", "likelihood", "sdeviance"), ...) {
    type <- match.arg(type)
    p <- object$fitted.values
    nobs <- length(p) 
    vone <- rep(1, nobs) 
    vnmax <- as.vector(object$vnmax) 
    vnmax1 <- vnmax + vone
    wts <- weights(object$model)
    if (object$data.type==TRUE) { 
       y <- if (is.null(object$y)) { model.response(model.frame(object)) 
                                   } else { object$y }
                            } else {
                       y <- object$list.data 
                                   } # end of data.type=0
    if (is.null(wts)) { wts <- vone }
    scale.factor <- predict(object, type = "scale.factor")
# phi is the phi of the beta distribution as in beta regression
    phi <- vnmax / scale.factor - 1   
    distribution.parameters <- predict(object, type = "distribution.parameters")
    if (object$model.name=="binomial") { distribution.parameters[2] <- vone }

    if (object$data.type==TRUE) { 
       wk.y     <- y
       wk.p     <- p
       wk.scale.factor <- scale.factor
       wk.phi   <- phi
       wk.vnmax <- vnmax
       wk.wts   <- wts
       wk.vone  <- vone
       wk.resid <- wk.y - wk.p
       total.ninlist <- nobs

                            } else { 

       vninlist <- c(rep(0,length(object$list.data)))
       vninlist <- sapply(1:length(object$list.data), function(ilist) 
                        vninlist[ilist] <- sum(object$list.data[[ilist]]) )
       total.ninlist <- sum(vninlist)
       wk.p     <- c(rep(0,total.ninlist))
       wk.scale.factor <- c(rep(1,total.ninlist))
       wk.y     <- c(rep(0,total.ninlist))
       wk.phi   <- c(rep(0,total.ninlist))
       wk.vnmax <- c(rep(0,total.ninlist))
       wk.wts   <- c(rep(1,total.ninlist))
       wk.vone  <- c(rep(1,total.ninlist))
       wk.resid <- c(rep(0,total.ninlist))
       nstart <- 1
       nend <- 0
       for (ilist in 1:length(object$list.data)) { 
          ninlist <- sum(object$list.data[[ilist]])
          for ( i in 1:length(object$list.data[[ilist]])) { 
             nt <- object$list.data[[ilist]][i] 
             if (nt>0) {
                nend <- nend + nt
                wk.p[nstart:nend]     <- p[ilist]
                wk.scale.factor[nstart:nend] <- scale.factor[ilist]
                wk.phi[nstart:nend]   <- phi[ilist]
                wk.vnmax[nstart:nend] <- vnmax[ilist]
                wk.y[nstart:nend]     <- (i - 1) / vnmax[ilist] 
                wk.wts[nstart:nend]   <- as.vector(wts[ilist])
                wk.resid[nstart:nend] <- wk.y[nstart:nend] - wk.p[nstart:nend]
                nstart <- nstart + nt } } # end of for i loop
                                } # end of for ilist loop 
                                   } # end of if object$data.type  

    wk.resp.resid <- wk.resid
    if ((type=="pearson") | (type=="spearson") | (type=="likelihood")) {
       wk.resid <- wk.resid * sqrt( (wk.wts) / ( wk.p * (wk.vone - wk.p) * wk.scale.factor / wk.vnmax)) }
    if ((type=="deviance") | (type=="sdeviance") | (type=="likelihood")) {
        if ((object$model.name=="binomial") | (object$model.name=="generalized binomial")) {
            wk.va <- - log(wk.vone - wk.y) } # end of if object$model.name
        wk.y <- sapply(1:total.ninlist, function(i) { if ((wk.y[i]==0) | (wk.y[i]==1)) { 
                     wk.y[i] <- (wk.y[i]*(wk.vnmax[i]-1) + 0.5) / wk.vnmax[i]
                     } else { wk.y[i] <- wk.y[i] }} ) # end of sapply
        ll.obs <- c(rep(0,total.ninlist))
        ll.fit <- c(rep(0,total.ninlist))

        if (object$data.type==TRUE) { 
           ll.obs <- sapply(1:nobs, function(i) { 
              if ((object$model.name=="binomial") | (object$model.name=="generalized binomial")) {
                 probability <- GBprob(twoparameter=c(wk.va[i], distribution.parameters[[2]][i]), wk.vnmax[i])
                                                    } # end of if binomial or gen binomial
              if (object$model.name=="beta binomial") {
                 output.prob <- BBprob(twoparameter=c(wk.y[i], distribution.parameters[[2]][i]), wk.vnmax[i])
                 probability <- output.prob$probability
                                              } # end of beta binomial
              if (object$model.name=="correlated binomial") {
                 output.prob <- CBprob(twoparameter=c(wk.y[i], distribution.parameters[[2]][i]), wk.vnmax[i])
                 probability <- output.prob$probability
                                                    } # end of correlated binomial
                 wks <- sum(probability*as.numeric(object$list.data[[i]]>0))
                 if (wks>0) { ll.obs[i] <- log(wks)
                     } else { ll.obs[i] <- 0 }       
#	         ll.obs[i] <- log(sum(probability*as.numeric(object$list.data[[i]]>0))) 
                                   } ) # end of sapply
           ll.fit <- sapply(1:nobs, function(i) { 
              if ((object$model.name=="binomial") | (object$model.name=="generalized binomial")) {
                 probability <- GBprob(twoparameter=c(distribution.parameters[[1]][i],
                                                   distribution.parameters[[2]][i]), vnmax[i])
                                                    } # end of if binomial or gen binomial
              if (object$model.name=="beta binomial") {
                 output.prob <- BBprob(twoparameter=c(distribution.parameters[[1]][i],
                                                   distribution.parameters[[2]][i]), vnmax[i])
                 probability <- output.prob$probability
                                                    } # end of beta binomial
              if (object$model.name=="correlated binomial") {
                 output.prob <- CBprob(twoparameter=c(distribution.parameters[[1]][i],
                                                   distribution.parameters[[2]][i]), vnmax[i])
                 probability <- output.prob$probability
                                                    } # end of correlated binomial
                 wks <- sum(probability*as.numeric(object$list.data[[i]]>0))
                 if (wks>0) { ll.fit[i] <- log(wks)
                     } else { ll.fit[i] <- 0 }       
                                        } ) # end of sapply 
                             } else {
           wk.ind   <- c(rep(0,total.ninlist))
           wk.dp1   <- c(rep(0,total.ninlist))
           wk.dp2   <- c(rep(0,total.ninlist))
           nstart <- 1
           nend <- 0
           for (ilist in 1:length(object$list.data)) { 
              ninlist <- sum(object$list.data[[ilist]])
              for ( i in 1:length(object$list.data[[ilist]])) { 
                 nt <- object$list.data[[ilist]][i] 
                 if (nt>0) {
                    nend <- nend + nt
                    wk.ind[nstart:nend]   <- i  
                    wk.dp1[nstart:nend]   <- distribution.parameters[[1]][ilist]
                    wk.dp2[nstart:nend]   <- distribution.parameters[[2]][ilist]
                    nstart <- nstart + nt } } # end of for i loop
                                     } # end of for ilist loop 
        ll.obs <- sapply(1:total.ninlist, function(i) { 
           if ((object$model.name=="binomial") | (object$model.name=="generalized binomial")) {
              probability <- GBprob(twoparameter=c(wk.va[i], wk.dp2[i]), wk.vnmax[i])
                                                 } # end of if binomial or gen binomial
           if (object$model.name=="beta binomial") {
              output.prob <- BBprob(twoparameter=c(wk.y[i], wk.dp2[i]), wk.vnmax[i])
              probability <- output.prob$probability } # end of beta binomial
           if (object$model.name=="correlated binomial") {
              output.prob <- CBprob(twoparameter=c(wk.y[i], wk.dp2[i]), wk.vnmax[i])
              probability <- output.prob$probability } # end of correlated binomial
	        ll.obs[i] <- log(probability[wk.ind[i]]) } ) # end of sapply ll.obs
        ll.fit <- sapply(1:total.ninlist, function(i) { 
           if ((object$model.name=="binomial") | (object$model.name=="generalized binomial")) {
              probability <- GBprob(twoparameter=c(wk.dp1[i], wk.dp2[i]), wk.vnmax[i])
                                                  } # end of if binomial or gen binomial
           if (object$model.name=="beta binomial") {
              output.prob <- BBprob(twoparameter=c(wk.dp1[i], wk.dp2[i]), wk.vnmax[i])
              probability <- output.prob$probability } # end of beta binomial
           if (object$model.name=="correlated binomial") {
              output.prob <- CBprob(twoparameter=c(wk.dp1[i], wk.dp2[i]), wk.vnmax[i])
              probability <- output.prob$probability } # end of correlated binomial
	   ll.fit[i] <- log(probability[wk.ind[i]]) } ) # end of sapply ll.fit
                                   } # end of if object$data.type
        wk.resid.dev <- sqrt(wk.wts)*sign(wk.resid)*sqrt(2*abs(ll.obs - ll.fit)) 
                                   } # end of if deviance or sdeviance

    res <- switch(type, pearson = {
        wk.resid
    }, response = {
        wk.resp.resid
    }, deviance = {
        wk.resid.dev
    }, likelihood = {
# likelihood residuals
          sqrt(wk.wts)*sign(wk.resid.dev)*sqrt( hatvalues(object)*wk.resid*wk.resid  
                    + wk.resid.dev*wk.resid.dev ) 
    }, sdeviance = {
          wkv <- rep(0,length(wk.resid.dev)) 
# deviance residuals standardized to have an asymptotic variance of 1
          wkv <- sapply(1:length(wk.resid.dev), function(i) { 
              if (wk.scale.factor[i]<=0) { wkv[i] <- wkv[i]
                              } else { wkv[i] <-
                  wk.resid.dev[i] / sqrt( wk.scale.factor[i] * (1 - hatvalues(object)[i]) )
                                     } } ) 
          wkv
    }, spearson = {
# Pearson residuals standardized to have an asymptotic variance of 1
          wk.resid / sqrt( wk.vone - hatvalues(object) )
    }) # end of switch
    return(res) }
