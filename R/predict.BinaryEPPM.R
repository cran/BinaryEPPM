predict.BinaryEPPM <-
function (object, newdata = NULL, type = c("response", "linear.predictor", "p",
                      "scale.factor", "mean", "variance", "distribution", "distribution.parameters"),
                      na.action = na.pass, ...) {

    type <- match.arg(type)

    if (missing(newdata)) {
       nobs <- nrow(object$covariates.matrix.p) 
       ntrials <- list(rep(c(0),nobs))
       ntrials <- lapply(1:nobs, function(i) {
          nmax <- object$vnmax[i]
          ntrials[[i]] <- c(rep(nmax,nmax+1)) } )
       output.model  <- Model.Binary(parameter=object$optim$par,
                          model.type=object$model.type,model.name=object$model.name,
                          link=object$link,ntrials=ntrials,
                          covariates.matrix.p=object$covariates.matrix.p,
                          covariates.matrix.scalef=object$covariates.matrix.scalef,
                          offset.p=object$offset.p,
                          offset.scalef=object$offset.scalef)

        rval <- switch(type, response = {
            object$fitted.values
        }, linear.predictor = {
             vlp <- as.vector(object$covariates.matrix.p %*% object$coefficients$p.est + object$offset.p)
             vlp
        }, p = {
             vlp <- as.vector(object$covariates.matrix.p %*% object$coefficients$p.est + object$offset.p)
             vone <- rep(1,length(vlp))
# inverse of link function
             vp <- attr(object$link, which="p")$linkinv(vlp)
             vp
        }, scale.factor = { scalef.prob <- rep(1,nobs)
             if (object$model.type=="p only") {
                scalef.prob <- sapply(1:nobs, function(i) {
                              probability    <- output.model$probabilities[[i]]
                              nmax           <- object$vnmax[i]
                              fmean          <- t(probability)%*%c(0:nmax) 
                              fvariance      <- t(probability)%*%((c(0:nmax)-as.vector(fmean))^2) 
                              scalef.prob[i] <- fvariance / (fmean*(1 - fmean/nmax)) })}
             if (object$model.type=="p and scale-factor") {
                scalef.prob <- exp(as.vector(object$covariates.matrix.scalef %*% 
                                   object$coefficients$scalef.est + object$offset.scalef)) }
             scalef.prob
        }, mean = { mean.prob <- rep(0,nobs) 
             mean.prob <- sapply(1:nobs, function(i) {
                              probability  <- output.model$probabilities[[i]]
                              nmax         <- object$vnmax[i]
                              mean.prob[i] <- t(probability)%*%c(0:nmax) } )
             mean.prob
        }, variance = { variance.prob <- rep(0,nobs) 
             variance.prob <- sapply(1:nobs, function(i) {
                              probability      <- output.model$probabilities[[i]]
                              nmax             <- object$vnmax[i]
                              fmean            <- t(probability)%*%c(0:nmax) 
                              variance.prob[i] <- t(probability)%*%((c(0:nmax)-as.vector(fmean))^2) } )
             variance.prob  
        }, distribution = { probabilities <- output.model$probabilities
             probabilities  
        }, distribution.parameters = {  
             output.model$Dparameters } )
        return(rval)

     } else {

        mf <- model.frame(delete.response(object$terms[["p"]]), 
                          data = newdata, na.action = na.action, xlev = object$levels[["p"]])
        newdata <- newdata[rownames(mf), , drop = FALSE]
        offset <- list(p = rep.int(0, nrow(mf)), scale.factor = rep.int(0, nrow(mf)))
        X <- model.matrix(delete.response(object$terms[["p"]]), 
                   mf, contrasts = object$contrasts$p)
        if (!is.null(object$call$offset)) 
            offset[[1L]] <- offset[[1L]] + eval(object$call$offset, newdata)
        if (!is.null(off.num <- attr(object$terms$p, "offset"))) {
            for (j in off.num) offset[[1L]] <- offset[[1L]] + 
              eval(attr(object$terms$mean, "variables")[[j + 1L]], newdata) }

        Z <- matrix(c(rep(1,nrow(newdata))),ncol=1)
        if (object$model.type=="p and scale-factor") {
            mf <- model.frame(delete.response(object$terms[["scale.factor"]]), 
                      data = newdata, na.action = na.action, xlev = object$levels[["scale.factor"]])
            Z <- model.matrix(delete.response(object$terms$scale.factor),
                       mf, contrasts = object$contrasts[["scale.factor"]])
            if (!is.null(off.num <- attr(object$terms$scale.factor, "offset"))) {
                for (j in off.num) offset[[2L]] <- offset[[2L]] + 
                  eval(attr(object$terms$scale.factor, "variables")[[j + 1L]], newdata) }
                                                     } # end of if p and scale-factor

          nobs <- nrow(newdata) 
          ntrials <- list(rep(c(0),nobs))
          if (object$data.type==TRUE) {
# extracting name of numerator in response variable for use 
             denom.name <- attr(object$vnmax, which="names")[1]
             vnmax <- newdata[as.character(denom.name)]
             ntrials <- lapply(1:nobs, function(i) { nmax <- vnmax[[1]][i]
                               ntrials[[i]] <- c(rep(nmax,nmax+1)) } ) # end of lapply
                            } else {
             ntrials <- lapply(1:nobs, function(i) { nmax1 <- object$vnmax[i] + 1
                   ntrials[[i]] <- c(rep(object$vnmax[i],nmax1)) } ) # end of lapply
                                  } # end of if data.type

          output.model  <- Model.Binary(parameter=object$optim$par,
                             model.type=object$model.type,model.name=object$model.name,
                             link=object$link,ntrials=ntrials,
                             covariates.matrix.p=X,covariates.matrix.scalef=Z,
                             offset.p=offset[[1L]],offset.scalef=offset[[2L]])

        rval <- switch(type, response = {
             vlp <- drop(X %*% object$coefficients$p.est + offset[[1L]])
             vone <- rep(1,length(vlp))
# inverse of link function
             vp <- attr(object$link, which="p")$linkinv(vlp)
             vp
        }, linear.predictor = {
             vlp <- drop(X %*% object$coefficients$p.est + offset[[1L]])
             vlp
        }, p = {
             vlp <- drop(X %*% object$coefficients$p.est + offset[[1L]])
             vone <- rep(1,length(vlp))
# inverse of link function
             vp <- attr(object$link, which="p")$linkinv(vlp)
             vp
        }, scale.factor = { scalef.prob <- rep(1,nobs)
             if (object$model.type=="p only") {
                scalef.prob <- sapply(1:nobs, function(i) {
                              probability    <- output.model$probabilities[[i]]
                              nmax           <- length(ntrials[[i]]) - 1
                              fmean          <- t(probability)%*%c(0:nmax) 
                              fvariance      <- t(probability)%*%((c(0:nmax)-as.vector(fmean))^2) 
                              scalef.prob[i] <- fvariance / (fmean*(1 - fmean/nmax)) })}
             if (object$model.type=="p and scale-factor") {
                veta <- drop(Z %*% object$coefficients$scalef.est + offset[[2L]])
                scalef.prob <- exp(veta) }
             scalef.prob
        }, mean = { mean.prob <- rep(0,nobs) 
             mean.prob <- sapply(1:nobs, function(i) {
                              probability  <- output.model$probabilities[[i]]
                              nmax           <- length(ntrials[[i]]) - 1
                              mean.prob[i] <- t(probability)%*%c(0:nmax) } )
             mean.prob
        }, variance = { variance.prob <- rep(0,nobs) 
             variance.prob <- sapply(1:nobs, function(i) {
                              probability      <- output.model$probabilities[[i]]
                              nmax           <- length(ntrials[[i]]) - 1
                              fmean            <- t(probability)%*%c(0:nmax) 
                              variance.prob[i] <- t(probability)%*%((c(0:nmax)-as.vector(fmean))^2) } )
             variance.prob  
        }, distribution = { probabilities <- output.model$probabilities
             probabilities  
        }, distribution.parameters = {  
             output.model$Dparameters } )
        return(rval) } }
