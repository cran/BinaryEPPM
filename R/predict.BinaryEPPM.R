predict.BinaryEPPM <-
function (object, newdata = NULL, type = c("response", "linear.predictor.p",
   "linear.predictor.scale.factor", "p", "scale.factor", "scale.factor.limits",
   "mean", "variance", "distribution", "distribution.parameters"), na.action = na.pass, ...) {

    type <- match.arg(type)

    if (missing(newdata)) {
       nobs <- nrow(object$covariates.matrix.p) 
       vzero <- rep(0,nobs)
       vone  <- rep(1,nobs)
       ntrials <- list(rep(c(0),nobs))
       ntrials <- lapply(1:nobs, function(i) {
          nmax <- object$vnmax[i]
          ntrials[[i]] <- c(rep(nmax,nmax+1)) } )
       output.model  <- Model.Binary(parameter=object$optim$par,
                          model.type=object$model.type,model.name=object$model.name,
                          link=object$link,ntrials=ntrials,
                          covariates.matrix.p=object$covariates.matrix.p,
                          covariates.matrix.scalef=object$covariates.matrix.scalef,
                          offset.p=object$offset.p,offset.scalef=object$offset.scalef)
        rval <- switch(type, response = {
            object$fitted.values
        }, linear.predictor.p = {
             vlp <- as.vector(object$covariates.matrix.p %*% object$coefficients$p.est
                      + object$offset.p)          
             vlp
        }, linear.predictor.scale.factor = {
             if (object$model.type=="p and scale-factor") { 
                vlp  <-as.vector(object$covariates.matrix.scalef %*% object$coefficients$scalef.est
                         + object$offset.scalef)
                                      } else { vlp <- NULL}
             vlp
        }, p = { vp <- rep(0,nobs)
             vp <- sapply(1:nobs, function(i) {
                          probability <- output.model$probabilities[[i]]
                          nmax        <- object$vnmax[i]
                          vp[i]       <- t(probability)%*%c(0:nmax) / nmax } )
             vp
        }, scale.factor = { scalef.prob <- rep(1,nobs)
                scalef.prob <- sapply(1:nobs, function(i) {
                              probability    <- output.model$probabilities[[i]]
                              nmax           <- object$vnmax[i]
                              fmean          <- t(probability)%*%c(0:nmax) 
                              fvariance      <- t(probability)%*%((c(0:nmax)-as.vector(fmean))^2) 
                              scalef.prob[i] <- fvariance / (fmean*(1 - fmean/nmax)) })
             scalef.prob
        }, scale.factor.limits = { 
# limits on scale factor    
             if (object$model.name=="binomial") { vlimits <- NULL }
             if (object$model.name=="EPPM extended binomial") {
                vlp <- as.vector(object$covariates.matrix.p %*% object$coefficients$p.est + object$offset.p)
# inverse of link function
                vp <- attr(object$link, which="p")$linkinv(vlp)
                vlimits <- data.frame(vzero, vone/(vone-vp)) }
             if (object$model.name=="beta binomial") {
                 vlimits <- data.frame(vone + output.model$Dparameters$lower.limit*(object$vnmax-1) / 
                                              (vone + output.model$Dparameters$lower.limit),
                                       object$vnmax) }
             if (object$model.name=="correlated binomial") {
                 vlimits <- data.frame(vone + output.model$Dparameters$lower.limit*(object$vnmax-1), 
                                       vone + output.model$Dparameters$upper.limit*(object$vnmax-1)) }
             names(vlimits) <- c("lower","upper")
             vlimits
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
        }, distribution = { 
             output.model$probabilities
        }, distribution.parameters = {  
             output.model$Dparameters } )
        return(rval)

     } else {

# Regardless of whether the original data is a list or a data frame, 
# newdata must be a data frame. 
        if (is.data.frame(newdata)==FALSE) { 
           warning("\n","Input newdata is not a data frame.")
           return(object=NULL)
                                    } else {
           if (object$data.type==FALSE) {
              if (("vnmax"%in%names(newdata))==FALSE) {
                 warning("\n","No variable named vnmax in newdata.")
              return(object=NULL) }} } # end of check data.frame

        mf <- model.frame(delete.response(object$terms[["p"]]), 
                          data = newdata, na.action = na.action, xlev = object$levels[["p"]])
        newdata <- newdata[rownames(mf), , drop = FALSE]
        offset <- list(p = rep.int(0, nrow(mf)), scale.factor = rep.int(0, nrow(mf)))

        X <- model.matrix(delete.response(object$terms[["p"]]), 
                   mf, contrasts = object$contrasts$p)
        if (!is.null(off.num <- attr(object$terms$p, "offset"))) {
            for (j in off.num) offset[[1L]] <- offset[[1L]] + 
              eval(attr(object$terms$p, "variables")[[j + 1L]], newdata) } # end of !is.null

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

          if (object$data.type==TRUE) {
# extracting name of numerator in response variable for use 
             nobs <- nrow(newdata) 
             ntrials <- list(rep(c(0),nobs))
             denom.name <- attr(object$vnmax, which="names")[1]
             vnmax <- newdata[as.character(denom.name)]
                            } else {
             nobs <- length(newdata[[1]])
             ntrials <- list(rep(c(0),nobs))
             vnmax <- newdata["vnmax"] } # end of if data.type
          ntrials <- lapply(1:nobs, function(i) { nmax <- vnmax[[1]][i]
                            ntrials[[i]] <- c(rep(nmax,nmax+1)) } ) # end of lapply
          vzero <- rep(0,nobs)
          vone  <- rep(1,nobs)
          output.model  <- Model.Binary(parameter=object$optim$par,
                             model.type=object$model.type,model.name=object$model.name,
                             link=object$link,ntrials=ntrials,
                             covariates.matrix.p=X,covariates.matrix.scalef=Z,
                             offset.p=offset[[1L]],offset.scalef=offset[[2L]])
        rval <- switch(type, response = {
             vlp <- drop(X %*% object$coefficients$p.est + offset[[1L]])
# inverse of link function
             vp <- attr(object$link, which="p")$linkinv(vlp)
             vp
        }, linear.predictor.p = {
             vlp <- drop(X %*% object$coefficients$p.est + offset[[1L]])
             vlp
        }, linear.predictor.scale.factor = {
             if (object$model.type=="p and scale-factor") {
                vlp  <-drop(Z %*% object$coefficients$scalef.est + offset[[2L]])
                                                   } else { vlp <- NULL }
             vlp
        }, p = { vp <- rep(0,nobs)
             vp <- sapply(1:nobs, function(i) {
                          probability <- output.model$probabilities[[i]]
                          nmax        <- length(ntrials[[i]]) - 1
                          vp[i]       <- t(probability)%*%c(0:nmax) / nmax } )
             vp
        }, scale.factor = { scalef.prob <- rep(1,nobs)
                scalef.prob <- sapply(1:nobs, function(i) {
                              probability    <- output.model$probabilities[[i]]
                              nmax           <- length(ntrials[[i]]) - 1
                              fmean          <- t(probability)%*%c(0:nmax) 
                              fvariance      <- t(probability)%*%((c(0:nmax)-as.vector(fmean))^2) 
                              scalef.prob[i] <- fvariance / (fmean*(1 - fmean/nmax)) })
             scalef.prob
        }, scale.factor.limits = { 
# limits on scale factor    
             if (object$model.name=="binomial") { vlimits <- NULL }
             if (object$model.name=="EPPM extended binomial") {
                vp <- rep(0,nobs)
                vp <- sapply(1:nobs, function(i) {
                          probability <- output.model$probabilities[[i]]
                          nmax        <- length(ntrials[[i]]) - 1
                          vp[i]       <- t(probability)%*%c(0:nmax) / nmax } )
                vlimits <- data.frame(vzero, vone/(vone-vp)) }
             if (object$model.name=="beta binomial") {
                 vlimits <- data.frame(vone + output.model$Dparameters$lower.limit*(vnmax-1) / 
                                            (vone + output.model$Dparameters$lower.limit),
                                       vnmax) }
             if (object$model.name=="correlated binomial") {
                 vlimits <- data.frame(vone + output.model$Dparameters$lower.limit*(vnmax-1), 
                                       vone + output.model$Dparameters$upper.limit*(vnmax-1)) }
             names(vlimits) <- c("lower","upper")
             vlimits
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
