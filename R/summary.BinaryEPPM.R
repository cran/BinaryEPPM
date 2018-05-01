summary.BinaryEPPM <-
function(object, ...) {

       nobs <- nrow(object$covariates.matrix.p) 
       p.par        <- rep(0,nobs)
       scalef.par   <- rep(1,nobs)
       vone         <- rep(1,nobs)
       scalef.limit <- rep(0,nobs)
       npar   <- length(object$optim$par)
# Calculation of p and means from parameter estimates and design matrices
       lp.p <- object$covariates.matrix.p%*%object$coefficients$p.est + 
                    object$offset.p

# inverse of link function
       p.par <- attr(object$link, which="p")$linkinv(lp.p)

       if (object$model.type=="p only") { 
           wk.object <- object
           if (object$model.name=="generalized binomial") { 
# Changing to test of b=1 for the generalized binomial
# error in version 2.0 corrected in version 2.1
               wk.object$coefficients$scalef.est <- 
                    wk.object$coefficients$scalef.est - 1 } # end of generalized binomial 
           coeff.table.p <- coeftest(wk.object)
           coeff.table.scalef <- NULL 
                     } else {         
           wk.object <- object
           wk.object$coefficients <- object$coefficients[1]
           wk.object$vcov <- vcov(object,model="p")
           coeff.table.p <- coeftest(wk.object)
           wk.object <- object
           wk.object$coefficients <- object$coefficients[2]
           wk.object$vcov <- vcov(object,model="scale.factor")
           coeff.table.scalef <- coeftest(wk.object) 
                                         } # end if (model.type=="p only")

       object <- list(data.type=object$data.type, call=object$call, formula=object$formula, 
                    model.type=object$model.type, model.name=object$model.name, 
                    link=object$link, offset.p=object$offset.p,offset.scalef=object$offset.scalef,
                    coeff.table.p=coeff.table.p, coeff.table.scalef=coeff.table.scalef, loglik=object$loglik, 
                    n=object$nobs, nobs=object$nobs, df.null=object$df.null, df.residual=object$df.residual,
                    vnmax=object$vnmax, weights=object$weights, converged=object$converged, 
                    method=object$method,pseudo.r.squared=object$pseudo.r.squared,
                    optim=object$optim,control=object$control, fitted.values=object$fitted.values,
                    y=object$y, terms=object$terms, npar=npar) 

      class(object) <- "summaryBinaryEPPM"
      object }
