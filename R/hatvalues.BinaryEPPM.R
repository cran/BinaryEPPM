hatvalues.BinaryEPPM <-
function (model, ...) {
    x <- model$covariates.matrix.p
    if (model$model.type=="p only") { z <- matrix(rep.int(1, NROW(x)), ncol=1)
                                } else {
       z <- model$covariates.matrix.scalef } # end if model.name
    if (NCOL(x) < 1L) 
        return(structure(rep.int(0, NROW(x)), .Names = rownames(x)))
    if (is.null(model$offset$p)) { offset.p <- rep(0, NROW(x)) }
    if (is.null(model$offset$scalef)) { offset.scalef <- rep(0, NROW(z)) }
    wts <- weights(model$model)
    if (is.null(wts)) { wts <- rep(1, NROW(x)) }
    beta <- model$coefficients$p.est   
    if (model$model.name=="binomial") { gamma <- matrix(rep.int(1, 1), ncol=1)
                                } else {
       gamma <- model$coefficients$scalef.est } # end if model.name
    eta <- as.vector(x %*% beta + offset.p)
    scale.factor_eta <- as.vector(z %*% gamma + offset.scalef)
    vone <- rep(1,length(eta))
# inverse of link function
    p <- attr(model$link, which="p")$linkinv(eta)
    if (model$model.name=="binomial") { scale.factor <- z
                               } else {
# modeling scalefactor
          scale.factor <- exp(scale.factor_eta) }
     vvariance <- p*(vone-p)*scale.factor/model$vnmax
# first differential of p w.r.t. eta (linear predictor)
    dmu <- attr(model$link, which="p")$mu.eta(eta)
# constructing a vector of cases from the list of groups
    if (model$data.type==TRUE) { 
          wk.p    <- p
          wk.wts  <- wts
          wk.vvariance <- vvariance
          wk.dmu <- dmu
          wk.x    <- x
                            } else {
       wks <- length(x[1,])
       nend <- 0
       for (ilist in 1:length(model$list.data)) { 
          for ( i in 1:length(model$list.data[[ilist]])) { 
             nt <- model$list.data[[ilist]][i] 
             if (nt>0) { 
                if (nend==0) {
                    wk.x <- as.vector(c(rep(x[ilist,],nt)))
                    wk.p <- c(rep(p[ilist],nt))
                    wk.wts <- c(rep(wts[ilist],nt))
                    wk.vvariance <- c(rep(vvariance[ilist],nt))
                    wk.dmu <- c(rep(dmu[ilist],nt))
                             } else {
                    wk.x <- append(wk.x, as.vector(c(rep(x[ilist,],nt))), after=(nend*wks))
                    wk.p <- append(wk.p, c(rep(p[ilist],nt)), after=nend)
                    wk.wts <- append(wk.wts, c(rep(wts[ilist],nt)), after=nend)
                    wk.vvariance <- append(wk.vvariance, c(rep(vvariance[ilist],nt)), after=nend)
                    wk.dmu <- append(wk.dmu, c(rep(dmu[ilist],nt)), after=nend) }
                nend <- nend + nt } } # end of for i loop 
                                 } # end of for ilist loop 
       wk.x <- t(matrix(wk.x, nrow=wks))
                                } # end of if model$data.type
    wkm <- sqrt(as.vector(wk.wts * wk.dmu * wk.dmu / wk.vvariance)) * wk.x
    xwx1 <- chol2inv(qr.R(qr(wkm)))
    as.vector(diag(wkm %*% xwx1 %*% t(wkm))) }
