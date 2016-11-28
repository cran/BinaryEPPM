coef.BinaryEPPM <-
function(object, prtpar = c("full", "p", "scale.factor"), ...) {
      if (missing(prtpar)) { prtpar <- c("full") } 
# Checking for correct prtpar option
      if ((prtpar!="full") & (prtpar!="p") & (prtpar!="scale.factor")) {
         cat("\n","unknown prtpar option","\n")
         coefficients <- NULL
                                   } else {
         if (prtpar=="full") { 
            if (object$model.name=="binomial") {
               coefficients <- object$coefficients$p.est
                                        } else {
               coefficients <- c(object$coefficients$p.est,
                                 object$coefficients$scalef.est)
                                               } # end of if binomial
             } else { npar.p <- ncol(object$covariates.matrix.p)
               if (prtpar=="p") { coefficients <- object$coefficients$p.est
                                }
               if (prtpar=="scale.factor") { 
                  if (is.null(object$coefficients$scalef.est)==TRUE) {
                     coefficients <- NULL
                                           } else {
                     coefficients <- object$coefficients$scalef.est } } # end of if scale.factor
                                    } } # end of if prtpar!= full, etc.
      return(coefficients) }
