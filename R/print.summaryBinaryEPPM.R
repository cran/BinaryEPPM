print.summaryBinaryEPPM <-
function(x, ...) {

       if (x$data.type==TRUE) {
          cat("\n","Dependent variable a vector of numerator / denominator.","\n")
                        } else {
          cat("\n","Dependent variable is a list of binomial frequency distributions","\n")
                               } # end of data.type==

       if (is.null(x$converged)==FALSE) {
          cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 
             0.85)), "", sep = "\n")
          cat("Model type        :",x$model.type,"\n")
          cat("Model name        :",x$model.name,"\n")
          if (x$link=="powerlogit") { 
             cat("Link p            :",x$link,"power",attr(x$link, which="power"),"\n",sep=" ")
                                   } else {
             cat("Link p            :",x$link,"\n")
                                   } # end if link
          if (x$model.type=="p and scale-factor") { 
                         cat("Link scale-factor : log","\n") }
          offsetid.p      <- sum(x$offset.p)
          offsetid.scalef <- sum(x$offset.scalef)
          if ((offsetid.p!=0) | (offsetid.scalef!=0)) { 
              cat("non zero offsets in linear predictors","\n") }  


          if (x$link=="powerlogit") { 
              cat(paste("Coefficients (model for p with", x$link, 
                  "link power",attr(x$link, which="power"),"):\n", sep = " "))
                                           } else {
              cat(paste("Coefficients (model for p with", x$link, 
                  " link):\n", sep = " ")) } # end if link
          if ((x$model.type=="p only") & (x$model.name=="generalized binomial")) {
             wks <- length(x$optim$par) 
             cat(paste("Coefficient of",names(x$optim$par)[wks],"has 1 subtracted from it\n", sep=" "))
             cat(paste("so the test is against 1 i.e., a binomial.\n")) }
          print(x$coeff.table.p)      
          if (is.null(x$coeff.table.scalef)==FALSE) {
             if (x$model.type=="p and scale-factor") {
                cat(paste("Coefficients (model for scale-factor with log link)\n"))
                print(x$coeff.table.scalef) } } # end if is.null

          if (is.null(x$weights)==FALSE) {
             cat("\n","Maximum weighted likelihood regression.")
             if (x$data.type==TRUE) {
                cat("\n","Vector of weights used.","\n")
                             } else {
                cat("\n","List of weights used.","\n") }
             if (is.null(attr(x$weights, which="normalize"))==FALSE) {
                if (attr(x$weights, which="normalize")==TRUE) {
                   cat("Normalization to a value of",
                    attr(x$weights, which="norm.to.n"),".\n", sep = " ") }}
                                    } # end of is.null(weights)

          if (is.na(x$loglik)==TRUE) { cat("Log-likelihood is NA","\n")
                                          } else { 
             cat("\n","Type of estimator: ML (maximum likelihood)")
             cat("\n","Log-likelihood:",x$loglik,"on",length(x$optim$par),"Df", sep=" ")
             cat("\n","Pseudo R-squared:",x$pseudo.r.squared, "type", 
                       attr(x$pseudo.r.squared, which="names"), sep=" ")
             if (length(x$optim$par)==1) { 
                cat("\n","Single parameter binomial so no iteration", sep=" ")
                                                    } else {
                optim.method <- x$method
                if (optim.method=="Nelder-Mead") { 
                   cat("\n","Number of iterations:",x$optim$counts[1],"of optim method",optim.method, sep=" ") 
                                         } else { gradient.method <- attr(x$method,which="grad.method")
                   cat("\n","Number of iterations:",x$optim$counts[1],"of optim method",optim.method, 
                       "gradient method",gradient.method,sep=" ") }
                                                           } # end of if length(x$optim$par)=1
             code <- list(c("successful"),
                          c("iteration limit max has been reached"),
                          c(" "),c(" "),c(" "),c(" "),
                          c(" "),c(" "),c(" "),c(" "),
                          c("degeneracy of the Nelder-Mead"))
             wks <- attr(x$converged, which="code") + 1
             cat("\n","return code",attr(x$converged, which="code"),code[[as.numeric(wks)]],"\n", sep=" ")     
                 } # end of if is.na(x$loglik

                  } else { 
                       cat("\n","Failure of checks on entry arguments to BinaryEPPM")
                       cat("\n","or numerical derivative calculations failed.")
                         } # end of if (is.null(converged)

         }
