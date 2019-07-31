vcov.BinaryEPPM <-
function (object, model = c("full", "p", "scale.factor"), ...) {
    vc <- object$vcov
    k <- length(object$coefficients$p.est)
    m <- length(object$coefficients$scalef.est)
    match.arg(model)
    if (missing(model)) { model <- c("full") } 
# Checking for correct model option
      if ((model!="full") & (model!="p") & (model!="scale.factor")) {
         warning("\n","unknown model option")
         vc <- NULL
                                   } else {
         switch(model, full = { vc 
                  }, p = {
           vc <- vc[seq.int(from = 1, to = k, by = 1), 
              seq.int(from = 1, to = k, by = 1), drop = FALSE]
                  }, scale.factor = {
           if (m==0) { warning("\n","the scale-factor model has no elements")
              vc <- NULL
                     } else { 
              vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + 
                  k, drop = FALSE]
              colnames(vc) <- rownames(vc) <- names(object$coefficients$scalef.est)
              vc }} ) } # end of if model!="full", etc.,
    return(vc) }
