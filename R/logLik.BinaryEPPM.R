logLik.BinaryEPPM <-
function (object, ...) {
      structure(object$loglik, df = length(object$optim$par),
                nobs = nrow(object$covariates.matrix.p),
                class = "logLik") }
