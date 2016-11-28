EPPMprob <-
function (vlambda) {
    vwork <- vlambda[1:(length(vlambda)-1)]
    tran  <- cbind(rbind(diag(-vlambda, length(vlambda)))) +
             cbind(0,rbind(diag(vwork, (length(vlambda)-1)),0))
    tran  <- expm(tran)
    prob  <- round(tran[1,], digits=20)
    return(prob)                }
