doubrecip <-
function()
{ 
   linkfun <- function(mu) { leng <- length(mu)
                 sapply(1:leng, function(i) {
                    if (mu[i]<0.5) { 1 - 0.5/mu[i]
                            } else { 0.5/(1-mu[i]) - 1 }} ) }
   linkinv <- function(eta) { (1+sign(eta))/2 - sign(eta)*(1/(1+abs(eta)))/2 }
   mu.eta <- function(eta) { (1/(1+abs(eta))^2)/2  }
   valideta <- function(eta) TRUE
   link <- paste0("doubrecip(",")")
   structure (list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
