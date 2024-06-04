doubrecip <-
function()
{
   linkfun <- function(mu) { ifelse( (mu<0.5), 1 - 0.5/mu,  0.5/(1-mu) - 1 ) }  
   linkinv <- function(eta) { (1+sign(eta))/2 - sign(eta)*(1/(1+abs(eta)))/2 }
   mu.eta <- function(eta) { (1/(1+abs(eta))^2)/2  }
   valideta <- function(eta) TRUE
   link <- paste0("doubrecip(",")")
   structure (list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
