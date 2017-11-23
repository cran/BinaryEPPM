negcomplog <-
function()
{ 
   linkfun <- function(mu) { -log(1 - mu) }
   linkinv <- function(eta) { pmax(pmin(1-exp(-eta), 
                              1-.Machine$double.eps), .Machine$double.eps) }
   mu.eta <- function(eta) { eta <- pmin(eta, 700)
                             pmax(exp(-eta), .Machine$double.eps) }
   valideta <- function(eta) TRUE
   link <- paste0("negcomplog(",")")
   structure (list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
