#' @title new.mcmc
#'
#' @description Converts an object from \code{coda.samples} to an
#' \code{mcmc.list} object complete with the names of the saved parameters.
#'
#' @param x A list of samples created from the \code{coda.samples} function.
#' @return An \code{mcmc.list} object containing the samples
#' @export

new.mcmc <- function(x){
  n.chains <- length(x)
  n.var <- coda::nvar(x)
  newobjects <- vector("list", length = n.chains)

  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var,
                              dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- x[[i]]
    newobjects[[i]] <- mcmc(newobjects[[i]])
  }
  mcmc.list(newobjects)
}
