#' @title find.max.gelman
#'
#' @description Calculates the maximum Gelman and Rubin's convergence diagnostic
#' for a sample generated from \code{coda.samples}.
#'
#' @param x A list of samples created from the \code{coda.samples} function.
#' @return Returns the maximum Gelman and Rubin's convergence diagnostic for the
#' samples.
#' @export

find.max.gelman <- function(samples){
  samples <- lapply(samples, function(x) { x[,colSums(abs(x)) != 0] })
  gr.uni.stats <- gelman.diag(samples, multivariate = FALSE)$psrf[,1]
  largest.gr <- max(gr.uni.stats[!is.nan(gr.uni.stats)])
  out = list(largest.gr)
  names(out) <- names(which(gr.uni.stats == largest.gr))
  return(out)
}
