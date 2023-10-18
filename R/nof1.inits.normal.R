#' @title nof1.inits.normal
#'
#' @description The nof1.inits.normal function generates initial parameter
#' values to be passed to the \code{jags.fit} function.
#'
#' @param data A list containing Y (the outcome), vectors of 0s and 1s in
#' the model matrix for each \emph{Treat.name}, indicating whether or not the
#' response was recorded under the treatment, and base spline information
#' @return A list of the initial values for each treatment-specific intercept
#' and spline term to be used in the \code{jags.fit} function
#' @export

nof1.inits.normal <- function(data,
                              Treat.name,
                              bs.trend,
                              bs_df,
                              n.chains){
  Y <- data$Y
  inits <- tryCatch({
    Treat.matrix <- NULL
    for(i in Treat.name){
      Treat.matrix <- cbind(Treat.matrix, data[[paste0("Treat_", i)]])
    }

    if (bs.trend){
      for (i in 1:bs_df){
        Treat.matrix <- cbind(Treat.matrix, data[[paste0("bs", i)]])
      }
    }

    model <- lm(Y ~ Treat.matrix - 1)
    co <- coef(summary(model))

    # Generate initial values
    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    for(i in 1:n.chains){

      for(j in 1:length(Treat.name)){
        initial.values[[i]][[paste0("beta_", Treat.name[j])]] <- co[j,1] +
          rnorm(1) * co[j,2]

        # Take care of bounded beta values if uniform prior prespecified
        if (beta.prior[[1]] == "dunif") {
          initial.values[[i]][[paste0("beta_", Treat.name[j])]] <-
            max(initial.values[[i]][[paste0("beta_",
                                            Treat.name[j])]],
                beta.prior[[2]])
          initial.values[[i]][[paste0("beta_", Treat.name[j])]] <-
            min(initial.values[[i]][[paste0("beta_", Treat.name[j])]],
                beta.prior[[3]])
        }
      }

      if(bs.trend){
        for(j in 1:bs_df){
          initial.values[[i]][[paste0("eta_", j)]] <-
            co[length(Treat.name)+j, 1] +
            rnorm(1, 0, co[length(Treat.name)+j, 2])
        }
      }
    }

    if(any(is.nan(unlist(initial.values)))) initial.values <- NULL

    #if initial value generated is too big (i.e. model didn't work because
    # of sparse data), just set initial value to be NULL
    if(!is.null(initial.values)){
      if(any(abs(unlist(initial.values)) > 100)) initial.values <- NULL
    }
    initial.values
  }, error = function(err){
    print(paste("Initial value not working: ",err))
    return(NULL)
  }, warning = function(warning){
    print(paste("Warning: ", warning))
    return(NULL)
  })
  return(inits)
}
