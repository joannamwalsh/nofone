#' @title indiv.summary
#'
#' @description Creates density plots for each parameter and a summary table
#' results.table creates an overall data frame of organized results for everyone

#' Creates kernel density plots and summary tables for each parameter from an
#' individual analysis
#'
#' @param samples MCMC samples of each parameter
#' @return A list containing a list of kernel density plots for each parameter
#' and a summary table of the parameter distributions
#' @import ggplot2
#' @import stats
#' @export

indiv.summary <- function(samples) {
  #create kernel density plots of each parameter
  sampledata <- as.data.frame(samples)
  plotlist <- list()
  for (col in 1:ncol(sampledata)) {
    col.name <- names(sampledata)[col]
    plotlist[[col]] <- ggplot(data = sampledata) +
      geom_density(aes(x = .data[[col.name]])) +
      geom_vline(aes(xintercept = mean(.data[[col.name]])),
                 linetype = "dashed") +
      theme_bw()
  }

  #create a parameter summary table
  extraparamtable <- data.frame(value = names(sampledata))
  for (row in 1:nrow(extraparamtable)) {
    extraparamtable$mean[row] <-
      mean(samples[, extraparamtable$value[row]])
    extraparamtable$median[row] <-
      median(samples[, extraparamtable$value[row]])
    extraparamtable$ci2.5[row] <-
      quantile(samples[, extraparamtable$value[row]], 0.025)
    extraparamtable$ci50[row] <-
      quantile(samples[, extraparamtable$value[row]], 0.5)
    extraparamtable$ci97.5[row] <-
      quantile(samples[, extraparamtable$value[row]], 0.975)
  }
  return(list(plotlist, extraparamtable))
}

#' @title results.table
#'
#' @description Computes and presents a summary of the results from an
#' individual analysis.
#'
#' @param indiv ID of the individual
#' @param treatment The name of the treatment variable in the data set. This
#' variable should have the same length as the outcome variable and can be
#' character or numeric.
#' @param groupvars The names of grouping variables that the user would like to
#' be included in the summary table (not including the stratifying variables)
#' @param stratvars The names of variables that the user would like the plots to
#'  be faceted by. A separate plot will be created for each level of the first
#' stratifying variable, and these plots will be arranged next to each other. If
#' a second stratifying variable is specified, this variable will be the facet
#' on the x axis for each of the plots.
#' @param data A dataframe containing all of the N of 1 data
#' @param data.ind  A data frame containing all of the individual's N of 1 data
#' @param Y The name of the outcome variable in the dataset. This variable
#' should be a vector with NA's included in time order.
#' @param nof1 Modeling result of class \code{nof1.data} produced by
#' \code{nof1.data}
#' @param samples MCMC samples of each parameter
#' @param clinicaldiff The lowest possible practical clinical difference between
#' treatment coefficients
#' @param alpha The desired alpha value for confidence intervals. The default is
#' 0.05.
#' @return Returns a list of summary statistics of the raw data and the
#' fitted model for each treatment comparison for the individual.
#' \item{raw.y.mean}{The raw mean of the outcome for each treatment}
#' \item{raw.y.median}{The raw median of the outcome for each treatment}
#' \item{post.coef.mean}{The posterior mean of the coefficient for each
#' treatment}
#' \item{post.coef.median}{The posterior median of the coefficient for each
#' treatment}
#' \item{post.y.mean}{The posterior mean of the outcome for each treatment}
#' \item{post.y.median}{The posterior median of the outcome for each treatment}
#' \item{post.coef.ci}{The credible interval of the coefficient for each
#' treatment}
#' \item{post.y.ci}{The credible interval of the outcome for each treatment}
#' \item{comp.treat.post.coef}{The posterior quantiles of one coefficient minus
#' the other when comparing two treatments}
#' \item{p.comp.coef.greater.0}{The posterior probability of one coefficient
#' minus the other is greater than 0}
#' \item{comp.treat.post.y}{The posterior quantiles of one outcome minus the
#' other when comparing two treatments}
#' @import dplyr
#' @import utils
#' @import stats
#' @export

results.table <- function(indiv,
                          treatment,
                          groupvars,
                          stratvars,
                          data,
                          data.ind,
                          Y,
                          nof1,
                          samples,
                          clinicaldiff,
                          alpha = 0.05) {

  Treat.name =  as.vector(sort(unique(
    data[[treatment]])))

  n            <- NULL
  raw.y.mean   <- NULL
  raw.y.median <- NULL
  raw.y.sd     <- NULL
  post.coef.mean   <- NULL
  post.coef.median <- NULL
  post.y.mean   <- NULL
  post.y.median <- NULL
  post.coef.ci <- NULL
  post.y.ci    <- NULL

  for (i in Treat.name){

    if (i %in% data.ind[[treatment]]) {
      n <- c(n, sum(!is.na(nof1$data$Y[nof1$Treat == i])))
      raw.y.mean   <- c(raw.y.mean, mean(nof1$data$Y[nof1$Treat == i],
                                         na.rm = TRUE))
      raw.y.median <- c(raw.y.median, median(nof1$data$Y[nof1$Treat == i],
                                             na.rm = TRUE))
      raw.y.sd     <- c(raw.y.sd, sd(nof1$data$Y[nof1$Treat == i],
                                     na.rm = TRUE))
    } else {
      n <- c(n, NA)
      raw.y.mean <- c(raw.y.mean, NA)
      raw.y.median <- c(raw.y.median, NA)
      raw.y.sd <- c(raw.y.sd, NA)
    }


    col.treat.name <- paste0("beta_", i)

    if (i %in% data.ind[[treatment]]) {
      post.coef <- samples[, col.treat.name]
    } else {
      post.coef <- NA
    }

    post.coef.mean   <- c(post.coef.mean, mean(post.coef))
    post.coef.median <- c(post.coef.median, median(post.coef))

    post.y <- post.coef
    post.y.mean   <- c(post.y.mean, mean(post.y))
    post.y.median <- c(post.y.median, median(post.y))

    post.coef.ci <- rbind(post.coef.ci, quantile(post.coef,
                                                 c(alpha/2, 1-alpha/2),
                                                 na.rm = T))
    post.y.ci    <- rbind(post.y.ci, quantile(post.y, c(alpha/2, 1-alpha/2),
                                              na.rm = T))
  }

  names(n) <-
    names(raw.y.mean) <-
    names(raw.y.median) <-
    names(raw.y.sd) <-
    names(post.coef.mean) <-
    names(post.coef.median) <-
    names(post.y.mean) <-
    names(post.y.median) <-
    Treat.name

  rownames(post.coef.ci) <- rownames(post.y.ci) <- Treat.name

  comp.treat.name <- t(utils::combn(Treat.name, 2))
  comp.treat.post.coef   <- NULL
  p.comp.coef.greater.0 <- NULL
  for (i in 1:nrow(comp.treat.name)){

    col.treat.name.1 <- paste0("beta_", comp.treat.name[i, 1])
    col.treat.name.2 <- paste0("beta_", comp.treat.name[i, 2])
    if (comp.treat.name[i, 1] %in% data.ind[[treatment]]) {
      post.coef.1 <- samples[, col.treat.name.1]
    } else {
      post.coef.1 <- NA
    }

    if (comp.treat.name[i, 2] %in% data.ind[[treatment]]) {
      post.coef.2 <- samples[, col.treat.name.2]
    } else {
      post.coef.2 <- NA
    }

    comp.treat.post.coef <- rbind(comp.treat.post.coef,
                                  quantile(post.coef.2 - post.coef.1,
                                           c(alpha/2, 0.5, 1-alpha/2),
                                           na.rm = T))

    p.comp.coef.greater.0 <- c(p.comp.coef.greater.0,
                               mean((post.coef.2 - post.coef.1) > 0))
  }

  rownames(comp.treat.post.coef) <- paste(comp.treat.name[, 2],
                                          comp.treat.name[, 1],
                                          sep = "_minus_")
  names(p.comp.coef.greater.0)  <- paste(comp.treat.name[, 2],
                                         comp.treat.name[, 1],
                                         sep = "_minus_")

  summ <- list(n            = n,
               raw.y.mean   = raw.y.mean,
               raw.y.median = raw.y.median,
               raw.y.sd     = raw.y.sd,
               post.coef.mean   = post.coef.mean,
               post.coef.median = post.coef.median,
               post.y.mean   = post.y.mean,
               post.y.median = post.y.median,
               post.coef.ci  = post.coef.ci,
               post.y.ci     = post.y.ci,
               comp.treat.post.coef  = comp.treat.post.coef,
               p.comp.coef.greater.0 = p.comp.coef.greater.0)

  comp.treat.name <- t(utils::combn(Treat.name, 2))
  names = paste(comp.treat.name[, 2], comp.treat.name[, 1], sep = "_minus_")

  res <- c(indiv)
  resnames <- c("ID")
  if (!is.null(groupvars)) {
    for (j in 1:length(groupvars)) {
      res <- c(res, unique(data.ind[, groupvars[j]]))
    }
    resnames <- c(resnames, groupvars)
  }

  if (!is.null(stratvars)) {
    for (k in 1:length(stratvars)) {
      res <- c(res, as.character(unique(data.ind[[stratvars[k]]])))
    }
    resnames <- c(resnames, stratvars)
  }

  for (trt in Treat.name) {
    if (!is.na(summ$post.coef.mean[trt])) {
      res <- c(res, summ$post.coef.mean[names(summ$post.coef.mean) == trt],
               summ$post.coef.median[names(summ$post.coef.median) == trt])
    }
    else {
      res <- c(res, rep(NA, 2))
    }
    resnames <- c(resnames, paste(trt, "_mean", sep = ""),
                  paste(trt, "_median", sep = ""))
  }

  for (i in names) {
    if (!is.na(summ$p.comp.coef.greater.0[i])) {
      if ((sum(!is.na(data.ind[[Y]][data.ind[, treatment] ==
                                    sub('.*minus_', '', i)])) != 0) &
          (sum(!is.na(data.ind[[Y]][data.ind[, treatment] ==
                                    sub('_minus_.*', '', i)])) != 0)) {
        res <- c(res, 1)
      } else {
        res <- c(res, 0)
      }
      res <-
        c(res,
          summ$comp.treat.post.coef[rownames(summ$comp.treat.post.coef) ==
                                      i, ],
          summ$p.comp.coef.greater.0[names(summ$p.comp.coef.greater.0) ==
                                       i])
      res <- c(res,
               mean((samples[, paste("beta", sub('_minus_.*', '', i),
                                     sep = "_")] -
                       samples[, paste("beta", sub('.*minus_', '', i),
                                       sep = "_")]) < -1*clinicaldiff))
      res <- c(res,
               mean((samples[, paste("beta", sub('_minus_.*', '', i),
                                     sep = "_")] -
                       samples[, paste("beta", sub('.*minus_', '', i),
                                       sep = "_")]) > clinicaldiff))
    } else {
      res <- c(res,
               rep(NA, 7))
    }
    resnames <- c(resnames,
                  "show",
                  "comp2.5,",
                  "comp50",
                  "comp97.5",
                  "compprob",
                  "secondbetter",
                  "firstbetter")
  }
  res <- setNames(res, resnames)
  summ.ind <- rbind(res)

  ncol = ncol(summ.ind)
  Treat.name =  as.vector(sort(unique(
    data[[treatment]])))
  ninfo = length(groupvars) + length(stratvars) + 1
  ncomp = choose(length(as.vector(unique(data[[treatment]]))), 2)
  extraparam = 2*length(as.vector(unique(data[[treatment]])))

  comp.treat.name <- t(utils::combn(Treat.name, 2))
  names = paste(comp.treat.name[, 2], comp.treat.name[, 1], sep = "_minus_")

  summ.columns <- list()
  for (i in 1:ncomp) {
    summ.columns[[i]] <- c(1:(ninfo + extraparam),
                           (ninfo + extraparam + 1 + 7*(i-1)):
                             (ninfo + extraparam + 7 + 7*(i-1)))
  }

  res.ind <- summ.ind[, summ.columns[[1]]]

  if (length(summ.columns) == 1) {
    res.ind = t(res.ind)
  }

  if (length(summ.columns) > 1) {
    for (i in 2:length(summ.columns)) {
      newrow <- summ.ind[, summ.columns[[i]]]
      res.ind <- rbind(res.ind, newrow)
    }
  }

  res.ind <- data.frame(res.ind)

  meanmednames <- c()
  for (trt in Treat.name) {
    meanmednames <- c(meanmednames, paste(trt, "_mean"), paste(trt, "_median"))
  }

  colnames(res.ind)  <- c("ID",
                          groupvars,
                          stratvars,
                          meanmednames,
                          "show",
                          "comp2.5",
                          "comp50",
                          "comp97.5",
                          "compprob",
                          "secondbetter",
                          "firstbetter")

  res.ind <- res.ind %>%
    mutate(ID = as.character(.data$ID))

  if (!is.null(groupvars)) {
    for (i in 1:length(groupvars)) {
      res.ind[[groupvars[i]]] <- as.factor(unlist(res.ind[[groupvars[i]]]))
    }
  }

  if (!is.null(stratvars)) {
    for (i in 1:length(stratvars)) {
      res.ind[[stratvars[i]]] <- as.factor(unlist(res.ind[[stratvars[i]]]))
    }
  }

  res.ind <- res.ind %>%
    mutate(comp2.5 = as.numeric(.data$comp2.5),
           comp50 = as.numeric(.data$comp50),
           comp97.5 = as.numeric(.data$comp97.5),
           compprob  = as.numeric(.data$compprob),
           secondbetter = as.numeric(.data$secondbetter),
           firstbetter = as.numeric(.data$firstbetter))

  res.ind <- res.ind %>%
    mutate(comp2.5 = ifelse(.data$show == 0, NA, comp2.5),
           comp50 = ifelse(.data$show == 0, NA, comp50),
           comp97.5 = ifelse(.data$show == 0, NA, comp97.5),
           compprob  = ifelse(.data$show == 0, NA, compprob),
           secondbetter = ifelse(.data$show == 0, NA, secondbetter),
           firstbetter = ifelse(.data$show == 0, NA, firstbetter))

  res.ind <- res.ind %>%
    mutate(neitherbetter = 1 - secondbetter - firstbetter)

  res.ind <- res.ind %>%
    mutate(comp = rep(gsub("_minus_", "-", names), each = nrow(summ.ind)))

  res.ind <- res.ind %>%
    mutate(comp = factor(.data$comp))

  noshow.ids <-
    names(tapply(res.ind$secondbetter,
                 paste0(res.ind$ID),
                 function(x)
                   sum(is.na(x))))[tapply(res.ind$secondbetter,
                                          paste0(res.ind$ID),
                                          function(x)
                                            sum(is.na(x))) == ncomp]
  res.ind <- res.ind %>%
    dplyr::filter(!(ID %in% noshow.ids)) %>%
    select(-.data$show)

  rownames(res.ind) <- NULL

  return(res.ind)
}
