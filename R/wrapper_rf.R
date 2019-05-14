#' @title Wrapper function to call ranger
#'
#' @description Provides an interface to the random forest algorithm implemented
#' in \code{\link[ranger]{ranger}}. For more details about random forests
#' parameters see \code{\link[ranger]{ranger}}.
#'
#' @param x Matrix or data.frame of predictor variables with variables in
#' columns and samples in rows (Note: missing values are not allowed).
#' @param y Vector with values of outcome variable (Note: will be converted
#' to factor if classification mode is used).
#' @param num.trees Number of trees. Default is 500.
#' @param mtry.prop Proportion of variables that should be used at each split.
#' Default value is the number of variables^(3/4) as recommended by
#' Ishwaran et al. (2011) for minimal depth importance. (Note: This is a different
#' default value than in the \code{\link[ranger]{ranger}} function).
#' @param type Type of random forest ('regression' or 'classification'). Default
#' is 'regression'.
#' @param importance Variable importance mode ('none', 'impurity',
#' 'impurity_corrected' or 'permutation'). Default is 'impurity_corrected'.
#' @param seed Random seed. Default is NULL, which generates the seed from R.
#' Set to 0 to ignore the R seed.
#' @param ... Additional parameters used in the ranger function.
#'
#' @import ranger
#'
#' @return Object of class ranger
#'
#' @references
#' \itemize{
#'   \item Ishwaran H, Kogalur UB, Chen X, Minn AJ (2011). Random survival
#'   forests for high-dimensional data. Stat Anal Data Min 4:115–32.
#'   \url{ https://doi.org/10.1002/sam.10103}.
#'   \item Wright, MN & Ziegler, A (2017). ranger: A fast implementation of
#'   random forests for high dimensional data in C++ and R. J Stat Softw
#'   77: 1-17. \url{https://doi.org/10.18637/jss.v077.i01}.
#'   }
#' @export
wrapper.rf <-
  function(x = NULL,
           y = NULL,
           num.trees = 500,
           mtry.prop = NULL,
           type = "regression",
           importance = "impurity_corrected",
           # sgu 25.04.2019: removed smd code
           seed = NULL,
           ...) {
    # check data
    if (is.null(x) || is.null(y)) {
      stop("Training data is NULL")
    }

    # convert outcome vector based on type
    switch(type, regression = {
      if (is.factor(y)) {
        warning(paste0(
          "type is regression but outcome is of type factor. Outcome will be",
          " converted to integer with x-1."))
        y <- as.numeric(y) - 1
      }
    }, classification = {
      # warning if type == 'classification' and y has more than 2 different
      # values
      if (nlevels(y) > 2) {
        warning(paste0(
          "type is classification but outcome has more than 2 different",
          " values. Will proceed with classification."))
      }
      y <- factor(y)
    }, {
      stop("type must be regression or classification")
    })


    # set mtry
    no.var <- ncol(x)

    if (is.null(mtry.prop)) {
      mtry <- floor(no.var ^ (3 / 4))
    } else {
      mtry <- floor(mtry.prop * no.var)
    }

    if (mtry <= 0) {
      mtry <- 1
    }

    # run ranger
    res <-
      ranger::ranger(
        data = data.frame(y, x),
        num.trees = num.trees,
        mtry = mtry,
        importance = importance,
        write.forest = TRUE,
        probability = FALSE,
        scale.permutation.importance = FALSE,
        keep.inbag = TRUE,
        dependent.variable.name = "y",
        seed = seed,
        ...
      )
    if (anyNA(res$predictions)) {
      warning(paste0(
        "missing data in predictions. Advice: increase number of trees.",
        " Default number of trees is 500."
        ))
    }

    return(res)
  }
