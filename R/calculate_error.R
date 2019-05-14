#' @title Prediction error calculation
#' @description Calculates prediction errors by comparing predicted values with
#' the true outcome values. For regression and probability mode, it will give
#' mean squared error (mse) and pseudo R-squared (rsq). For classification mode,
#' overall accuracy (acc), overall error (err), Matthews correlation coefficient
#' (mcc), sensitivity (sens) and specificity (spec) are returned.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}.
#' @param true Vector with true value for each sample.
#' @param test.set Matrix or data.frame with predictor variables of test set
#' with variables in columns and samples in rows (Note: missing values are not
#' allowed).
#'
#' @return numeric vector with two elements for regression and probability
#' estimation (rmse, rsq) and five elements for classification (acc, err, mcc,
#' sens, spec)
#'
#' @export
#'
calculate.error <- function(rf, true, test.set = NULL) {
  if (!is.null(test.set)) {
    pred <- predict(rf, data = test.set)$predictions
  } else {
    pred <- rf$predictions
  }
  if (rf$treetype == "Probability estimation") {
    pred <- pred[, 2]
  }

  if (rf$treetype == "Classification") {
    conf.matrix <- table(pred = pred, true = true)
    tp <- conf.matrix[2, 2]
    tn <- conf.matrix[1, 1]
    fn <- conf.matrix[2, 1]
    fp <- conf.matrix[1, 2]

    # accuracy
    acc <- (tp + tn) / sum(conf.matrix)

    # Matthews correlation coefficient
    mcc <-
      (tp * tn - fp * fn) / sqrt((tp + fn) * (tn + fp) * (tp + fp) * (tn +  fn))

    # sensitivity
    sens <- tp / (tp + fn)

    # specificity
    spec <- tn / (fp + tn)

    error <- c(
      err = 1 - acc,
      acc = acc,
      mcc = mcc,
      sens = sens,
      spec = spec
    )
  } else {
    mse <- sum((pred - true) ^ 2, na.rm = TRUE) / sum(!is.na(pred))

    # pseudo R-squared uses sum of squared differences divided by n
    # instead of variance!
    v <- sum((true - mean(true)) ^ 2) / length(true)
    rsq <- 1 - mse / v
    error <- c(mse = mse, rsq = rsq)
  }

  return(error)
}
