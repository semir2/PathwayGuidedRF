#' Prepare input function
#'
#' Prepare input and check for correctness.
#'
#' @param x Matrix or data.frame of predictor variables with variables in
#' columns and samples in rows (Note: missing values are not allowed).
#' @param y Vector with values of outcome variable (Note: will be converted
#' to factor if classification mode is used).
#' @param info.pw Matrix or data.frame assigning variables to pathways with
#' variables in rows and pathways in columns (0: variable does not belong to
#' pathway, 1: variable belongs to pathway).
#' @param min.no.var Minimal number of variables per pathway (Note: pathway
#' with less variables are removed or moved into 'not_associated_genes_pw').
#' Default is 5.
#' @param not.assoc.pw Combines all variables not associated to any pathway
#' or pathway with less then <min.no.var> into pathway
#' 'not_associated_genes_pw'. Default is TRUE.
#'
#' @return List with the following components:
#' \itemize{#' \item \code{info.pw} variable pathway table containing all
#' samples from <x> and a not associated vector if <not.assoc.pw>
#' \item \code{x} input data with possible corrected names variable names
#' }
#'
#' @noRd
prepare.input <-
  function(x = NULL,
           y = NULL,
           info.pw = NULL,
           min.no.var = 5,
           not.assoc.pw = TRUE) {
    # abbort if type <info.pw> is not data.frame or matrix
    if (!(is.data.frame(info.pw) || is.matrix(info.pw))) {
      stop(
        "info.pw must be of type matrix or better data.frame"
        )
    }
    # warning if type <info.pw> is matrix
    if (is.matrix(info.pw)) {
      warning("converted info.pw into data.frame")
      info.pw <- as.data.frame(info.pw)
    }

    # check for syntactically valid names
    nm <- colnames(x)
    cor_nm <- make.names(nm, unique = T, allow_ = F)
    idx_cor_nm <- which(!cor_nm %in% nm, useNames = F)
    if (length(idx_cor_nm != 0)) {
      df_nm <-
        data.frame(
          index = idx_cor_nm,
          original = nm[idx_cor_nm],
          corrected = cor_nm[idx_cor_nm]
          )
      warning(paste0("Replaced",
        length(idx_cor_nm),
        "variable names with syntactically valid names:",
        print.and.capture(df_nm)))
      colnames(x) <- cor_nm
      row.names(info.pw) <- cor_nm
    }

    # NA check
    if (anyNA(y) || anyNA(x) || anyNA(info.pw)) {
      stop("NAs in dataset, outcome or pathways. Missing values are not  allowed")
    }
    # entries <y>, <x> check
    if (nrow(x) != length(y)) {
      stop(
        "dimensions of dataset and outcome differ. Outcome vector and dataset",
        " must be the same size")
    }

    # reduce <info.pw> [variable pathway table] to variables in <x> [sample
    # variable  table]
    if (!all(rownames(info.pw) %in% colnames(x))) {
      info.pw <-
        info.pw[rownames(info.pw) %in% colnames(x), , drop = FALSE]
      warning("reduced variables in pathways to dataset variables")
    }

    # abbort if <info.pw> is empty
    if (nrow(info.pw) == 0 || ncol(info.pw) == 0) {
      stop(
        "pathways do not contain any variables from dataset"
        )
    }

    # remove small pathways with less than <min.no.var> variables
    # also removes pathways with no or zero variables assigned
    if (any(colSums(info.pw != 0) < min.no.var)) {
      info.pw <-
        info.pw[, min.no.var <= colSums(info.pw != 0), drop = FALSE]
      warning("pathways are removed because of small numbers of variables")
    }

    # add variables from <x> to <info.pw> that are not in <info.pw>
    if (not.assoc.pw) {
      rnames <-
        colnames(x)[!(colnames(x) %in% rownames(info.pw)), drop = FALSE]
      if (0 < length(rnames)) {
        not_associated_genes <-
          as.data.frame(
            matrix(
              0,
              ncol = ncol(info.pw),
              nrow = length(rnames)
              )
            )
        rownames(not_associated_genes) <- rnames
        colnames(not_associated_genes) <- colnames(info.pw)

        # circumvent zero-row and zero-column drop
        info.pw <-
          as.data.frame(
            rbind(
              as.matrix(info.pw),
              as.matrix(not_associated_genes)
              )
            )
      }

      # define additional pathway to contain all variables not assigned to any
      # pathway or assigned to small pathways
      not_associated_genes_pw <- rep(0, nrow(info.pw))
      not_associated_genes_pw[which(rowSums(info.pw) == 0)] <- 1

      if (0 < sum(not_associated_genes_pw)) {
        info.pw <- cbind(info.pw, remaining_genes = not_associated_genes_pw)
      } else {
        warning(paste0("no additional pathway was added, as all variables are",
            "already associated to a pathway"))
      }
    }
    return(list(info.pw = info.pw, x = x))
  }

#' convert x into a string
#' @param obj object to print
#' @param row.names print row names of the obj.  Default is false.
#' @noRd
print.and.capture <- function(obj, row.names = FALSE) {
  paste(utils::capture.output(print(obj, row.names = row.names)),
        collapse = "\n")
}

#' calculate pathway scores and p-value for pw.rf.hunt()
#' @param imp variable importance for all variables
#' @param info.pw matrix or data.frame assigning variables to pathways with
#' variables in rows and pathways in columns (0: variable does not belong to
#' pathway, 1: variable belongs to pathway)
#' @param low logical; if TRUE probabilities are P[X ≤ x] otherwise, P[X > x].
#' @noRd
calculate.pw.score.pval <- function(imp, info.pw, low) {
  ## importance for each pathway
  imp.pw <- lapply(seq_len(ncol(info.pw)), function(i) {
    imp[rownames(info.pw)[which(info.pw[, i] != 0)]]
  })
  len <- vapply(imp.pw, length, FUN.VALUE = numeric(1))

  ## pathway enrichment score
  x <- vapply(imp.pw, mean, FUN.VALUE = numeric(1))

  ## standardized pathway enrichment score
  mu <- mean(imp)
  p <- length(imp)
  sigma.1 <-
    vapply(len, function(m)
      (p - m) / (m * (p - 1)), FUN.VALUE = numeric(1))
  sigma.2 <- mean(imp ^ 2) - mu ^ 2
  var <- sigma.1 * sigma.2
  z <- (x - mu) / sqrt(var)

  ## p-value based on standard normal distribution
  pval <- pnorm(z, lower.tail = low)

  return(
    data.frame(
      x.score = x,
      z.score = z,
      pval = pval
      )
    )
}

#' calculate empirical p-value
#' @import stats
#' @param test permuted data
#' @param x original data
#' @param y number of permutations
#' @param alternative alternative hypothesis
#' @noRd
calculate.var.weighting <- function(test, x, y, alternative) {
  switch(test,
         "wilcox.test" = {
           wt <- wilcox.test(
             x = x,
             y = y,
             alternative = alternative
             )
         },
         "t.test" = {
           wt <- t.test(
             x = x,
             y = y,
             alternative = alternative
             )
         },
         "ks.test" = {
           # alternative = 'greater' includes distributions
           #for which x is stochastically smaller than y
           if (alternative == "greater") {
             alternative <- c("less")
           } else {
             alternative <- c("greater")
           }
           wt <- ks.test(
             x = x,
             y = y,
             alternative = alternative
             )
         },
         {
           stop(paste0("test available: wilcox.test, t.test or ks.test. received(",
                test,
                ")"))
         })
  return(wt)
}

#' calculate empirical p-value
#' @param perm permuted data
#' @param org original data
#' @param no.p number of permutations
#' @noRd
calculate.empirical.pvalue <- function(perm, org, no.p) {
  #check if perm has dimensions
  if (is.null(dim(perm)))
    stop(paste0(
      "calculate empirical p-value must have perm data with dimension bigger 1. got '",
      dim(perm),
      "'"))

  pval.emp <-
    vapply(
      seq_len(nrow(perm)),
      function(x) {
        temp <- perm[x, seq_len(no.p)]
        no.NA <- sum(is.na(temp))
        (sum(temp < org[x], na.rm = TRUE) + 1) / (no.p - no.NA + 1)
        },
      FUN.VALUE = numeric(1)
      )
  names(pval.emp) <- rownames(perm)
  return(pval.emp)
}

#' Variable selection using Boruta function.
#' Variable selection using the Boruta function in the R package
#' \code{\link[Boruta]{Boruta}}.
#' This function selects only variables that are confirmed based on Boruta
#' implementation.
#' For more details see \code{\link[Boruta]{Boruta}}.
#' Note that this function uses an internally defined getImp function that is
#' based on the wrapper.rf function so that
#' user specified importance scores and mtry.prop values can be employed.
#' @param min.node.size parameter for ranger - Minimal node size. Default 1 for
#' classification, 5 for regression, 3 for survival, and 10 for probability.
#' @param num.threads parameter for ranger - Number of threads. Default is
#' number of CPUs available. Default is 500.
#' @inheritParams wrapper.rf
#' @inheritParams Boruta::Boruta
#' @return List with the following components:
#' \itemize{#' \item \code{info} data.frame with information of each variable
#' \itemize{#' \item run.x: original variable importance (VIM) in run x
#' \item decision: Boruta decision (Confirmed, Rejected or Tentative)
#' \item selected: variable has been selected
#' }
#' \item \code{var} vector of selected variables
#' }
#' @references
#' \itemize{
#'   \item Kursa, Miron B., and Witold R. Rudnicki. "Feature selection with the
#'   Boruta package." J Stat Softw 36.11 (2010): 1-13.
#'   }
#' @noRd
var.sel.boruta <-
  function(x = NULL,
           y = NULL,
           pValue = 0.01,
           maxRuns = 100,
           num.trees = 500,
           mtry.prop = NULL,
           min.node.size = 1,
           num.threads = 1,
           type = "regression",
           importance = "impurity_corrected",
           ...) {
    # modified version of getImpRfRaw function to enable user defined mtry
    # values
    get.imp.r.f.raw.mtry <- function(x, y, ...) {
      rf <- wrapper.rf(x = x, y = y, ...)
      return(rf$variable.importance)
    }

    # variable selection using Boruta function
    res.boruta <-
      Boruta::Boruta(
        x = x,
        y = y,
        pValue = pValue,
        maxRuns = maxRuns,
        num.trees = num.trees,
        min.node.size = min.node.size,
        mtry.prop = mtry.prop,
        importance = importance,
        num.threads = num.threads,
        getImp = get.imp.r.f.raw.mtry,
        ...
      )
    # select variables
    dec <- res.boruta$finalDecision
    ind.sel <- rep(0, ncol(x))
    ind.sel[dec == "Confirmed"] <- 1
    info.sel <- data.frame(decision = dec, selected = ind.sel)

    # info about variables
    info.var <- t(res.boruta$ImpHistory)
    colnames(info.var) <-
      paste(
        "run",
        seq_len(ncol(info.var)),
        sep = "."
        )
    info.var <-
      info.var[-grep("shadow", rownames(info.var)), , drop = FALSE]
    if (all.equal(rownames(info.var), rownames(info.sel))) {
      info.var <- cbind(info.var, info.sel)
    } else {
      info.var <-
        merge(info.var, info.sel, by.x = "row.names", by.y = "row.names")
    }

    return(
      list(
        info = info.var,
        var = sort(rownames(info.var)[info.var$selected == 1])
        )
      )
  }

#' train pathway specific RFs
#' Random forest is build by the ranger function in the R package
#' \code{\link[ranger]{ranger}}.
#' @param x sample variable table - data of class data.frame. Default is 'NULL'.
#' @param y vector with values of phenotype variable (Note: will be converted to
#' factor if classification mode is used). Default is 'NULL'.
#' @param info.pw variable pathway table. Default is 'NULL'.
#' @param return.pred return predictions of synthetic features based on ranger.
#' Default is false. This will set <type> to 'regression' and <importance> to
#' 'none'.
#' @param return.pred.error return predictions errors of pathways. Default is
#' true.
#' @param type mode of prediction ('regression' or 'classification'). Default
#' is 'regression'. In case of <return.pred> it is set to 'regression'
#' @param importance variable importance mode ('none', 'impurity',
#' 'impurity_corrected' or 'permutation'). Default is 'impurity_corrected'.
#' In case of <return.pred> it is set to 'none'
#' @param ... additional parameters used in the ranger function.
#' @noRd
pw.specific.rfs <-
  function(x = NULL,
           y = NULL,
           info.pw = NULL,
           type = "regression",
           return.pred = FALSE,
           return.pred.error = TRUE,
           importance = "impurity_corrected",
           ...) {
    if (return.pred) {
      type <- "regression"
      importance <- "none"
    }

    ## RF for each pathway
    pw.rf.l <- lapply(seq_len(ncol(info.pw)), function(i) {
      x.pw <- x[, rownames(info.pw)[which(info.pw[, i] != 0)]]
      wrapper.rf(
        x = x.pw,
        y = y,
        importance = importance,
        type = type,
        ...
      )
    })
    names(pw.rf.l) <- colnames(info.pw)

    if (return.pred) {
      pred.pw <- vapply(
        pw.rf.l,
        function(r) {
          r$predictions
          },
        FUN.VALUE = numeric(length(y))
        )
    } else {
      pred.pw <- NULL
    }

    if (return.pred.error) {
      pred.error.pw <- vapply(
        pw.rf.l,
        function(r) {
          r$prediction.error
          },
        FUN.VALUE = numeric(1)
        )
    } else {
      pred.error.pw <- NULL
    }

    return(
      list(
        pred.error.pw = pred.error.pw,
        pred.pw = pred.pw)
      )
  }
