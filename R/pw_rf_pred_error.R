#' @title Pathway-guided random forest based on prediction error
#'
#' @description Trains a separate random forest (RF) for each pathway and
#' selects pathways with prediction error significantly lower than expected by
#' chance. An empirical P value is calculated by comparing the observed
#' prediction error of the RF trained on the original data with the distribution
#' of prediction errors of RFs trained on data sets with permuted outcome. For
#' more details about random forests parameters see \code{\link{wrapper.rf}} and
#' for more details about the method see Pang et al. (2006).
#'
#' @details The method specific parameters \code{no.perm} and \code{adaptive}
#' can be adapted. Note that we introduced the latter parameter to reduce run
#' time, which was not used in the original implementation (Pang et al. 2006).
#'
#' @inheritParams wrapper.rf
#' @param info.pw Matrix or data.frame assigning variables to pathways with
#' variables in rows and pathways in columns (0: variable does not belong to
#' pathway, 1: variable belongs to pathway).
#' @param min.no.var Minimal number of variables per pathway (Note: pathways
#' with less variables are removed or moved into 'not_associated_genes_pw').
#' Default is 5.
#' @param not.assoc.pw Combines all variables not associated to any pathway
#' or pathways with less then <min.no.var> into pathway
#' 'not_associated_genes_pw'. Default is TRUE.
#' @param no.perm Number of permutations used to determine the P value. Default
#' is 1000.
#' @param adaptive Logical; this parameter can be set to TRUE to reduce run
#' time. When set to TRUE the number of permutations are reduced for pathways
#' that have high P values (>= 0.2). Default is TRUE.
#' @param ... Additional parameters used in the ranger or the associated
#' wrapper function.
#'
#' @return List with the following components:
#' \itemize{
#' \item \code{results.pw} data.frame with information about each pathway
#' \itemize{
#' \item \code{id} pathway identifier
#' \item \code{no.var} number of variables in pathway
#' \item \code{pred.error} prediction error of pathway specific RF
#' \item \code{pred.error.perm.min}, \code{pred.error.perm.median},
#' \code{pred.error.perm.max} summary information about prediction errors after
#' permuting the outcome
#' \item \code{no.perm.used} number of permutations calculated for each pathway
#' \item \code{pval} empirical P value
#' \item \code{pval.adj} Bonferroni adjusted empirical P value
#' \item \code{selected} 0 or 1 denoting if pathway is selected
#' }
#' \item \code{pw.sel} vector with ids of selected pathways
#' \item \code{info.pw} actually used variables to pathways assignment
#' }
#'
#' @references
#' Pang, H et al. (2006) Pathway analysis using random forests
#' classification and regression. Bioinformatics 22: 2028-2036.
#' \url{https://doi.org/10.1093/bioinformatics/btl344}.
#'
#' @examples
#' # define pathway parameters
#' info.par = data.frame(pw = paste0("pw", 1:3),
#'                       no.genes = rep(20, 3),
#'                       prop.de = c(0, 0.5, 1),
#'                       cor = rep(0.6, 3),
#'                       stringsAsFactors = FALSE)
#'
#' # simulate 30 samples (including 15 cases)
#' sim.data.l = sim.data.study.1(info.par = info.par,
#'                               no.samples = 30,
#'                               no.cases = 15,
#'                               seed = 42,
#'                               gamma.range = c(0.5, 1.5))
#'
#' # Note: no.perm of 100 only used for illustration
#' res.pe = pw.rf.pred.error(x = sim.data.l$data.train[, -1],
#'                           y = sim.data.l$data.train[, 1],
#'                           info.pw = sim.data.l$info.pw,
#'                           type = "classification",
#'                           not.assoc.pw = FALSE,
#'                           no.perm = 100)
#'
#' @export
pw.rf.pred.error <-
  function(x = NULL,
           y = NULL,
           info.pw = NULL,
           type = "regression",
           min.no.var = 5,
           not.assoc.pw = TRUE,
           no.perm = 1000,
           adaptive = TRUE,
           ...) {
    pr_inp <-
      prepare.input(
        x = x,
        y = y,
        info.pw = info.pw,
        min.no.var = min.no.var,
        not.assoc.pw = not.assoc.pw
      )
    info.pw <- pr_inp$info.pw
    x <- pr_inp$x

    # pathway specific prediction errors for original outcome values
    res.pw.org <-
      pw.specific.rfs(
        x = x,
        y = y,
        info.pw = info.pw,
        type = type,
        return.pred = FALSE,
        return.pred.error = TRUE,
        ...
      )
    pred.error.org <- res.pw.org$pred.error.pw

    # pathway specific prediction errors after permuting the outcome
    pred.error.perm <- matrix(nrow = ncol(info.pw), ncol = no.perm)
    rownames(pred.error.perm) <- colnames(info.pw)

    # determine checkpoints for adaptive P values
    if (adaptive) {
      s <- seq_len(log10(no.perm))
      cp <- sort(c(10 ^ s[-1], 5 * 10 ^ s[-1]))
      cp <- cp[which(cp <= no.perm)]
      if (length(cp) == 0) {
        adaptive <- FALSE
        warning("No adaptive P value calculation since no.perm is too low")
      }
    }

    pw.use <- colnames(info.pw)
    for (p in seq_len(no.perm)) {
      # permute outcome
      y.per <- sample(y, length(y), replace = FALSE)

      # compute prediction error per pathway
      res.pw.perm <-
        pw.specific.rfs(
          x = x,
          y = y.per,
          info.pw = info.pw[, pw.use, drop = FALSE],
          type = type,
          return.pred = FALSE,
          return.pred.error = TRUE,
          ...
        )

      # save prediction error
      prd.err.pw <- res.pw.perm$pred.error.pw
      pred.error.perm[names(prd.err.pw), p] <- prd.err.pw

      if (adaptive) {
        if (p %in% cp) {
          # compute pval
          pval.temp <-
            calculate.empirical.pvalue(
              perm = pred.error.perm[pw.use, ,drop = FALSE],
              org = pred.error.org,
              no.p = p
              )

          # remove pathways with large P values from further calculations
          pw.use <- setdiff(pw.use, names(pval.temp)[which(pval.temp > 0.2)])
          if (length(pw.use) == 0) {
            warning("Early stopping of permutations because all P values > 0.2!"
                    )
          }
        }
      }
    }

    # summarize permutation based prediction errors
    info.pred.error.perm <-
      t(apply(
        pred.error.perm,
        1,
        quantile,
        probs = c(0, 0.5,
                  1),
        na.rm = TRUE,
        names = FALSE
      ))
    colnames(info.pred.error.perm) <-
      paste(
        "pred.error.perm",
        c(
          "min",
          "median",
          "max"
          ),
        sep = "."
        )
    no.perm.used <- apply(
      pred.error.perm,
      1,
      function(i) {
        sum(!is.na(i))
        }
      )

    # calculate empirical p-value
    pval <-
      calculate.empirical.pvalue(
        perm = pred.error.perm,
        org = pred.error.org,
        no.p = no.perm
        )
    pval.adj <- p.adjust(pval, method = "bonferroni")

    results.pw <-
      data.frame(
        id = colnames(info.pw),
        no.var = apply(info.pw, 2, sum),
        pred.error = pred.error.org,
        info.pred.error.perm,
        no.perm.used = no.perm.used,
        pval = pval,
        pval.adj = pval.adj,
        selected = as.numeric(pval.adj < 0.05),
        stringsAsFactors = FALSE
      )

    pw.sel <- rownames(results.pw)[which(results.pw$selected == 1)]

    return(
      list(
        results.pw = results.pw,
        pw.sel = pw.sel,
        info.pw = info.pw
      )
    )
  }
