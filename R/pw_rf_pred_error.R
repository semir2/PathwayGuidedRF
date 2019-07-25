#' @title Pathway-guided random forest based on prediction error
#'
#' @description Trains a separate random forest (RF) for each pathway and
#' selects pathways with out-of-bag (OOB) prediction error significantly lower
#' than expected by chance using a permutation approach proposed by Hedinger et
#' al. (2019). OOB prediction errors are calculated for RFs trained on data sets
#' with permuted outcome and the mean and standard deviation of all those
#' prediction errors (across pathways and repetitions) are estimated. The
#' P values are then derived using a normal distribution. For
#' more details about random forests parameters see \code{\link{wrapper.rf}} and
#' for more details about the method see Pang et al. (2006).
#'
#' @details The method specific parameter \code{no.perm} can be adapted.
#' However, it is recommended to use a sufficiently large number so that
#' no.perm * number of pathways > 1000.
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
#' @param p.adjust.method Method for multiple testing adjusting (see
#' \code{\link{p.adjust}}). Default is "BH" which is different from the
#' adjustment method used prior to version 0.4.0 ("bonferroni").
#' @param no.perm Number of permutations used to determine the P value. Default
#' is 20.
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
#' \item \code{pval} permutation based P value
#' \item \code{pval.adj} adjusted permutation based P value
#' \item \code{selected} 0 or 1 denoting if pathway is selected
#' }
#' \item \code{pw.sel} vector with ids of selected pathways
#' \item \code{info.pw} actually used variables to pathways assignment
#' }
#'
#' @references
#' \itemize{
#' \item Hediger, S et al. (2019) On the use of random forest for two-sample testing.
#' arXiv. \url{https://arxiv.org/abs/1903.06287v2}
#' \item Pang, H et al. (2006) Pathway analysis using random forests
#' classification and regression. Bioinformatics 22: 2028-2036.
#' \url{https://doi.org/10.1093/bioinformatics/btl344}.
#' }
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
#' # Note: no.perm of 20 only used for illustration (should be larger for
#' # analysis of a small number of pathways)
#' res.pe = pw.rf.pred.error(x = sim.data.l$data.train[, -1],
#'                           y = sim.data.l$data.train[, 1],
#'                           info.pw = sim.data.l$info.pw,
#'                           type = "classification",
#'                           not.assoc.pw = FALSE,
#'                           no.perm = 20)
#'
#' @export
pw.rf.pred.error <- function(x = NULL,
                             y = NULL,
                             info.pw = NULL,
                             type = "regression",
                             min.no.var = 5,
                             not.assoc.pw = TRUE,
                             p.adjust.method = "BH",
                             no.perm = 20,
                             ...) {

    pr_inp <- prepare.input(x = x,
                            y = y,
                            info.pw = info.pw,
                            min.no.var = min.no.var,
                            not.assoc.pw = not.assoc.pw)
    info.pw <- pr_inp$info.pw
    x <- pr_inp$x

    ## observed OOB prediction error for each pathway
    res.pw.obs <- pw.specific.rfs(x = x,
                                  y = y,
                                  info.pw = info.pw,
                                  type = type,
                                  return.pred = FALSE,
                                  return.pred.error = TRUE,
                                  ...)
    pred.error.obs <- res.pw.obs$pred.error.pw

    ## prediction errors on permuted data
    pred.error.perm <- matrix(nrow = ncol(info.pw),
                              ncol = no.perm)
    rownames(pred.error.perm) <- colnames(info.pw)
    pw.use <- colnames(info.pw)
    for (p in seq_len(no.perm)) {
        y.per <- sample(y, length(y), replace = FALSE)
        res.pw.perm <- pw.specific.rfs(x = x,
                                       y = y.per,
                                       info.pw = info.pw[, pw.use,
                                                         drop = FALSE],
                                       type = type,
                                       return.pred = FALSE,
                                       return.pred.error = TRUE,
                                       ...)
        pred.error.perm[, p] <- res.pw.perm$pred.error.pw
    }
    info.pred.error.perm <- data.frame(t(apply(pred.error.perm, 1, quantile,
                                               probs = c(0, 0.5, 1),
                                               na.rm = TRUE,
                                               names = FALSE)))
    colnames(info.pred.error.perm) <- paste("pred.error.perm",
                                            c("min", "median", "max"),
                                            sep = ".")

    ## estimate P value based on mean and sd of permuted prediction errors
    pval <- pnorm(q = pred.error.obs,
                 mean = mean(as.numeric(pred.error.perm)),
                 sd = sd(as.numeric(pred.error.perm)),
                 lower.tail = TRUE)

    pval.adj <- p.adjust(pval, method = p.adjust.method)
    results.pw <- data.frame(id = colnames(info.pw),
                             no.var = apply(info.pw, 2, sum),
                             pred.error = pred.error.obs,
                             info.pred.error.perm,
                             pval = pval,
                             pval.adj = pval.adj,
                             selected = as.numeric(pval.adj < 0.05),
                             stringsAsFactors = FALSE)
    pw.sel <- rownames(results.pw)[which(results.pw$selected == 1)]
    return(list(results.pw = results.pw,
                pw.sel = pw.sel,
                info.pw = info.pw))
}
