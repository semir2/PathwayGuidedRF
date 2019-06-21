#' @title Pathway-guided random forest based on pathway importance score
#'
#' @description Trains a random forest (RF) using all predictor variables and
#' estimates variable specific importance. A pathway importance score
#' (enrichment score) is calculated using the mean of the importance scores of
#' the pathway members. A standardized score is then computed for each pathway
#' which is used to obtain P values under the assumption that the standardized
#' scores are standardnormally distributed. For more details about random
#' forests parameters see \code{\link{wrapper.rf}} and for more details about
#' the method see Chen et al. (2012).
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
#' @param ... Additional parameters used in the \code{\link{wrapper.rf}}
#' function.
#'
#' @return List with the following components:
#' \itemize{
#' \item \code{results.pw} data.frame with information about each pathway
#' \itemize{
#' \item \code{id} pathway identifier
#' \item \code{no.var} number of variables in pathway
#' \item \code{x.score} the pathway enrichment score of the pathway
#' \item \code{z.score} the standardized pathway enrichment score
#' \item \code{pval} P value based on standard normal distribution
#' \item \code{pval.adj} Bonferroni adjusted P value
#' \item \code{selected} 0 or 1 denoting if pathway is selected
#' }
#' \item \code{pw.sel} vector with ids of selected pathways
#' \item \code{info.pw} actual used variables to pathways assignment
#' \item \code{importance.genes} gene specific importances
#' }
#'
#' @references
#' Chen, X. & Ishwaran, H. (2012) Pathway hunting by random survival
#' forests. Bioinformatics 29:99–105.
#' \url{https://doi.org/10.1093/bioinformatics/bts643}.
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
#' res.hunt = pw.rf.hunt(x = sim.data.l$data.train[, -1],
#'                       y = sim.data.l$data.train[, 1],
#'                       info.pw = sim.data.l$info.pw,
#'                       type = "classification",
#'                       not.assoc.pw = FALSE)
#'
#' @export
pw.rf.hunt <-
  function(x = NULL,
           y = NULL,
           info.pw = NULL,
           type = "regression",
           importance = "impurity_corrected",
           min.no.var = 5,
           not.assoc.pw = TRUE,
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

    # train RF with all variables
    rf <-
      wrapper.rf(
        x = x,
        y = y,
        type = type,
        importance = importance,
        ...
      )

    # calculate gene specific importance
    imp <- rf$variable.importance

    res <-
      calculate.pw.score.pval(
        imp = imp,
        info.pw = info.pw,
        low = FALSE
        )

    # final result
    pval.adj <- p.adjust(res$pval, method = "bonferroni")
    results.pw <-
      data.frame(
        id = colnames(info.pw),
        no.var = apply(info.pw, 2, sum),
        res,
        pval.adj = pval.adj,
        selected = as.numeric(pval.adj < 0.05),
        stringsAsFactors = FALSE
      )

    pw.sel <- rownames(results.pw)[which(results.pw$selected == 1)]

    return(
      list(
        results.pw = results.pw,
        pw.sel = pw.sel,
        info.pw = info.pw,
        importance.genes = imp
        )
      )
  }
