#' @title Pathway-guided random forest based on synthetic features
#'
#' @description Trains a separate random forest (RF) for each pathway and
#' determines predictions for each of them. These predictions (called synthetic
#' features) are subsequently analyzed with the variable selection method
#' \code{\link[Boruta]{Boruta}} (Kursa et al. 2010) to select
#' important synthetic features and, hence, important pathways. For
#' more details about random forests parameters see \code{\link{wrapper.rf}} and
#' for more details about the method see Pan et al. (2014).
#'
#' @details In the original implementation (Pan et al. 2014) the RF based on the
#' synthetic features is only evaluated regarding its prediction performance.
#' We additionally implemented the Boruta method to select important synthetic
#' features, i.e. pathways.
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
#' @param ... Additional parameters used in the ranger or the associated
#' wrapper function.
#'
#' @import Boruta
#'
#' @return List with the following components:
#' \itemize{
#' \item \code{results.pw} data.frame with information about each pathway
#' \itemize{
#' \item \code{id} pathway identifier
#' \item \code{no.var} number of variables in pathway
#' \item \code{run.x} importances of attributes gathered in each importance
#' source run
#' \item \code{decision} Boruta decision (Confirmed, Rejected or Tentative)
#' \item \code{selected} 0 or 1 denoting if pw is selected
#' }
#' \item \code{pw.sel} vector with ids of selected pathways
#' \item \code{info.pw} actual used variables to pathways assignment
#' }
#'
#' @references
#' \itemize{
#' \item Pan, Q, et al. (2014) A system‐level pathway‐phenotype association
#'   analysis using synthetic feature random forest. Genetic Epi 38: 209-219.
#'   \url{https://dx.doi.org/10.1002/gepi.21794}.
#' \item Kursa, MB, & Rudnicki, WR. (2010) Feature selection with the
#'   Boruta package. J Stat Softw 36: 1-13.
#'   \url{http://dx.doi.org/10.18637/jss.v036.i11}.
#'   }
#'
#' #' @examples
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
#' res.sf = pw.rf.synthetic.features(x = sim.data.l$data.train[, -1],
#'                                   y = sim.data.l$data.train[, 1],
#'                                   info.pw = sim.data.l$info.pw,
#'                                   type = "classification",
#'                                   not.assoc.pw = FALSE)

#' @export
pw.rf.synthetic.features <-
  function(x = NULL,
           y = NULL,
           info.pw = NULL,
           type = "regression",
           importance = "impurity_corrected",
           min.no.var = 5,
           not.assoc.pw = TRUE,
           ...) {
    results.pw <- NULL
    pw.sel <- NULL

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

    res.pw.sf <-
      pw.specific.rfs(
        x = x,
        y = y,
        info.pw = info.pw,
        type = "regression",
        return.pred = TRUE,
        return.pred.error = FALSE,
        ...
      )
    x.sf <- res.pw.sf$pred.pw

    # sgu 25.04.2019: removed smd code
    # train RF based on synthetic features and select important variables using
    # Boruta
    res.boruta <-
      var.sel.boruta(
        x = x.sf,
        y = y,
        type = type,
        importance = importance,
        ...
      )

    res.imp <- res.boruta$info[colnames(info.pw),]

    results.pw <-
      data.frame(
        id = colnames(info.pw),
        no.var = apply(info.pw, 2, sum),
        res.imp,
        stringsAsFactors = FALSE
      )

    pw.sel <- rownames(results.pw)[which(results.pw$selected == 1)]

    # sgu 25.04.2019: removed smd code
    return(
      list(
        results.pw = results.pw,
        pw.sel = pw.sel,
        info.pw = info.pw
        )
      )
  }
