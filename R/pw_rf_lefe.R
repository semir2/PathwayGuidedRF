#' @title Pathway-guided random forest based on Learner of Functional Enrichment
#'
#' @description Trains a random forest (RF) for each pathway using the pathway
#' variables and a proportionally sized set of randomly selected additional
#' variables that are not contained in the pathway. Subsequently, a statistical
#' test is applied to compare the importance scores of the pathway variables
#' with the additional variables. These steps are repeated multiple times and
#' significant pathways are selected based on median P values. For more details
#' about random forests parameters see \code{\link{wrapper.rf}} and for more
#' details about the method see Eichler et al. (2007).
#'
#' @details The default values of the method specific parameters
#' \code{sample.factor} and \code{samples.runs} are based on Eichler et. al.
#' (2007). The default statistical test that is used in this implementation is
#' the Wilcoxon test unlike in the original publication where a permutation test
#' based on the statistic of the t-test is used. (see Eichler et. al. (2007)).
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
#' @param sample.factor Factor to determine the size of the set of
#' additional variables. Default is 6 meaning that six times more additional
#' variables than pathway variables are used.
#' @param sample.runs Number of sample runs that are conducted. The median
#' value of the P values of the sample runs is used to select pathways. Default
#' is 75.
#' @param test The statistical test that is applied. ('wilcox.test', 't.test'
#' or 'ks.test') Default is 'wilcox.test'.
#' @param ... Additional parameters used in the ranger or the associated wrapper
#' function.
#'
#' @return List with the following components:
#' \itemize{
#' \item \code{results.pw} data.frame with information about each pathway
#' \itemize{
#' \item \code{id} pathway identifier
#' \item \code{no.var} number of variables in pathway
#' \item \code{pval} P value
#' \item \code{pval.adj} Bonferroni adjusted P value
#' \item \code{selected} 0 or 1 denoting if pathway is selected
#' }
#' \item \code{pw.sel} vector with ids of selected pathways
#' \item \code{info.pw} actually used variables to pathways assignment
#' }
#'
#' @references
#' Eichler, GS, et al. (2007) The LeFE algorithm: embracing the
#' complexity of gene expression in the interpretation of microarray data.
#' Genome Biol 8:R187. \url{https://doi.org/10.1186/gb-2007-8-9-r187}
#'
#' @examples
#' # define pathway parameters
#' info.par = data.frame(pw = paste0("pw", 1:10),
#'                       no.genes = rep(20, 10),
#'                       prop.de = c(rep(0, 8), 0.5, 1),
#'                       cor = 0.6,
#'                       stringsAsFactors = FALSE)
#'
#' # simulate 30 samples (including 15 cases)
#' sim.data.l = sim.data.study.1(info.par = info.par,
#'                               no.samples = 30,
#'                               no.cases = 15,
#'                               seed = 42,
#'                               gamma.range = c(0.5, 1.5))
#'
#' # Note: sample.runs of 10 only used for illustration
#' res.lefe = pw.rf.lefe(x = sim.data.l$data.train[, -1],
#'                      y = sim.data.l$data.train[, 1],
#'                      info.pw = sim.data.l$info.pw,
#'                      type = "classification",
#'                      not.assoc.pw = FALSE,
#'                      sample.runs = 10)
#'
#' @export
pw.rf.lefe <-
  function(x = NULL,
           y = NULL,
           info.pw = NULL,
           type = "regression",
           importance = "impurity_corrected",
           min.no.var = 5,
           not.assoc.pw = TRUE,
           sample.factor = 6,
           sample.runs = 75,
           test = "wilcox.test",
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

    pvalues <- array(as.double(NA), dim = dim(info.pw)[2])
    names(pvalues) <- colnames(info.pw)
    impm.pw <- list()
    for (pw in colnames(info.pw)) {
      E <- intersect(rownames(info.pw[which(info.pw[, pw] != 0), ,
                                      drop = FALSE]),
                     colnames(x))
      if (length(E) > 0) {
        work.E <- x[, E]
        nE <- ncol(work.E)
        notE <- setdiff(colnames(x), E)
        work.notE <- x[, notE]
        nnE <- ncol(work.notE)
        if (nnE >= nE * sample.factor) {
          impE <- array(0, dim = c(sample.runs, nE))
          colnames(impE) <- colnames(work.E)
          pv_vec <- vector()
          for (i in seq_len(sample.runs)) {
            work.sample <-
              cbind(
                work.E,
                work.notE[, sample(seq_len(nnE), nE * sample.factor)]
                )
            work.sample <- as.data.frame(work.sample)
            work.sample$response <- y

            rf <-
              wrapper.rf(
                x = work.sample[,!(names(work.sample) %in% c("response"))],
                y = work.sample[, "response"],
                type = type,
                importance = importance,
                ...
              )
            imp <- ranger::importance(rf)
            alternative <- c("greater")
            impE[i,] <- imp[seq_len(nE)]
            wt <-
              calculate.var.weighting(
                test = test,
                x = imp[seq_len(nE)],
                y = imp[(nE + 1):(nE * (1 + sample.factor))],
                alternative = alternative
              )
            pv_vec[i] <- wt$p.value
          }
          impm <- apply(impE, 2, median)
          names(impm) <- colnames(work.E)
          impm.pw[[pw]] <- impm
          pvalues[pw] <- median(pv_vec)
        } else {
            warning(paste0("too few genes for sampling without category: ",
                           pw,
                           " (",
                           nE,
                           "*",
                           sample.factor,
                           " =< ",
                           nnE,
                           ")"))
        }
      } else {
        warning(paste0("empty gene set for category ", pw))
      }
    }
    pval.adj <- p.adjust(pvalues, method = "bonferroni")

    results.pw <-
      data.frame(
        id = colnames(info.pw),
        no.var = apply(info.pw, 2, sum),
        pval = pvalues,
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
