
#' @title Simulate data sets for simulation study 1
#'
#' @description Simulates pathways with specified number of genes, amount of
#' differential expression and pairwise correlation.
#'
#' @details This function uses a modifed version of R code provided by Poisson
#' et al. (2011) as supplementary material.
#'
#' @param info.par Information about pathways as data.frame with columns
#' no.genes (number of genes in pathway), prop.de (proportion of genes to be
#' simulated as differentially expressed) and cor (pairwise correlation between
#' genes).
#' @param no.samples Total number of individuals to simulate.
#' @param no.cases Number of cases to simulate.
#' @param seed Random seed. Default is NULL, which generates the seed from R.
#' @param dd1 Probability of up regulation.
#' @param gamma.range Range of absolute effect sizes.
#'
#' @return List with the following components:
#' \itemize{
#' \item \code{data.train} data.frame of training data with samples in rows and
#' variables in columns; the first column, called y, denotes the case-control
#' status (0: control, 1: case)
#' \item \code{data.test} data.frame of test data with samples in rows and
#' variables in columns; the first column, called y, denotes the case-control
#' status (0: control, 1: case)
#' \item \code{info.pw} data.frame assigning variables to pathways with
#' variables in rows and pathways in columns (0: variable does not belong to
#' pathway, 1: variable belongs to pathway).
#' }
#'
#' @references
#' Poisson LM, Taylor JM, Ghosh D. (2011) Integrative set enrichment testing
#' for multiple omics platforms. BMC Bioinformatics. 12:459.
#' \url{https://doi.org/10.1186/1471-2105-12-459}.
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
#' @export

sim.data.study.1 <- function(info.par,
                             no.samples,
                             no.cases,
                             seed,
                             dd1 = 0.5,
                             gamma.range) {

    no.genes.total = sum(info.par$no.genes)

    set.seed(seed)
    expr.1 = NULL
    expr.2 = NULL
    info.pw = NULL
    ind.gene = 0
    for (p in 1:nrow(info.par)) {
        par.pw = info.par[p,]
        data.p.l = data.gen.DJsim.mod(n.samp = no.samples,
                                      n.case = no.cases,
                                      n.path = 1,
                                      n.assoc.path = 1,
                                      n.gene.path = par.pw$no.genes,
                                      d1 = par.pw$prop.de,
                                      dd1 = dd1,
                                      gamma.range = gamma.range,
                                      corr.gg = par.pw$cor)
        expr.1 = cbind(expr.1, t(data.p.l[[1]]))
        expr.2 = cbind(expr.2, t(data.p.l[[2]]))

        info = rep(0, no.genes.total)
        info[(ind.gene + 1):(ind.gene + par.pw$no.genes)] = 1
        info.pw = cbind(info.pw, info)
        ind.gene = ind.gene + par.pw$no.genes
    }

    dimnames(info.pw) = list(paste0("gene", 1:no.genes.total),
                             paste0("pw", 1:ncol(info.pw)))
    dimnames(expr.1) = dimnames(expr.2) = list(paste0("ind", 1:no.samples),
                                               paste0("gene", 1:no.genes.total))

    ## phenotype data
    outcome = rep(0, no.samples)
    outcome[1:no.cases] = 1
    expr.1 = data.frame(y = outcome, expr.1)
    expr.2 = data.frame(y = outcome, expr.2)

    return(list(data.train = expr.1,
                data.test = expr.2,
                info.pw = data.frame(info.pw)))
}
