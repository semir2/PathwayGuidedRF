
#' @title Simulate data sets for simulation study 2
#'
#' @description Simulates pathway based gene expression data using given
#' covariance matrix. A proportion of genes of a single pathway is randomly
#' selected and assigned a random selection of given effect sizes.
#'
#' @details Randomly selected effect sizes are assigned to a random subset of
#' the genes in the selected pathway. Two data sets for this scenario are then
#' simulated and denoted as training and test data. This function uses the
#' \href{https://cran.r-project.org/web/packages/Umpire/index.html}{Umpire}
#' package for efficient simulation of multivariate normally distributed
#' random variables.
#'
#' @param mvn Object of class \code{\link[Umpire]{MVN}} for efficient simulation
#' of data from a multivariate normal distribution (Note: either mvn or
#' cov.matrix need to be specified).
#' @param cov.matrix Matrix with covariances, as estimated with the
#' \code{\link[stats]{cov}} function (Note: either mvn or
#' cov.matrix need to be specified).
#' @param gene.names Vector with names of all genes that should be simulated
#' (Note: needs to have the same dimension as cov.matrix).
#' @param no.samples Total number of individuals to simulate. A balanced number
#' of cases and controls will be simulated.
#' @param effects Vector with effect sizes (mean differences) for differentially
#' expressed genes.
#' @param genes.pw Vector with names of genes in selected
#' pathway.
#' @param prop.de Proportion of genes in selected pathway to be differentially
#' expressed.
#' @param seed Random seed. Default is NULL, which generates the seed from R.
#'
#' @import Umpire
#'
#' @return List with the following components:
#' \itemize{
#' \item \code{data.train} data.frame of training data with samples in rows and
#' variables in columns; the first column, called y, denotes the case-control
#' status (0: control, 1: case)
#' \item \code{data.test} data.frame of test data with samples in rows and
#' variables in columns; the first column, called y, denotes the case-control
#' status (0: control, 1: case)
#' \item \code{info.sim} named vector with simulated effect sizes for each of
#' the genes in the selected pathway (Note: list of genes might have been
#' reduced to the genes available in gene.names)
#' }
#'
#' @references
#' Zhang J, Roebuck PL, Coombes KR. (2012) Simulating gene expression data to
#' estimate sample size for class and biomarker discovery. Int J Adv Life Sci.
#' 4:44–51. \url{https://doi.org/10.1186/1471-2105-13-S13-S1}.
#'
#' @examples
#' # set seed to make analysis reproducible
#' set.seed(42)
#'
#' # define effect sizes of differentially expressed genes
#' effects = c(runif(10, -3, -0.5),
#'             runif(10, 0.5, 3))
#'
#' # random covariance matrix for 100 genes
#' no.genes = 100
#' m = matrix(rnorm(no.genes * no.genes, sd = 0.5),
#'            nrow = no.genes)
#' cov.matrix = t(m) %*% m
#'
#' # simulate one replicate with 30 samples and first 20 genes belonging to
#' # selected pathway
#' sim.data.l = sim.data.study.2(cov.matrix = cov.matrix,
#'                               gene.names = paste0("gene", seq_len(no.genes)),
#'                               genes.pw = paste0("gene", 1:20),
#'                               effects = effects,
#'                               prop.de = 0.5,
#'                               seed = 42,
#'                               no.samples = 30)
#'
#' # alternatively use MVN object (might be more efficient if many replicates
#' # are simulated)
#' mvn = Umpire::MVN(mu = rep(0, nrow(cov.matrix)),
#'                   Sigma = cov.matrix, tol = 1e-06)
#' sim.data.l = sim.data.study.2(mvn = mvn,
#'                               gene.names = paste0("gene", seq_len(no.genes)),
#'                               genes.pw = paste0("gene", 1:20),
#'                               effects = effects,
#'                               prop.de = 0.5,
#'                               seed = 42,
#'                               no.samples = 30)
#'
#' @export

sim.data.study.2 = function(mvn = NULL,
                            cov.matrix = NULL,
                            gene.names,
                            no.samples,
                            effects,
                            genes.pw,
                            prop.de,
                            seed = NULL) {

    ## if only covariance matrix is provided
    if (is.null(mvn)) {
        mvn = Umpire::MVN(mu = rep(0, nrow(cov.matrix)),
                          Sigma = cov.matrix, tol = 1e-06)
    }

    if (length(gene.names) != nrow(mvn)) {
        stop("different number of genes in gene.names and mvn!")
    }

    genes.pw.avail = intersect(gene.names, genes.pw)
    if (length(genes.pw.avail) < length(genes.pw)) {
        warning("reduced genes in target pathway to genes available in
                gene.names")
    }
    if (length(genes.pw.avail) < 10) {
        stop("less than 10 genes found in target pathway")
    }

    no.causal = ceiling(prop.de * length(genes.pw.avail))
    set.seed(seed)

    ## randomly select causal genes
    genes.causal = sample(genes.pw.avail,
                          size = no.causal,
                          replace = TRUE)

    ## randomly assign effect sizes
    effects.causal = sample(effects,
                            size = no.causal,
                            replace = FALSE)

    ## training data
    data.train = simulate.data.set(mvn = mvn,
                                   no.samples = no.samples,
                                   gene.names = gene.names,
                                   genes.causal = genes.causal,
                                   effects.causal = effects.causal)

    ## test data
    data.test = simulate.data.set(mvn = mvn,
                                  no.samples = no.samples,
                                  gene.names = gene.names,
                                  genes.causal = genes.causal,
                                  effects.causal = effects.causal)

    ## return values
    info.sim = rep(0, length(genes.pw.avail))
    names(info.sim) = genes.pw.avail
    info.sim[genes.causal] = effects.causal
    return(list(data.train = data.train,
                data.test = data.test,
                info.sim = info.sim))
}
