
## internal function for simulation study 2
## simulate omics data set (rows: genes, columns: samples)
simulate.data.set = function(mvn,
                             no.samples,
                             gene.names,
                             genes.causal,
                             effects.causal) {

    x = Umpire::rand(mvn, no.samples)
    rownames(x) = gene.names

    ## modify means in cases
    no.cases = floor(no.samples / 2)
    x[genes.causal, 1:no.cases] = x[genes.causal, 1:no.cases] +
        effects.causal

    ## simulate outcome
    y = c(rep(1, no.cases),
          rep(0, no.samples - no.cases))

    data = data.frame(y = y, t(x))
    return(data)
}

## modified version of function data.gen.DJsim (implementation provided by
## Poisson et al. 2011 BMC Bioinform)
##
## remove metabolites
## simulate two data sets for a single replicate (with fixed mean differences
## and variances per gene)
## removed d0 argument
##
## documentation extracted from original script:
##########################################
### DATA PARAMETERS ###
# n.samp, number of samples
# n.case, number of cases
# n.path, number of pathways
# n.assoc.path, number of pathways associated with case status
# n.gene.path, number of genes per pathway
# d1, probability of differential gene expression within associated pathway
# dd1, probability of up regulation
# gamma.range (a,b), draw gamma from uniform U(a,b), a>=0, b>0
# corr.gg, correlation between two genes
############################################

data.gen.DJsim.mod <- function(n.samp,
                               n.case,
                               n.path,
                               n.assoc.path,
                               n.gene.path,
                               d1,
                               dd1,
                               corr.gg,
                               gamma.range) {

    n.control <- n.samp - n.case
    W.i <- c(rep(1, n.case), rep(0, n.control))

    ## GENE EXPRESSION ##
    n.gene <- n.path * n.gene.path
    #alpha = 0
    #beta.j
    mean.gene.var <- geoR::rinvchisq(n = n.gene,
                                     df = 4)
    beta.j <- stats::rnorm(n = n.gene,
                           mean = 0,
                           sd = 2*sqrt(mean.gene.var))
    #omega.j
    abs.omega.j <- stats::runif(n = n.gene,
                                min = gamma.range[1],
                                max = gamma.range[2])
    bin.omega.j <- stats::rbinom(n = n.gene,
                                 size = 1,
                                 prob = dd1)
    omega.j <- abs.omega.j * ifelse(bin.omega.j == 1,1,-1)
    #Dj
    n.assoc.gene <- n.assoc.path*n.gene.path
    Dj <- stats::rbinom(n = n.gene,
                        size = 1,
                        prob = d1)

    ## Generate two data sets ##
    data.l = NULL
    for (i in 1:2) {
        ## Draw multivariate normal vectors per pathway ##
        y.ij.base <- vector()
        for (l in 1:n.path) {
            start.g <- (l - 1) * n.gene.path + 1
            stop.g <- l * n.gene.path
            # mean vector #
            mean.vector <- c(beta.j[start.g:stop.g])
            # variance-covariance matrix #
            var.g <- mean.gene.var[start.g:stop.g]
            cov.mat <- corr.gg * (sqrt(var.g) %*% t(sqrt(var.g))) #Ng x Ng
            diag(cov.mat) <- c(var.g)
            y.base <- MASS::mvrnorm(n = n.samp,
                                    mu = mean.vector,
                                    Sigma = cov.mat)
            y.ij.base <- rbind(y.ij.base,
                               t(y.base)[1:n.gene.path,])
        }
        y.ij <- y.ij.base + (omega.j*Dj) %*% t(W.i)

        data.l <- c(data.l, list(y.ij))
    }

    return(data.l)
}
