PathwayGuidedRF
================

This R package provides functions for the identification of important pathways or gene sets using multiple pathway guided random forest (RF) approaches. Furthermore, it includes functions to simulate pathway based gene expression data under two different scenarios.

Please cite the following manuscript if you use the package:
S Seifert, S Gundlach, O Junge and S Szymczak (2019) Integrating biological knowledge and omics data using pathway guided random forests: a benchmarking study. Submitted to Genome Biology.

Installation
------------

<!-- install.packages(c("ranger", "Boruta", "Umpire", "geoR", "MASS")) -->
You can install the most recent version of PathwayGuidedRF from GitHub as follows:

``` r
library(devtools)
install_github("szymczak-lab/PathwayGuidedRF")
```

Simulated data
--------------

We will use an example data set simulated under the scenario of simulation study 1 as described in the paper.

First, we load the package:

``` r
library(PathwayGuidedRF)
```

We then specify some characteristics of the three pathways we would like to simulate. Note that the first two pathways contain some differentially expressed genes, while the third pathway has no signal. A well performing pathway analysis approach should thus identify the first two pathways.

``` r
# define pathway parameters
info.par = data.frame(pw = paste0("pw", 1:3),
                      no.genes = c(100, rep(20, 2)),
                      prop.de = c(1, 0.5, 0),
                      cor = rep(0.6, 3),
                      stringsAsFactors = FALSE)
```

We finally simulate 30 individuals, of which 15 should be cases. The parameter gamma.range defines the the range of absolute effect sizes for the differential expression.

``` r
sim.data.l = sim.data.study.1(info.par = info.par,
                              no.samples = 30,
                              no.cases = 15,
                              seed = 42,
                              gamma.range = c(0.5, 1.5))
```

The resulting list sim.data.l contains a training and a test data set simulated under the same setting (i.e. same genes and effect sizes). We can check the dimensions (30 individuals and 140 genes). Note that the first column contains the outcome coded as 0 (= control) and 1 (= case).

``` r
data.train = sim.data.l$data.train
dim(data.train)
```

    ## [1]  30 141

``` r
head(data.train[, 1:5])
```

    ##      y     gene1      gene2      gene3      gene4
    ## ind1 1 -1.768058 -2.8195900 -1.8894840  0.3253353
    ## ind2 1 -1.424678 -1.1272396 -0.9619253  0.7730726
    ## ind3 1 -1.769282 -1.8992853 -1.7081115  2.2135723
    ## ind4 1 -1.222623 -0.8822824 -1.4896662  1.9478316
    ## ind5 1 -2.343630 -2.8215748 -1.6992956 -0.1723490
    ## ind6 1 -1.777972 -2.3467923 -1.5735073 -0.4089688

``` r
table(data.train$y)
```

    ## 
    ##  0  1 
    ## 15 15

Pathway analysis using RF
-------------------------

We can perform RF based pathway analysis using several approaches that are implemented in the package. The results have a similar structure but also contain some method specific information.

``` r
# synthetic feature approach
res.sf = pw.rf.synthetic.features(x = sim.data.l$data.train[, -1],
                                  y = sim.data.l$data.train[, 1],
                                  info.pw = sim.data.l$info.pw,
                                  type = "classification",
                                  not.assoc.pw = FALSE)

# selected pathways
res.sf$pw.sel
```

    ## [1] "pw1" "pw2"

``` r
# more information for each pathway
res.sf$results.pw
```

    ##      id no.var       run.1        run.2       run.3    run.4       run.5
    ## pw1 pw1    100  2.48109919  2.375730425  2.29551243 2.324955  2.42981470
    ## pw2 pw2     20  2.23084332  2.561339322  2.37880148 2.636408  2.31525519
    ## pw3 pw3     20 -0.00690834 -0.005310959 -0.06775715 0.044172 -0.02914315
    ##          run.6      run.7      run.8      run.9    run.10     run.11
    ## pw1  2.5302109  2.2607940 2.31483019 2.36222675 2.4729116 1.99022825
    ## pw2  2.5049536  2.2452267 2.27746236 2.80113734 2.1682347 2.43528680
    ## pw3 -0.0889709 -0.4628929 0.03927959 0.02257604 0.0612325 0.01819423
    ##           run.12  decision selected
    ## pw1  2.480991105 Confirmed        1
    ## pw2  2.254337086 Confirmed        1
    ## pw3 -0.007816118  Rejected        0

``` r
# hunting approach
res.hunt = pw.rf.hunt(x = sim.data.l$data.train[, -1],
                      y = sim.data.l$data.train[, 1],
                      info.pw = sim.data.l$info.pw,
                      type = "classification",
                      not.assoc.pw = FALSE)
res.hunt$pw.sel
```

    ## character(0)

``` r
res.hunt$results.pw
```

    ##      id no.var      x.score    z.score       pval  pval.adj selected
    ## pw1 pw1    100 0.1194429386  1.5803931 0.05700844 0.1710253        0
    ## pw2 pw2     20 0.1219263095  0.4687985 0.31960685 0.9588206        0
    ## pw3 pw3     20 0.0001805263 -2.5090771 0.99394765 1.0000000        0

For reduced run time we use only a small number of permutations for the prediction error method. However, the default value of 1000 should be used for a real analysis.

``` r
# prediction error approach
res.pe = pw.rf.pred.error(x = sim.data.l$data.train[, -1],
                                  y = sim.data.l$data.train[, 1],
                                  info.pw = sim.data.l$info.pw,
                                  type = "classification",
                                  not.assoc.pw = FALSE,
                                  no.perm = 100)
res.pe$pw.sel
```

    ## [1] "pw1" "pw2"

``` r
res.pe$results.pw
```

    ##      id no.var pred.error pred.error.perm.min pred.error.perm.median
    ## pw1 pw1    100 0.00000000           0.2666667              0.5666667
    ## pw2 pw2     20 0.03333333           0.3000000              0.5666667
    ## pw3 pw3     20 0.63333333           0.3333333              0.5333333
    ##     pred.error.perm.max no.perm.used       pval   pval.adj selected
    ## pw1           0.8333333          100 0.00990099 0.02970297        1
    ## pw2           0.7666667          100 0.00990099 0.02970297        1
    ## pw3           0.7666667          100 0.77227723 1.00000000        0

A warning is issued for the LeFE approach since not enough genes outside of the first pathway are available in our small example data set. Thus, no P value could be calculated for this pathway.

``` r
# LeFE approach
res.lefe = pw.rf.lefe(x = sim.data.l$data.train[, -1],
                      y = sim.data.l$data.train[, 1],
                      info.pw = sim.data.l$info.pw,
                      type = "classification",
                      not.assoc.pw = FALSE)
```

    ## Warning in pw.rf.lefe(x = sim.data.l$data.train[, -1], y = sim.data.l
    ## $data.train[, : too few genes for sampling without category: pw1 (100*6 =<
    ## 40)

``` r
res.lefe$pw.sel
```

    ## character(0)

``` r
res.lefe$results.pw
```

    ##      id no.var      pval pval.adj selected
    ## pw1 pw1    100        NA       NA       NA
    ## pw2 pw2     20 0.7360294        1        0
    ## pw3 pw3     20 0.9999367        1        0
