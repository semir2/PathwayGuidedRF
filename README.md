PathwayGuidedRF
================

This R package provides functions for the identification of important pathways or gene sets using multiple pathway guided random forest (RF) approaches. Furthermore, it includes functions to simulate pathway based gene expression data under two different scenarios.

<!-- Please cite the following manuscript if you use the package:   -->
<!-- S Seifert, S Gundlach, O Junge and S Szymczak (2019) Integrating biological  -->
<!-- knowledge and omics data using pathway guided random forests: a benchmarking  -->
<!-- study. Submitted to Genome Biology. -->
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

    ## P values are now by default adjusted for multiple testing using the Benjamini-Hochberg procedure!

We then specify some characteristics of the three pathways we would like to simulate. Note that the first two pathways contain some differentially expressed genes, while the third pathway has no signal. A well performing pathway analysis approach should thus identify the first two pathways.

``` r
# define pathway parameters
info.par = data.frame(pw = paste0("pw", 1:3),
                      no.genes = c(100, rep(20, 2)),
                      prop.de = c(1, 0.5, 0),
                      cor = rep(0.6, 3),
                      stringsAsFactors = FALSE)
```

We finally simulate 30 individuals, of which 15 should be cases. In order to be able to reproduce the data and the results we set the random number generator to the specified seed. The parameter gamma.range defines the the range of absolute effect sizes for the differential expression.

``` r
set.seed(12345)
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
set.seed(12345)
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

    ##      id no.var     run.1       run.2       run.3       run.4      run.5
    ## pw1 pw1    100  2.389563 2.446646702  2.28263315  2.14467113  2.3111720
    ## pw2 pw2     20  2.359169 2.516721639  2.19719289  2.53212286  2.5810207
    ## pw3 pw3     20 -0.233480 0.002432354 -0.04562368 -0.03925563 -0.1017036
    ##           run.6      run.7       run.8       run.9  decision selected
    ## pw1  2.47805308 2.06059033  2.52851785  2.40617344 Confirmed        1
    ## pw2  2.07525720 2.40633176  2.32963284  2.27268405 Confirmed        1
    ## pw3 -0.04376433 0.02010458 -0.09289833 -0.09862246  Rejected        0

``` r
# hunting approach
set.seed(12345)
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

    ##      id no.var   x.score    z.score       pval  pval.adj selected
    ## pw1 pw1    100 0.1190678  1.4793224 0.06952709 0.2085813        0
    ## pw2 pw2     20 0.1284275  0.6100654 0.27090925 0.4063639        0
    ## pw3 pw3     20 0.0000000 -2.5198624 0.99412996 0.9941300        0

For illustration we use the default of 20 permutations for the prediction error method. However, since we analyse only three pathways the total number of permutations is 60 which is too small for reliable results. In a real analysis the total number of permutation should be larger than 1000 which is often already achieved with 20 permutations per pathway since many more pathways are usually tested.

``` r
# prediction error approach
set.seed(12345)
res.pe = pw.rf.pred.error(x = sim.data.l$data.train[, -1],
                          y = sim.data.l$data.train[, 1],
                          info.pw = sim.data.l$info.pw,
                          type = "classification",
                          not.assoc.pw = FALSE,
                          no.perm = 20)
res.pe$pw.sel
```

    ## [1] "pw1" "pw2"

``` r
res.pe$results.pw
```

    ##      id no.var pred.error pred.error.perm.min pred.error.perm.median
    ## pw1 pw1    100 0.00000000                 0.3              0.5166667
    ## pw2 pw2     20 0.03333333                 0.4              0.5500000
    ## pw3 pw3     20 0.53333333                 0.3              0.5500000
    ##     pred.error.perm.max         pval     pval.adj selected
    ## pw1           0.6666667 5.633729e-07 1.690119e-06        1
    ## pw2           0.7666667 2.534241e-06 3.801361e-06        1
    ## pw3           0.6333333 5.122125e-01 5.122125e-01        0

A warning is issued for the LeFE approach since not enough genes outside of the first pathway are available in our small example data set. Thus, no P value could be calculated for this pathway.

``` r
# LeFE approach
set.seed(12345)
res.lefe = pw.rf.lefe(x = sim.data.l$data.train[, -1],
                      y = sim.data.l$data.train[, 1],
                      info.pw = sim.data.l$info.pw,
                      type = "classification",
                      not.assoc.pw = FALSE)
```

    ## Warning in pw.rf.lefe(x = sim.data.l$data.train[, -1], y =
    ## sim.data.l$data.train[, : too few genes for sampling without category: pw1
    ## (100*6 =< 40)

``` r
res.lefe$pw.sel
```

    ## character(0)

``` r
res.lefe$results.pw
```

    ##      id no.var      pval  pval.adj selected
    ## pw1 pw1    100        NA        NA       NA
    ## pw2 pw2     20 0.7276402 0.9999449        0
    ## pw3 pw3     20 0.9999449 0.9999449        0
