---
title: "MetaStudies"
output: 
  html_document: 
    keep_md: yes
---



### Overview

This package implements one of the two methods from [Andrews and Kasy (2019)](https://www.aeaweb.org/articles?id=10.1257/aer.20180310) to estimate the degree of publication bias given a meta-study data set that contains estimates and standard errors (typically of reported regression results) from published studies.

The code mainly consists of Maximilan Kasy's code of his shiny app for MetaStudies on publication bias. I slightly modified the interfaces and added some documentation. The original code is here: 

https://github.com/maxkasy/MetaStudiesApp

A crucial (and strong) assumption of the implemented method from Andrews and Kasy (2019) is that in the latent distribution absent publication bias, the estimates and standard errors are statistically independent. This package adds to Maximilian Kasy's original code functions that compute correlations that could indicate possible violations of this assumption. This will be explained further below (see also Kranz and Pütz (2021)).

### Running the app

To run the app with a preloaded example data set just run inside RStudio:


```r
library(MetaStudies)
example.csv = system.file("data/IV.csv", package="MetaStudies")
MetaStudiesApp(example.csv, show.cor=TRUE)
```

You could also put the the code above into an `app.R` file to specify it as a app for Shiny server.

If you want to run the app in Maximilan Kasy's orginal version (without pre-loaded file and no tab for correlations) just type:


```r
library(MetaStudies)
MetaStudiesApp(show.cor=FALSE)
```

### Using the package without the shiny app

One can also use the relevant functions of the package directly without starting the shiny app.


```r
example.csv = system.file("data/IV.csv", package="MetaStudies")
dat = read.csv(example.csv)
head(dat)
```

```
##         X  sigma
## 1  98.592 72.412
## 2 124.907 36.827
## 3   1.070  0.290
## 4  91.604 25.373
## 5   2.150  0.810
## 6   1.050  0.270
```

This is an example data set extracted from the supplemental data of the meta study by [Brodeur et al. (2019)](https://www.aeaweb.org/articles?id=10.1257/aer.20190687). It is a subset of their data set containing estimated coefficients and their standard errors from a sample of economic articles that use an instrumental variable (IV) approach to estimate causal effects.

Let us run an analysies with the Andrews and Kasy (2019) method:

```r
ms = metastudies_estimation(X=dat$X,sigma=dat$sigma,model = "t",cutoffs = c(1.645, 1.96, 2.58),symmetric = TRUE)
# Show plot of results
estimates_plot(ms)
```

![](README_files/figure-html/estimate-1.png)<!-- -->

The upper plot shows the estimated relative publication probabilities for the different intervals of z-statistics specified by the cutoffs. All publication probabilities are relative to the publication probability given z > 2.58 which is normalized to 1.

The lower plot shows the estimated density of the true z-statistics. We can also get the information as a data frame:


```r
ms$est_tab
```

```
##                          µ           t         df [0, 1.645 ] ( 1.645 , 1.96 ]
## estimate       0.019289521 0.012066501 1.83396754  0.20943948        0.7008920
## standard error 0.003529987 0.002766197 0.04715736  0.01479393        0.0488445
##                ( 1.96 , 2.58 ]
## estimate            1.05342030
## standard error      0.05981187
```

The first 3 columns show the parameters of the generalized t-distribution and the next columns the estimated relative publication probabilities compared to z > 2.58.

The field `dat` contains a data frame that adds as columns the corresponding publication probability `pub.prob` and an index of the interval into which z falls:

```r
head(ms$dat)
```

```
##         X  sigma        z interval.ind  pub.prob
## 1  98.592 72.412 1.361542            1 0.2094395
## 2 124.907 36.827 3.391723            4 1.0000000
## 3   1.070  0.290 3.689655            4 1.0000000
## 4  91.604 25.373 3.610294            4 1.0000000
## 5   2.150  0.810 2.654321            4 1.0000000
## 6   1.050  0.270 3.888889            4 1.0000000
```

## Testing the independence assumption of Andrews and Kasy (2019)

A crucial (and strong) assumption of the implemented method from Andrews and Kasy (2019) is that in the latent distribution absent publication bias, the estimates (`X`) and standard errors (`sigma`) are statistically independent.

Of course, the problem is that we don't observe the latent distribution of `X` and `sigma` absent publication bias but only the resulting distribution distorted by publication bias.

As a first informal test, we can look at the correlations between `X` and `sigma` inside each interval of z-statistics for which we assume a constant publication probability. E.g. here for all non-significant observations with `abs(z)< 1.645`:


```r
fdat = ms$dat %>%
  filter(abs(z) <= 1.645, abs(sigma)>0, abs(X)>0)

cor.test(log(fdat$X), log(fdat$sigma))
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  log(fdat$X) and log(fdat$sigma)
## t = 94.053, df = 1582, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.9132008 0.9281756
## sample estimates:
##       cor 
## 0.9210279
```

We find that in the subsample of non-significant z-statistics the log's of estimated coefficients and standard errors are highly correlated with a 95% CI of [0.913, 0.928].

Kranz and Pütz (2021) state as intuition for that correlation that a lot of variation in estimated coefficient and standard error is driven by the scaling of the dependent and explanatory variables in the regressions. For example, if we would rescale all exlplanatory variables in a regression by multiplying them with a factor `1/m`, then the corresponding estimates and standard errors would both change by the factor `m`. This could lead to a strong positive correlation for the log-values even absent publication bias.

Yet, as Isaiah Andrews has pointed out in an email communication from a theoretic perspective that informal sub-sample test is less suited as a formal specification test since it has the theoretical problem that when taking a sub-sample conditional on z, this may in principle induce a statistical dependence between `X` and `sigma` even if both variables are statistically independent in the complete sample. 

Isaiah Andrews suggested as an alternative, to use all observations but weight them with the inverse of the estimated publication probabilities. Under the null hypothesis that all assumptions of the chosen specification of the Andrews and Kasy (2019) approach are satisfied, this inverse probability weighting allows to recover the correlation in the unobserved latent distribution of tests if no publication bias were present. Let's do it:


```r
dat.log = ms$dat %>%
  filter(X >0, sigma >0) %>%
  mutate(X = log(X), sigma=log(sigma))

cor.ipv = cov.wt(dat.log[,1:2],1/dat.log$pub.prob, cor=TRUE)$cor[1,2]
cor.ipv
```

```
## [1] 0.8782808
```

Also the inverse probability weight shows a high correlation between both variables.

The function `metastudy_X_sigma_cors` automatically computes those correlations and additional measures:


```r
metastudy_X_sigma_cors(ms)
```

```
## # A tibble: 10 x 9
##    mode         trans   cor conf.cor.low conf.cor.up  beta r.sqr conf.beta.low
##    <chr>        <chr> <dbl>        <dbl>       <dbl> <dbl> <dbl>         <dbl>
##  1 ipv          level 0.909       NA          NA     0.866 0.808         0.854
##  2 ipv          log   0.878       NA          NA     0.904 0.840         0.893
##  3 [0,1.645)    level 0.909        0.900       0.917 0.919 0.827         0.899
##  4 [0,1.645)    log   0.921        0.913       0.928 0.853 0.848         0.836
##  5 [1.645,1.96) level 1.00         1.00        1.00  0.598 1.00          0.597
##  6 [1.645,1.96) log   1.00         1.00        1.00  1.00  0.999         0.998
##  7 [1.96,2.58)  level 1.00         0.999       1.00  0.449 0.999         0.448
##  8 [1.96,2.58)  log   0.999        0.999       0.999 0.998 0.998         0.995
##  9 [2.58,Inf)   level 0.994        0.994       0.995 0.355 0.989         0.354
## 10 [2.58,Inf)   log   0.955        0.951       0.959 0.988 0.912         0.974
## # ... with 1 more variable: conf.beta.up <dbl>
```

You see the correlations in log and levels and also the coefficient of a linear regression of `sigma` on `X`.

### Confidence intervals for the inverse probability weighting approach

Note that correct standard errors and confidence intervals for the inverse probability weighting approach are not trivial to compute since they must also take into account the uncertainity when estimating the publication probabilities. Iasiah Andrews suggested that a GMM approach could yield correct standard errors. However, not beeing a theoretical econometrician, I don't find it trivial to implement the GMM approach. Luckiliy, one can oten substitute theoretical understanding with brute-force bootstrap computation. 

The function `bootstrap_specification_tests` computes confidence intervals using a bootstrap approach. In each run, we draw a bootstrap sample from the original data, run the Andrews and Kasy approach on the bootstrap sample and then compute the correlation using the estimated publication probabilities. This procedure can be very time consuming and is perhaps better run on a server with multiple cores:


```r
# Takes long to run unless you have a lot of cores
bst = bootstrap_specification_tests(dat$X, dat$sigma, B = 100, num.cores = 10)
```

The resulting object contains a field `sim` that contains raw data of all bootstrap results and a field `sum` that contains summary results.

I have included the bootstrap summary for the example data set in the MetaStudies package: 

```r
sum = readRDS(system.file("data/bootstrap_sum_IV.Rds", package="MetaStudies"))
filter(sum, trans=="log", mode=="ipv") %>%
  select(trans, mode, stat, cor)
```

```
## # A tibble: 4 x 4
## # Groups:   mode [1]
##   trans mode  stat     cor
##   <fct> <fct> <chr>  <dbl>
## 1 log   ipv   ci.low 0.862
## 2 log   ipv   ci.up  0.891
## 3 log   ipv   mean   0.877
## 4 log   ipv   median 0.878
```

We see from the rows where `stat` is `ci.low` and `ci.up` that the 95% CI of the correlation computed via inverse probabilty weighting is far away from 0. This suggests that the independence assumption of Andrews and Kasy (2019) seems violated in this example data set.

### References

- Andrews, Isaiah and Maximilian Kasy. 2019. “Identification of and correction for publication bias.” American Economic Review 109 (8): 2766-94.

- Brodeur, Abel, Nikolai Cook, and Anthony Heyes. 2020. “Methods Matter: p-Hacking and Publication Bias in Causal Analysis in Economics.” American Economic Review, 110 (11): 3634-60.

- Kranz, Sebastian and Peter Pütz. 2021 “Rounding and other pitfalls in meta-studies on p-hacking and publication bias. A comment on Brodeur et al. (2020)”, working paper.
