#' Computes correlations to test the independece assumption between
#' estimate and standard error of Andrews and Kasy (2019)
#'
#' A crucial assumption of Andrews and Kasy (2019) is that in the unobserved
#' latent distribution without publication error the estimate and its standard
#' error are statistically independent from each other.
#'
#' While the latent distribution cannot be observed, this function computes
#' some correlations that may indicate problems with respect to this assumption.
#' Results with \code{mode=="ipv"} use an inverse probability weighting approach.
#' More precisely, it weights inversely with the estimated publication probabilities to recover the correlation in the
#' unobserved latent distribution. (This approach was suggested in an email by
#' Isaiah Andrews and implemented by Sebastian Kranz).
#'
#' To compute standard errrors via bootstrap (very time consuming),
#' call the function \code{bootstrap_specification_tests}.
#'
#' In addition we show the correlations between the estimate and standard
#' errors separately for each interval inside which a constant publication
#' probability is assumed using directly the observed data without
#' inverse probability weighting.
#' While not a formal test, it would be reassuring if these correlations
#' would be close to zero.
#'
#' @param ms An object returned from the function \code{metastudies_estimation}
#' @returns a data frame with correlations, confidence intervals and also relevent results from a linear regression of standard errors on estimates.
metastudy_X_sigma_cors = function(ms, intervals = 1:NROW(ms$prob.df)) {
  #restore.point("metastudy_X_sigma_cors")

  dat = ms$dat

  intervals = intersect(intervals, unique(dat$interval.ind))

  dat.log = dat %>%
    filter(X >0, sigma >0) %>%
    mutate(X = log(X), sigma=log(sigma))

  # Compute inverse probability weighted correlations between X and sigma
  cor.ipv = cov.wt(dat[,1:2],1/dat$pub.prob, cor=TRUE)$cor[1,2]

  cor.log.ipv = cov.wt(dat.log[,1:2],1/dat.log$pub.prob, cor=TRUE)$cor[1,2]

  reg = lm(sigma ~ X, data=dat,weights = dat$pub.prob)
  beta.ipv = coef(reg)[2]
  names(beta.ipv) = NULL
  conf = confint(reg,2,level=0.95)
  conf.ipv.low = conf[1]; conf.ipv.up = conf[2]
  r2.ipv = summary(reg)$r.squared


  reg = lm(sigma ~ X, data=dat.log,weights = dat.log$pub.prob)
  beta.log.ipv = coef(reg)[2]
  names(beta.log.ipv) = NULL
  conf = confint(reg,2,level=0.95)
  conf.log.ipv.low = conf[1]; conf.log.ipv.up = conf[2]
  r2.log.ipv = summary(reg)$r.squared


  # Interval of z-closest to 0

  int.res = bind_rows(lapply(intervals, function(int) {
    dat.int = filter(dat, interval.ind == int)
    dat.log.int = filter(dat.log, interval.ind == int)

    cor.test = cor.test(dat.int[,1], dat.int[,2])
    cor = cor.test$estimate
    conf =  cor.test$conf.int
    conf.cor.low = conf[1]; conf.cor.up = conf[2]

    cor.test = cor.test(dat.log.int[,1], dat.log.int[,2])
    cor.log = cor.test$estimate
    conf =  cor.test$conf.int
    conf.cor.log.low = conf[1]; conf.cor.log.up = conf[2]

    reg = lm(sigma ~ X, data=dat.int)
    beta = coef(reg)[2]
    names(beta) = NULL
    conf = confint(reg,2,level=0.95)
    conf.low = conf[1]
    r2 = summary(reg)$r.squared

    reg = lm(sigma ~ X, data=dat.log.int)
    beta.log = coef(reg)[2]
    names(beta.log) = NULL
    conf.log = confint(reg,2,level=0.95)
    r2.log = summary(reg)$r.squared

    int.label = paste0("[",ms$prob.df$z.min[int],",",ms$prob.df$z.max[int],")")

    bind_rows(
      data.frame(
        mode = int.label, trans = "level",
        cor = cor, conf.cor.low = conf.cor.low, conf.cor.up = conf.cor.up,
        beta = beta, r.sqr = r2,
        conf.beta.low = conf[1], conf.beta.up =conf[2]
      ),
      data.frame(
        mode = int.label, trans = "log",
        cor = cor.log,
        conf.cor.low = conf.cor.log.low, conf.cor.up = conf.cor.log.up,
        beta = beta.log, r.sqr = r2.log,
        conf.beta.low = conf.log[1], conf.beta.up =conf.log[2]
      )
    )
  }))


  ipv.res = data.frame(
    mode = "ipv", trans = "level",
    cor = cor.ipv, conf.cor.low = NA, conf.cor.up = NA,
    beta = beta.ipv, r.sqr = r2.ipv,
    conf.beta.low =  conf.ipv.low, conf.beta.up =  conf.ipv.up
  )
  ipv.log.res = data.frame(
    mode = "ipv", trans = "log",
    cor = cor.log.ipv, conf.cor.low = NA, conf.cor.up = NA,
    beta = beta.log.ipv,  r.sqr = r2.log.ipv,
    conf.beta.low =  conf.log.ipv.low, conf.beta.up = conf.log.ipv.up
  )

  res = bind_rows(
    ipv.res,
    ipv.log.res,
    int.res
  ) %>% as_tibble()
  res
}

#' Call metastudy_X_sigma_cors repeatedly with new bootstrap samples to
#' compute standard errors for correlations using the inverse probability weighting approach
#'
#' Note that this procedure can take very long to run.
#' @param X vector of parameters
#' @param sigma vector of standard errors
#' @param B number of bootstrap replications
#' @param ... parameters passed to \code{metastudies_estimation}
#' @param num.cores number of cores used for parallel computation using \code{mclapply}.
bootstrap_specification_tests = function(X, sigma, B = 10, ..., num.cores = 1) {
  if (num.cores==1) {
    bs = bind_rows(lapply(1:B, function(i) {
      cat("\n",i)
      bootstrap_specification_tests_inner(X,sigma, ...)
    }))
  } else {
    library(parallel)
    bs = bind_rows(mclapply(1:B, mc.cores=num.cores, function(i) {
      bootstrap_specification_tests_inner(X,sigma, ...)
    }))
  }
  sum = summarize_bootstrap_specification_tests(bs)
  list(bs.sim=bs, bs.sum=sum)
}

summarize_bootstrap_specification_tests = function(bs.sim, cols = colnames(bs.sim)[-c(1:2)]) {
  bs = bs.sim

  bind_rows(
    bs %>%
      group_by(mode, trans) %>%
      summarise(
        stat = "ci.low",
        across(cols, ~quantile(.,0.025, na.rm=TRUE))
      ),
    bs %>%
      group_by(mode, trans) %>%
      summarise(
        stat = "ci.up",
        across(cols, ~quantile(.,0.975, na.rm=TRUE))
      ),

    bs %>%
      group_by(mode, trans) %>%
      summarise(
        stat = "median",
        across(cols, ~median(., na.rm=TRUE))
      ),

    bs %>%
      group_by(mode, trans) %>%
      summarise(
        stat = "mean",
        across(cols, ~mean(., na.rm=TRUE))
      )
  ) %>%
  arrange(trans, mode, stat) %>%
  select(trans, mode, stat, cor, beta, r.sqr, everything())
}

bootstrap_specification_tests_inner = function(X, sigma,...) {
  n = NROW(X)
  rows = sample.int(n,n, replace=TRUE)
  X.b = X[rows]
  sigma.b = sigma[rows]
  ms = metastudies_estimation(X.b, sigma.b, ...)
  #restore.point("bootstrap_specification_tests_inner2")
  metastudy_X_sigma_cors(ms)
}
