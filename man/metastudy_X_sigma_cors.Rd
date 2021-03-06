% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specification_tests.R
\name{metastudy_X_sigma_cors}
\alias{metastudy_X_sigma_cors}
\title{Computes correlations to test the independece assumption between
estimate and standard error of Andrews and Kasy (2019)}
\usage{
metastudy_X_sigma_cors(ms)
}
\arguments{
\item{ms}{An object returned from the function \code{metastudies_estimation}}
}
\value{
a data frame with correlations, confidence intervals and also relevent results from a linear regression of standard errors on estimates.
}
\description{
A crucial assumption of Andrews and Kasy (2019) is that in the unobserved
latent distribution without publication error the estimate and its standard
error are statistically independent from each other.
}
\details{
While the latent distribution cannot be observed, this function computes
some correlations that may indicate problems with respect to this assumption.
We use an inverse probability weighting approach.
More precisely, it weights inversely with the estimated publication probabilities to recover the correlation in the
unobserved latent distribution. (This approach was suggested in an email by
Isaiah Andrews and implemented by Sebastian Kranz).

To compute standard errrors via bootstrap (very time consuming),
call the function \code{bootstrap_specification_tests}.
}
