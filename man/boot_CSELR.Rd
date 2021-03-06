% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatCompR.R
\name{boot_CSELR}
\alias{boot_CSELR}
\alias{boot.selr}
\title{Use bootstrap method to compute the p_vaule of Composite semiparametric empirical likelihood ratio}
\usage{
boot.selr(dat, B = 999, n.cores = 1)
}
\arguments{
\item{dat}{the data}

\item{B}{is the number of bootstrap replication, with default B=999}

\item{n.cores}{is the number of cores used for parelell computing, with default n.cores=1}
}
\value{
the p_value of Composite semiparametric empirical likelihood ratio
}
\description{
use bootstrap method to compute the p value
}
\examples{
\dontrun{
bs<-boot.selr(data)
}
}
