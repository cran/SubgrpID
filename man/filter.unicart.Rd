% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filter.unicart}
\alias{filter.unicart}
\title{filter.unicart}
\usage{
filter.unicart(
  data,
  type,
  yvar,
  xvars,
  censorvar,
  trtvar,
  trtref = 1,
  pre.filter = length(xvars)
)
}
\arguments{
\item{data}{input data frame}

\item{type}{"c" continuous; "s" survival; "b" binary}

\item{yvar}{response variable name}

\item{xvars}{covariates variable name}

\item{censorvar}{censoring variable name 1:event; 0: censor.}

\item{trtvar}{treatment variable name}

\item{trtref}{code for treatment arm}

\item{pre.filter}{NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected}
}
\value{
selected covariates after rpart filtering
}
\description{
rpart filtering
}
\details{
rpart filtering (only for prognostic case)
}
