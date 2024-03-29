% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross_validation_functions.R
\name{kfold.cv}
\alias{kfold.cv}
\title{kfold.cv}
\usage{
kfold.cv(
  data,
  model.Rfunc,
  model.Rfunc.args,
  predict.Rfunc,
  predict.Rfunc.args,
  k.fold = 5,
  cv.iter = 50,
  strata,
  max.iter = 500
)
}
\arguments{
\item{data}{the CV data}

\item{model.Rfunc}{Name of the model function.}

\item{model.Rfunc.args}{List of input arguments to model.Rfunc.}

\item{predict.Rfunc}{Name of the prediction function, which takes the prediction rule returned by model.Rfunc along with any input data (not necessarily the input data to kfold.cv) and returns a TRUE-FALSE predictionvector specifying the positive and negative classes for the data.}

\item{predict.Rfunc.args}{List containing input arguments to predict.Rfunc, except for data and predict.rule.}

\item{k.fold}{Number of folds of the cross-validation.}

\item{cv.iter}{Number of iterations of the cross-validation. If model.Rfunc returns an error at any of the k.fold calls, the current iteration is aborted. Iterations are repeated until cv.iter successful iterations have occurred.}

\item{strata}{Stratification vector of length the number of rows of data, usually corresponding to the vector of events.}

\item{max.iter}{Function stops after max.iter iterations even if cv.iter successful iterations have not occurred.}
}
\value{
List of length 2 with the following fields:

cv.data - List of length cv.iter. Entry i contains the output of predict.Rfunc at the ith iteration.

sig.list - list of length cv.iter * k.fold, whose entries are the prediction.rules (signatures) returned by model.Rfunc at each k.fold iteration.
}
\description{
Perform k-fold cross-validation of a model.
}
\details{
Perform k-fold cross-validation of a model.
}
