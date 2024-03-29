% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/master_function.R
\name{SubgrpID}
\alias{SubgrpID}
\title{SubgrpID}
\usage{
SubgrpID(
  data.train,
  data.test = NULL,
  yvar,
  censorvar = NULL,
  trtvar = NULL,
  trtref = NULL,
  xvars,
  type = "c",
  n.boot = 25,
  des.res = "larger",
  min.sigp.prcnt = 0.2,
  pre.filter = NULL,
  filter.method = NULL,
  k.fold = 5,
  cv.iter = 20,
  max.iter = 500,
  mc.iter = 20,
  method = c("Seq.BT"),
  do.cv = FALSE,
  out.file = NULL,
  file.path = "",
  plots = FALSE
)
}
\arguments{
\item{data.train}{data frame for training dataset}

\item{data.test}{data frame for testing dataset, default = NULL}

\item{yvar}{variable (column) name for response variable}

\item{censorvar}{variable name for censoring (1: event; 0: censor), default = NULL}

\item{trtvar}{variable name for treatment variable, default = NULL (prognostic signature)}

\item{trtref}{coding (in the column of trtvar) for treatment arm}

\item{xvars}{vector of variable names for predictors (covariates)}

\item{type}{type of response variable: "c" continuous; "s" survival; "b" binary}

\item{n.boot}{number of bootstrap for batting procedure, or the variable selection procedure for PRIM; for PRIM, when n.boot=0, bootstrapping for variable selection is not conducted}

\item{des.res}{the desired response. "larger": prefer larger response. "smaller": prefer smaller response}

\item{min.sigp.prcnt}{desired proportion of signature positive group size for a given cutoff}

\item{pre.filter}{NULL (default), no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected}

\item{filter.method}{NULL (default), no prefiltering; "univariate", univaraite filtering; "glmnet", glmnet filtering; "unicart", univariate rpart filtering for prognostic case}

\item{k.fold}{cross-validation folds}

\item{cv.iter}{Algotithm terminates after cv.iter successful iterations of cross-validation, or after max.iter total iterations, whichever occurs first}

\item{max.iter}{total iterations, whichever occurs first}

\item{mc.iter}{number of iterations for the Monte Carlo procedure to get a stable "best number of predictors"}

\item{method}{current version only supports sequential-BATTing ("Seq.BT") for subgroup identification}

\item{do.cv}{whether to perform cross validation for performance evaluation. TRUE or FALSE (Default)}

\item{out.file}{Name of output result files excluding method name. If NULL no output file would be saved}

\item{file.path}{default: current working directory. When specifying a dir, use "/" at the end. e.g. "TEMP/"}

\item{plots}{default: FALSE. whether to save plots}
}
\value{
A list with SubgrpID output
\describe{
  \item{res}{list of all results from the algorithm}
  \item{train.stat}{list of subgroup statistics on training dataset}
  \item{test.stat}{list of subgroup statistics on testing dataset}
  \item{cv.res}{list of all results from cross-validation on training dataset}
  \item{train.plot}{interaction plot for training dataset}
  \item{test.plot}{interaction plot for testing dataset}
}
}
\description{
Exploratory Subgroup Identification main function
}
\details{
Function for SubgrpID
}
\examples{
# no run
n <- 40
k <- 5
prevalence <- sqrt(0.5)
rho<-0.2
sig2 <- 2
rhos.bt.real <- c(0, rep(0.1, (k-3)))*sig2
y.sig2 <- 1
yvar="y.binary"
xvars=paste("x", c(1:k), sep="")
trtvar="treatment"
prog.eff <- 0.5
effect.size <- 1
a.constent <- effect.size/(2*(1-prevalence))
set.seed(888)
ObsData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
                    sig2=sig2, y.sig2=y.sig2, rho=rho,
                    rhos.bt.real=rhos.bt.real, a.constent=a.constent)
TestData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
                     sig2=sig2, y.sig2=y.sig2, rho=rho,
                     rhos.bt.real=rhos.bt.real, a.constent=a.constent)
subgrp <- SubgrpID(data.train=ObsData$data,
                   data.test=TestData$data,
                   yvar=yvar,
                   trtvar=trtvar,
                   trtref="1",
                   xvars=xvars,
                   type="b",
                   n.boot=5, # suggest n.boot > 25, depends on sample size
                   des.res = "larger",
 #                 do.cv = TRUE,
 #                 cv.iter = 2, # uncomment to run CV
                   method="Seq.BT")
subgrp$res
subgrp$train.stat
subgrp$test.stat
subgrp$train.plot
subgrp$test.plot
#subgrp$cv.res$stats.summary #CV estimates of all results
}
