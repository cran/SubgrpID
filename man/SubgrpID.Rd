\name{SubgrpID}
\alias{SubgrpID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exploratory Subgroup Identification main function
}
\description{
Prognostic and predictive biomarker signature development for Exploratory Subgroup Identification in Randomized Clinical Trials
}
\usage{
SubgrpID(data.train, data.test=NULL, 
         yvar, 
         censorvar=NULL, 
         trtvar=NULL, 
         trtref=NULL, 
         xvars, 
         type="c",
         n.boot=ifelse(method=="PRIM",0,25),
         des.res="larger",
         min.sigp.prcnt=0.20, 
         pre.filter=NULL, 
         filter.method=NULL, 
         k.fold=5, 
         cv.iter=20, 
         max.iter=500, 
         mc.iter=20, 
         method=c("AIM.Rule"), 
         train.percent.prim=0.5, 
         do.cv=FALSE, 
         out.file=NULL, 
         file.path="", 
         plots=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.train}{
data frame for training dataset
}
  \item{data.test}{
data frame for testing dataset, default = NULL
}
  \item{yvar}{
variable (column) name for response variable
}
  \item{censorvar}{
variable name for censoring (1: event; 0: censor), default = NULL
}
  \item{trtvar}{
variable name for treatment variable, default = NULL (prognostic signature)
}
  \item{trtref}{
coding (in the column of trtvar) for treatment arm
}
  \item{xvars}{
vector of variable names for predictors (covariates)
}
  \item{type}{
type of response variable: "c" continuous; "s" survival; "b" binary
}
  \item{n.boot}{
number of bootstrap for batting procedure, or the variable selection procedure for PRIM; for PRIM, when n.boot=0, bootstrapping for variable selection is not conducted
}
  \item{des.res}{
the desired response. "larger": prefer larger response. "smaller": prefer smaller response
}
  \item{min.sigp.prcnt}{
desired proportion of signature positive group size for a given cutoff
}
  \item{pre.filter}{
NULL (default), no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
}
  \item{filter.method}{
NULL (default), no prefiltering; "univariate", univaraite filtering; "glmnet", glmnet filtering; "unicart", univariate rpart filtering for prognostic case
}
  \item{k.fold}{
cross-validation folds
}
  \item{cv.iter}{
Algotithm terminates after cv.iter successful iterations of cross-validation, or after max.iter total iterations, whichever occurs first
}
  \item{max.iter}{
total iterations, whichever occurs first
}
  \item{mc.iter}{
number of iterations for the Monte Carlo procedure to get a stable "best number of predictors"
}
  \item{method}{
algorithms performed for subgroup identification, one of the ("AIM", "AIM.Rule", "Seq.BT", "PRIM")
}
  \item{train.percent.prim}{
percentage of the sub-training set used only by PRIM method; if train.percent.prim=1, all data will be used both for sub-training and sub-testing inside the PRIM 
  }
  \item{do.cv}{
whether to perform cross validation for performance evaluation. TRUE or FALSE (Default)
}
  \item{out.file}{
Name of output result files excluding method name. If NULL no output file would be saved
}
  \item{file.path}{
default: current working directory. When specifying a dir, use "/" at the end. e.g. "TEMP/"
}
  \item{plots}{
default: F. whether to save plots
}
}
\details{ 
The function contains four algorithms for developing threshold-based multivariate (prognostic/predictive) biomarker signatures via resampled tree-based algorithms (Sequential BATTing),  Monte-Carlo variations of the Adaptive Indexing method (AIM and AIM-Rule)and  Patient Rule Induction Method. 
Variable selection is automatically built-in to these algorithms. Final signatures are returned with interaction plots for predictive signatures. Cross-validation performance evaluation and testing dataset results are also output.
}
\value{

  \item{res}{list of all results from the algorithm}
  \item{train.stat}{list of subgroup statistics on training dataset}
  \item{test.stat}{list of subgroup statistics on testing dataset}
  \item{cv.res}{list of all results from cross-validation on training dataset}
  \item{train.plot}{interaction plot for training dataset}
  \item{test.plot}{interaction plot for testing dataset}
}
\references{
  Huang X. et al. (2017)
  Patient subgroup identification for clinical drug development. \emph{Statistics in Medicine}, doi: 10.1002/sim.7236.

  Chen G. et al. (2015)
  A PRIM approach to predictive-signature development for patient stratification
  \emph{Statistics in Medicine}, \bold{34}, 317-342.
}
\author{
Xin Huang, Yan Sun, Saptarshi Chatterjee and Paul Trow
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  \dontrun{
    data(Sepsis.train)
    data(Sepsis.test)
    
    yvar="survival"
    xvars=names(Sepsis.train)[2:12]
    trtvar="THERAPY"
    
    set.seed(123)
    subgrp <- SubgrpID(data.train=Sepsis.train,
                      data.test=Sepsis.test,
                      yvar=yvar,
                      trtvar=trtvar,
                      trtref="active",
                      xvars=xvars,
                      type="b",
                      des.res = "smaller",
                      method="AIM.Rule")
    subgrp$res
    subgrp$train.stat
    subgrp$test.stat
    subgrp$train.plot
    subgrp$test.plot
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BATTing }
\keyword{ Sequential-BATTing }
\keyword{ AIM-Rule }
\keyword{ PRIM }
