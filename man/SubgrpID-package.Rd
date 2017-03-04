\name{SubgrpID-package}
\alias{SubgrpID-package}
\docType{package}
\title{
Patient subgroup identification for clinical drug development
}
\description{
Prognostic and predictive biomarker signature development for Exploratory Subgroup Identification in Randomized Clinical Trials
}
\details{
\tabular{ll}{
Package: \tab SubgrpID\cr
Type: \tab Package\cr
Version: \tab 0.10\cr
Date: \tab 2017-01-25\cr
License: \tab GPL-3\cr
}
}
\author{
Xin Huang, Yan Sun, Saptarshi Chatterjee and Paul Trow
Maintainer: Xin Huang <xhuang.stats@gmail.com>
}
\references{
  Huang X. et al. (2017)
  Patient subgroup identification for clinical drug development. \emph{Statistics in Medicine}, doi: 10.1002/sim.7236.
  
  Chen G. et al. (2015)
  A PRIM approach to predictive-signature development for patient stratification
  \emph{Statistics in Medicine}, \bold{34}, 317-342.
}
\keyword{ AIM }
%\seealso{
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
                      yvar=yvar,
                      trtvar=trtvar,
                      trtref="active",
                      xvars=xvars,
                      type="b",
                      des.res = "smaller",
                      method="AIM.Rule")
    subgrp$res
    subgrp$train.stat
    subgrp$train.plot
  }
}
