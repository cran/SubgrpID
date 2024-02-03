##################################################################################

######### !!! source required files before running the algorithms !!!##############

##################################################################################

#' Function for SubgrpID
#' @title SubgrpID
#' @description Exploratory Subgroup Identification main function
#' @param data.train data frame for training dataset
#' @param data.test data frame for testing dataset, default = NULL
#' @param yvar variable (column) name for response variable
#' @param censorvar variable name for censoring (1: event; 0: censor), default = NULL
#' @param trtvar variable name for treatment variable, default = NULL (prognostic signature)
#' @param trtref coding (in the column of trtvar) for treatment arm
#' @param xvars vector of variable names for predictors (covariates)
#' @param type type of response variable: "c" continuous; "s" survival; "b" binary
#' @param n.boot number of bootstrap for batting procedure, or the variable selection procedure for PRIM; for PRIM, when n.boot=0, bootstrapping for variable selection is not conducted
#' @param des.res the desired response. "larger": prefer larger response. "smaller": prefer smaller response
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff
#' @param pre.filter NULL (default), no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#' @param filter.method NULL (default), no prefiltering; "univariate", univaraite filtering; "glmnet", glmnet filtering; "unicart", univariate rpart filtering for prognostic case
#' @param k.fold cross-validation folds
#' @param cv.iter Algotithm terminates after cv.iter successful iterations of cross-validation, or after max.iter total iterations, whichever occurs first
#' @param max.iter total iterations, whichever occurs first
#' @param mc.iter number of iterations for the Monte Carlo procedure to get a stable "best number of predictors"
#' @param method current version only supports sequential-BATTing ("Seq.BT") for subgroup identification
#' @param do.cv whether to perform cross validation for performance evaluation. TRUE or FALSE (Default)
#' @param out.file Name of output result files excluding method name. If NULL no output file would be saved
#' @param file.path default: current working directory. When specifying a dir, use "/" at the end. e.g. "TEMP/"
#' @param plots default: FALSE. whether to save plots
#' @return A list with SubgrpID output
#' \describe{
#'   \item{res}{list of all results from the algorithm}
#'   \item{train.stat}{list of subgroup statistics on training dataset}
#'   \item{test.stat}{list of subgroup statistics on testing dataset}
#'   \item{cv.res}{list of all results from cross-validation on training dataset}
#'   \item{train.plot}{interaction plot for training dataset}
#'   \item{test.plot}{interaction plot for testing dataset}
#' }
#' @examples
#' # no run
#' n <- 40
#' k <- 5
#' prevalence <- sqrt(0.5)
#' rho<-0.2
#' sig2 <- 2
#' rhos.bt.real <- c(0, rep(0.1, (k-3)))*sig2
#' y.sig2 <- 1
#' yvar="y.binary"
#' xvars=paste("x", c(1:k), sep="")
#' trtvar="treatment"
#' prog.eff <- 0.5
#' effect.size <- 1
#' a.constent <- effect.size/(2*(1-prevalence))
#' set.seed(888)
#' ObsData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
#'                     sig2=sig2, y.sig2=y.sig2, rho=rho,
#'                     rhos.bt.real=rhos.bt.real, a.constent=a.constent)
#' TestData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
#'                      sig2=sig2, y.sig2=y.sig2, rho=rho,
#'                      rhos.bt.real=rhos.bt.real, a.constent=a.constent)
#' subgrp <- SubgrpID(data.train=ObsData$data,
#'                    data.test=TestData$data,
#'                    yvar=yvar,
#'                    trtvar=trtvar,
#'                    trtref="1",
#'                    xvars=xvars,
#'                    type="b",
#'                    n.boot=5, # suggest n.boot > 25, depends on sample size
#'                    des.res = "larger",
#'  #                 do.cv = TRUE,
#'  #                 cv.iter = 2, # uncomment to run CV
#'                    method="Seq.BT")
#' subgrp$res
#' subgrp$train.stat
#' subgrp$test.stat
#' subgrp$train.plot
#' subgrp$test.plot
#' #subgrp$cv.res$stats.summary #CV estimates of all results
#' @export
SubgrpID = function(data.train, data.test=NULL, yvar, censorvar=NULL, trtvar=NULL, trtref=NULL,
                    xvars, type="c",n.boot=25, des.res="larger",
                    min.sigp.prcnt=0.20, pre.filter=NULL, filter.method=NULL,k.fold=5, cv.iter=20,
                    max.iter=500, mc.iter=20, method=c("Seq.BT"),
                    do.cv=FALSE, out.file=NULL, file.path="", plots=FALSE)
{


  if(is.null(censorvar)) {
    censor.vec.train=NULL
    censor.vec.test=NULL
  }else{
    censor.vec.train=data.train[,censorvar,drop=FALSE]
    censor.vec.test=data.test[,censorvar,drop=FALSE]
  }
  if(is.null(trtvar)) {
    trt.vec.train=NULL
    trt.vec.test=NULL
  }else{
    trt.vec.train=data.train[,trtvar,drop=FALSE]
    trt.vec.test=data.test[,trtvar,drop=FALSE]
  }


  if(method=="Seq.BT"){
    res=seqlr.batting(y=data.train[,yvar,drop=FALSE], x=data.train[,xvars,drop=FALSE], censor.vec=censor.vec.train,
                      trt.vec=trt.vec.train, trtref=trtref, type=type, n.boot=n.boot,
                      des.res=des.res, min.sigp.prcnt=min.sigp.prcnt,
                      pre.filter=pre.filter, filter.method=filter.method)
    pred.data=pred.seqlr(data.train[,xvars,drop=FALSE], res)
    train.stat=evaluate.results(data.train[,yvar,drop=FALSE], pred.data, censor.vec=censor.vec.train,
                                trt.vec=trt.vec.train, trtref=trtref, type=type)
    if (!is.null(data.test)){
      pred.data=pred.seqlr(data.test[,xvars,drop=FALSE], res)
      test.stat=evaluate.results(data.test[,yvar,drop=FALSE], pred.data, censor.vec=censor.vec.test,
                                 trt.vec=trt.vec.test, trtref=trtref, type=type)
    }

    if (do.cv){
      cv.res=cv.seqlr.batting(y=data.train[,yvar,drop=FALSE],x=data.train[,xvars,drop=FALSE], censor.vec=censor.vec.train,
                              trt.vec=trt.vec.train, trtref=trtref, type=type, n.boot=n.boot,
                              des.res=des.res,  min.sigp.prcnt=min.sigp.prcnt,
                              pre.filter=pre.filter, filter.method=filter.method,
                              k.fold=k.fold, cv.iter=cv.iter, max.iter=max.iter)
    }


    if (!is.null(trtvar)){
      train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (Seq.BT, Train)", trt.lab=c("Trt.", "Ctrl."))
      if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (Seq.BT, Test)", trt.lab=c("Trt.", "Ctrl."))
    }
  }

  #############################################################################################
  # if (method == "PRIM") {
  #
  #   res = prim.train(data=data.train, yvar=yvar, censorvar=censorvar, trtvar=trtvar,
  #                    trtref=trtref, xvars=xvars, type=type, des.res=des.res,
  #                    alpha = c(0.05,0.06,0.07,0.08,0.09,0.10,0.20,0.30,0.40,0.50),
  #                    min.sigp.prcnt=min.sigp.prcnt,  training.percent = train.percent.prim,
  #                    n.boot=n.boot, pre.filter=pre.filter, filter.method=filter.method)
  #   pred.data=pred.prim(data.train[,xvars,drop=FALSE], res)
  #   train.stat=evaluate.results(data.train[,yvar,drop=FALSE], pred.data, censor.vec=censor.vec.train,
  #                               trt.vec=trt.vec.train, trtref=trtref, type=type)
  #
  #   if (!is.null(data.test)){
  #     pred.data=pred.prim(data.test[,xvars,drop=FALSE], res)
  #     test.stat=evaluate.results(data.test[,yvar,drop=FALSE], pred.data, censor.vec=censor.vec.test,
  #                                trt.vec=trt.vec.test, trtref=trtref, type=type)
  #   }
  #   if (do.cv){
  #     cv.res=prim.cv(data=data.train, yvar=yvar, censorvar=censorvar, trtvar=trtvar,
  #                    trtref=trtref, xvars=xvars, type=type, des.res=des.res,
  #                    alpha = c(0.05,0.06,0.07,0.08,0.09,0.10,0.20,0.30,0.40,0.50),
  #                    min.sigp.prcnt=min.sigp.prcnt,  training.percent = train.percent.prim,
  #                    n.boot=n.boot, pre.filter=pre.filter, filter.method=filter.method,
  #                    k.fold=k.fold, cv.iter=cv.iter, max.iter=max.iter)
  #
  #   }
  #   if (!is.null(trtvar)){
  #     train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (PRIM, Train)", trt.lab=c("Trt.", "Ctrl."))
  #     if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (PRIM, Test)", trt.lab=c("Trt.", "Ctrl."))
  #   }
  #
  # }
  #

  #############################################################################################
  # if(method=="AIM"){
  #   #        library(AIM)
  #   res=aim.batting(y=data.train[,yvar,drop=FALSE],x=data.train[,xvars,drop=FALSE], censor.vec=censor.vec.train,
  #                   trt.vec=trt.vec.train, trtref=trtref, type=type, n.boot=n.boot,
  #                   des.res=des.res, min.sigp.prcnt=min.sigp.prcnt, mc.iter=mc.iter,
  #                   pre.filter=pre.filter, filter.method=filter.method)
  #   pred.data=pred.aim(data.train[,xvars,drop=FALSE], res)
  #   train.stat=evaluate.results(data.train[,yvar,drop=FALSE], pred.data, censor.vec=censor.vec.train,
  #                               trt.vec=trt.vec.train, trtref=trtref, type=type)
  #
  #   if (!is.null(data.test)){
  #     pred.data=pred.aim(data.test[,xvars,drop=FALSE], res)
  #     test.stat=evaluate.results(data.test[,yvar,drop=FALSE], pred.data, censor.vec=censor.vec.test,
  #                                trt.vec=trt.vec.test, trtref=trtref, type=type)
  #   }
  #   if (do.cv){
  #     cv.res=cv.aim.batting(y=data.train[,yvar,drop=FALSE],x=data.train[,xvars,drop=FALSE], censor.vec=censor.vec.train,
  #                           trt.vec=trt.vec.train, trtref=trtref, type=type,
  #                           n.boot=n.boot, des.res=des.res, min.sigp.prcnt=min.sigp.prcnt,
  #                           mc.iter=mc.iter, pre.filter=pre.filter, filter.method=filter.method,
  #                           k.fold=k.fold, cv.iter=cv.iter, max.iter=max.iter)
  #
  #   }
  #
  #   if (!is.null(trtvar)){
  #     train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (AIM, Train)", trt.lab=c("Trt.", "Ctrl."))
  #     if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (AIM, Test)", trt.lab=c("Trt.", "Ctrl."))
  #   }
  #
  #
  # }
  #
  #
  # #############################################################################################
  # if(method=="AIM.Rule"){
  #   #        library(AIM)
  #   res=aim.rule.batting(y=data.train[,yvar,drop=FALSE], x=data.train[,xvars,drop=FALSE], censor.vec=censor.vec.train,
  #                        trt.vec=trt.vec.train, trtref=trtref, type=type, n.boot=n.boot,
  #                        des.res=des.res, min.sigp.prcnt=min.sigp.prcnt, mc.iter=mc.iter,
  #                        pre.filter=pre.filter, filter.method=filter.method)
  #
  #   pred.data=pred.aim(data.train[,xvars,drop=FALSE], res)
  #   train.stat=evaluate.results(data.train[,yvar,drop=FALSE], pred.data, censor.vec=censor.vec.train,
  #                               trt.vec=trt.vec.train, trtref=trtref, type=type)
  #
  #   if (!is.null(data.test)){
  #     pred.data=pred.aim(data.test[,xvars,drop=FALSE], res)
  #     test.stat=evaluate.results(data.test[,yvar,drop=FALSE], pred.data, censor.vec=censor.vec.test,
  #                                trt.vec=trt.vec.test, trtref=trtref, type=type)
  #
  #   }
  #
  #   if (do.cv){
  #     cv.res=cv.aim.rule.batting(y=data.train[,yvar,drop=FALSE],x=data.train[,xvars,drop=FALSE], censor.vec=censor.vec.train,
  #                                trt.vec=trt.vec.train, trtref=trtref, type=type,
  #                                n.boot=n.boot, des.res=des.res, min.sigp.prcnt=min.sigp.prcnt,
  #                                mc.iter=mc.iter, pre.filter=pre.filter, filter.method=filter.method,
  #                                k.fold=k.fold, cv.iter=cv.iter, max.iter=max.iter)
  #
  #   }
  #
  #   if (!is.null(trtvar)){
  #     train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (AIM.Rule, Train)", trt.lab=c("Trt.", "Ctrl."))
  #     if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (AIM.Rule, Test)", trt.lab=c("Trt.", "Ctrl."))
  #   }
  #
  #
  # }

  ############################ output files ########################################
  if (!is.null(out.file)){

    # if (method=="AIM") {
    #   res.temp=res$aim.rule
    #   res.temp=rbind(res.temp,data.frame(variable="Score", direction=">", cutoff=res$bat.cutoff))
    #   res=res.temp
    #   row.names(res)=NULL
    # }
    #
    # if (method=="AIM.Rule") {
    #   res=res$aim.rule
    #   row.names(res)=NULL
    # }
    #



    res.file=paste(file.path,out.file,".",method,".res.csv",sep="")
    train.pval.file=paste(file.path,out.file,".",method,".train.pval.csv",sep="")
    train.ratios.file=paste(file.path,out.file,".",method,".train.ratios.csv",sep="")
    train.gpstats.file=paste(file.path,out.file,".",method,".train.group.stats.csv",sep="")
    test.pval.file=paste(file.path,out.file,".",method,".test.pval.csv",sep="")
    test.ratios.file=paste(file.path,out.file,".",method,".test.ratios.csv",sep="")
    test.gpstats.file=paste(file.path,out.file,".",method,".test.group.stats.csv",sep="")
    cv.pval.file=paste(file.path,out.file,".",method,".cv.pval.csv",sep="")
    cv.ratios.file=paste(file.path,out.file,".",method,".cv.ratios.csv",sep="")
    cv.gpstats.file=paste(file.path,out.file,".",method,".cv.group.stats.csv",sep="")

    write.csv(res,file=res.file, row.names=FALSE)
    write.csv(train.stat$pval,file=train.pval.file, row.names=FALSE)
    if(type!="c") write.csv(train.stat$ratios,file=train.ratios.file, row.names=FALSE)
    write.csv(train.stat$group.stats,file=train.gpstats.file)

    if (!is.null(data.test)){
      write.csv(test.stat$pval,file=test.pval.file, row.names=FALSE)
      if(type!="c") write.csv(test.stat$ratios,file=test.ratios.file, row.names=FALSE)
      write.csv(test.stat$group.stats,file=test.gpstats.file)
    }

    if (do.cv){
      write.csv(cv.res$stats.summary$pvals,file=cv.pval.file, row.names=FALSE)
      if(type!="c") write.csv(cv.res$stats.summary$ratios,file=cv.ratios.file, row.names=FALSE)
      write.csv(cv.res$stats.summary$group.stats,file=cv.gpstats.file)
    }

    if (plots)
    {
      if (!is.null(trtvar))
      {
        train.plot.file=paste(file.path,out.file,".",method,".train.plot.jpg",sep="")
        ggsave(filename=train.plot.file,plot=train.plot, width=4, height=4,dpi=300)
      }

      if (!is.null(data.test))
      {
        if (!is.null(trtvar))
        {
          test.plot.file=paste(file.path,out.file,".",method,".test.plot.jpg",sep="")
          ggsave(filename=test.plot.file,plot=test.plot, width=4, height=4,dpi=300)
        }
      }

      if (do.cv)
      {
        if (!is.null(trtvar))
        {
          cv.plot.file=paste(file.path, out.file,".",method,".cv.plot.jpg",sep="")
          ggsave(filename=cv.plot.file,plot=cv.res$interplot, width=4, height=4,dpi=300)
        }
      }
    }



  }


  ###################################################
  if (is.null(trtvar)) {
    train.plot=NULL
    test.plot=NULL
  }

  if (is.null(data.test)){
    test.stat=NULL
    test.plot=NULL
  }

  if (!do.cv){
    cv.res=NULL
  }

  result=list(res=res, train.stat=train.stat, test.stat=test.stat, cv.res=cv.res,
              train.plot=train.plot, test.plot=test.plot)


  if (!is.null(out.file)){
    save(data.train, data.test, result, file = paste(file.path,out.file,".",method,".RData",sep=""))
  }

  result


}




#' Function for simulated data generation
#' @title data.gen
#' @description Function for simulated data generation
#' @param n Total sample size
#' @param k Number of markers
#' @param prevalence prevalence of predictive biomarkers with values above the cutoff
#' @param prog.eff  effect size \eqn{beta} for prognostic biomarker
#' @param sig2 standard deviation of each marker
#' @param y.sig2 Standard Deviation of the error term in the linear component
#' @param rho rho*sig2 is the entries for covariance matrix between pairs of different k markers
#' @param rhos.bt.real correlation between each prognostic and predictive markers
#' @param a.constent a constant is set such that there is no overall treatment effect
#' @return A list of simulated clinical trial data with heterogeneous prognostic and predictive biomarkers
#' @examples
#' n <- 500
#' k <- 10
#' prevalence <- sqrt(0.5)
#' rho<-0.2
#' sig2 <- 2
#' rhos.bt.real <- c(0, rep(0.1, (k-3)))*sig2
#' y.sig2 <- 1
#' prog.eff <- 0.5
#' effect.size <- 1
#' a.constent <- effect.size/(2*(1-prevalence))
#' ObsData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
#'                     sig2=sig2, y.sig2=y.sig2, rho=rho,
#'                     rhos.bt.real=rhos.bt.real, a.constent=a.constent)
#' @export
data.gen <- function(n, k, prevalence=sqrt(0.5), prog.eff=1, sig2, y.sig2, rho, rhos.bt.real, a.constent){
  covm <- matrix(rho*sig2,k,k)
  diag(covm) <- sig2
  covm[1,3:k] <- rhos.bt.real
  covm[2,3:k] <- rhos.bt.real
  covm[3:k,1] <- rhos.bt.real
  covm[3:k,2] <- rhos.bt.real
  dich.cutoff <- stats::qnorm(prevalence, sd = sqrt(sig2))
  x <- MASS::mvrnorm(n, rep(0,k), covm)
  x.dich <- 1*(x < dich.cutoff)
  w <- x.dich[,3:5]
  trt <- stats::rbinom(n, 1, prob=0.5)
  #trt <- rbinom(n, size = 1, prob = plogis(-1+0.05*w[,1] + 0.25*w[,2] + 0.6*w[,3] + 0.4*w[,1]*w[,3])) # confounding variables

  # predictive effect: x1, x2, prognostic effect x3
  prog.part <- prog.eff*w[,1] + stats::rnorm(n, sd=sqrt(y.sig2))
  pred.part <- a.constent*(-2*prevalence + x.dich[,1] + x.dich[,2]) # simulation based on binary biomarker
  #prog.part <- x[,3] + rnorm(n, sd=sqrt(y.sig2))
  #pred.part <- a.constent*(-2*prevalence + x[,1] + x[,2]) # simulation based on continous biomarker
  y.0 <- pred.part*0 + prog.part
  y.1 <- pred.part*1 + prog.part

  # continuous outcome
  y <- y.1*trt + y.0*(1-trt)

  # create binary outcome
  y.binary.0 <- 1*(stats::plogis(y.0)>prevalence)
  y.binary.1 <- 1*(stats::plogis(y.1)>prevalence)
  y.binary <- y.binary.1*trt + y.binary.0*(1-trt)

  # create time-to-event outcome
  surv.time.0 <- exp(y.0)
  surv.time.1 <- exp(y.1)
  cens.time <- exp(stats::rnorm(n, sd = 3))
  y.time.to.event.0 <- pmin(surv.time.0, cens.time)
  y.time.to.event.1 <- pmin(surv.time.1, cens.time)
  y.time.to.event <- y.time.to.event.1*trt + y.time.to.event.0*(1-trt)
  status.0 <- 1*(surv.time.0 <= cens.time)
  status.1 <- 1*(surv.time.1 <= cens.time)
  status <- status.1*trt + status.0*(1-trt)

  data <- cbind(y, y.binary, y.time.to.event, status,
                y.0, y.1, y.binary.0, y.binary.1,
                y.time.to.event.0, y.time.to.event.1,
                status.0, status.1, trt, x)
  colnames(data) <- c("y", "y.binary", "y.time.to.event", "status",
                      "y.0", "y.1", "y.binary.0", "y.binary.1",
                      "y.time.to.event.0", "y.time.to.event.1",
                      "status.0", "status.1","treatment", sapply(c(1:k), FUN=function(x) paste("x", x, sep="")))
  colnames(x) <- sapply(c(1:k), FUN=function(x) paste("x", x, sep=""))
  colnames(w) <- sapply(c(1:3), FUN=function(x) paste("w", x, sep=""))
  list(data=data, y=y, y.binary = y.binary,
       y.time.to.event = y.time.to.event, status = status,
       y.0=y.0, y.1=y.1, y.binary.0=y.binary.0, y.binary.1=y.binary.1,
       y.time.to.event.0=y.time.to.event.0, y.time.to.event.1=y.time.to.event.1,
       status.0=status.0, status.1=status.1, trt=trt, x=x.dich, x.continues=x,
       w=w, sigpos=(x.dich[,1]==1&x.dich[,2]==1))
}
