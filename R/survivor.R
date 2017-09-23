#' @title Estimating average effect of discrete treatment on survivors
#'
#' @description \code{surv} is used to estimate the mean outcome in the
#' subset of a population which survived having been given
#' levels of a dichotomous (unconfounded) treatment. Taken wholesale from
#' Ed Kennedy.
#'
#' @usage surv(y, a, s, x, nsplits=2, sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet",
#'   "SL.glm.interaction", "SL.mean","SL.ranger","rpart"))
#'
#' @param y outcome of interest.
#' @param a dichotomous treatment.
#' @param t indicator of survival.
#' @param x covariate matrix.
#' @param nsplits integer number of sample splits for nuisance estimation.
#' If nsplits=1, sample splitting is not used, and nuisance functions are estimated
#' on full sample (in which case validity of SEs/CIs requires empirical
#' process conditions). Otherwise must have nsplits>1.
#' @param sl.lib algorithm library for SuperLearner.
#' Default library includes "earth", "gam", "glm", "glmnet", "glm.interaction",
#' "mean", "ranger", "rpart.
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs/CIs/p-values for population means and relevant contrasts.}
#' \item{nuis}{ subject-specific estimates of nuisance functions (i.e., propensity score and outcome regression) }
#' \item{ifvals}{ matrix of estimated influence function values.}
#'
#' @examples
#'
#' @references Robins JM, Rotnitzky A (1995). Semiparametric efficiency in multivariate regression models with missing data. \emph{Journal of the American Statistical Association}.
#' @references Hahn J (1998). On the role of the propensity score in efficient semiparametric estimation of average treatment effects. \emph{Econometrica}.
#' @references van der Laan MJ, Robins JM (2003). \emph{Unified Methods for Censored Longitudinal Data and Causality} (Springer).
#' @references Tsiatis AA (2006). \emph{Semiparametric Theory and Missing Data} (Springer).
#' @references Robins JM, Li L, Tchetgen Tchetgen ET, van der Vaart A (2008). Higher order influence functions and minimax estimation of nonlinear functionals. \emph{Probability and Statistics: Essays in Honor of David A. Freedman}.
#' @references Zheng W, van der Laan (2010). Asymptotic theory for cross-validated targeted maximum likelihood estimation \emph{UC Berkeley Division of Biostatistics Working Paper Series}.
#' @references Chernozhukov V, Chetverikov V, Demirer M, et al (2016). Double machine learning for treatment and causal parameters.
#'
surv <- function(y,a,t,x, nsplits=2,
                sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.rpart")){

  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require("rpart")

  n <- dim(x)[1]
  avals <- names(table(a))
  tvals <- names(table(t))
  n.avals <- length(avals)
  if(n.avals!=2) {print('for now, can only handle dichotomous treatment'); stop()}
  n.tvals <- length(tvals)
  pb <- txtProgressBar(min=0, max=2*nsplits*n.avals, style=3)

  s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])

  muhat <- as.data.frame(matrix(NA,nrow=n,ncol=n.avals))
  colnames(muhat) <- paste("a",avals,sep="")
  pihat.a <- pihat.t <- muhat

  pbcount <- 0
  for (i in 1:n.avals){
    if (i==1){ Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }
    for (vfold in 1:nsplits){

      train <- s!=vfold; test <- s==vfold

      ttrain <- as.data.frame(cbind(a[train],x[train,]))
      ttest <- as.data.frame(cbind(a[test],x[test,]))
      names(ttrain) <- names(ttest) <- names(as.data.frame(cbind(a,x)))

      if (nsplits==1){ train <- test }

      # estimate propensity score of a and t for a = avals[1]
      if (i != n.avals){
        pifit.a <- SuperLearner(as.numeric(a==avals[i])[train],as.data.frame(x[train,]),
                              newX=as.data.frame(x[test,]), SL.library=sl.lib, family=binomial)
        pihat.a[test,i] <-pifit.a$SL.predict

        pifit.t <- SuperLearner(as.numeric(t==1)[train & a==avals[i]],ttest,newX=ttest, SL.library=sl.lib, family=binomial)
        pihat.t[test & a==avals[i],i] <-pifit.t$SL.predict
        }


      # estimate regression function
      xtrain <- as.data.frame(x[a==avals[i] & t == 1 & train,])
      xtest <- as.data.frame(x[test & t == 1,])
      names(xtrain)<- names(xtest)
      mufit <- SuperLearner(y[a==avals[i] & t == 1 & train],xtrain,newX=xtest, SL.library=sl.lib)
      muhat[test & t == 1,i] <- mufit$SL.predict

      Sys.sleep(0.1)
      setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
    }
    if (i == n.avals){
      pihat.a[,i] <- 1 - apply(pihat.a, 1, sum, na.rm = T)
      pifit.t <- SuperLearner(as.numeric(t==1)[train & a==avals[i]],ttest,newX=ttest, SL.library=sl.lib, family=binomial)
      pihat.t[test & a==avals[i],i] <-pifit.t$SL.predict
      }
  }

  amat <- matrix(rep(a,n.avals),nrow=n,byrow=F)
  alevel <- matrix(rep(1:n.avals, rep(n,n.avals)),nrow=n,byrow=F)
  tmat <- matrix(rep(t,n.tvals),nrow=n,byrow=F)
  ymat <- matrix(rep(y,n.avals),nrow=n,byrow=F)

  ifvals <- as.matrix( (amat==alevel)*(ymat-muhat)/(pihat.a*pihat.t) + muhat )
  ifvals[t==0,] <- 0

  est <- apply(ifvals,2,mean)
  se <- apply(ifvals,2,sd)/sqrt(n)
  ci.ll <- est-1.96*se; ci.ul <- est+1.96*se
  pval <- round(2*(1-pnorm(abs(est/se))),3)
  paste("E{Y(",avals,")}")
  res1 <- data.frame(parameter=paste("E{Y(",avals,")}",sep=""), est,se,ci.ll,ci.ul,pval)

  signdist <- function(x){ c(as.dist(outer(x,x,'-'))) }
  ifvals2 <- t(apply(ifvals,1,signdist))
  if (n.avals==2){ ifvals2 <- t(ifvals2)  }

  tmp <- expand.grid(1:n.avals,1:n.avals)
  tmp2 <- tmp[tmp[,1]>tmp[,2],]
  contrasts <- apply(cbind(avals[tmp2[,1]],avals[tmp2[,2]]),1,paste,collapse=")-Y(")
  contrasts <- paste("E{Y(",contrasts,")}",sep="")

  est2 <- apply(ifvals2,2,mean)
  se2 <- apply(ifvals2,2,sd)/sqrt(n)
  ci.ll2 <- est2-1.96*se2; ci.ul2 <- est2+1.96*se2
  pval2 <- round(2*(1-pnorm(abs(est2/se2))),3)
  res2 <- data.frame(parameter=contrasts,est=est2,se=se2,ci.ll=ci.ll2,ci.ul=ci.ul2,pval=pval2)

  res <- rbind(res1,res2); rownames(res) <- NULL

  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount)
  close(pb)

  nuis <- as.data.frame(cbind(pihat.a, pihat.t,muhat))
  colnames(nuis) <- paste(rep(c("piA", "piT","mu"), rep(n.avals,3)), colnames(nuis), sep="_")

  print(res)
  return(invisible(list(res=res, nuis=nuis, ifvals=as.data.frame(ifvals) )))

}


surv2 <- function(y,a,t,x, nsplits=2,
                  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.rpart")){

  # gonna try to just make this real simple
  # first column of everything is a = 1, second column is a = 0
  if(nsplits != 2){cat('Can only do 2 splits'); stop()}
  if(length(unique(a)) != 2){cat('Can only do dichotomous A'); stop()}

  pb <- txtProgressBar(min=0, max=2*nsplits*2, style=3)

  n = dim(x)[1]
  s0 = sample(1:n, n/2, replace = FALSE)
  s1 = c(1:n)[-s0]
  samples = cbind(s0,s1)

  piA.mat <- piT.mat <- muhat <- matrix(rep(NA,2*n),ncol = 2)

  for(i in 1:2){
    if (i==1){ Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }
    s = samples[,i]

    # a = 1 ----
    # P(A=1|X)
    pifit.a <- SuperLearner(as.numeric(a==1)[s],as.data.frame(x[s,]), newX=as.data.frame(x[-s,]), SL.library=sl.lib, family=binomial)
    piA.mat[-s,1] <-pifit.a$SL.predict

    Sys.sleep(0.1)
    setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

    # P(T=1|A=1,X)
    pifit.t <- SuperLearner(as.numeric(t==1)[s & (a==1)],as.data.frame(x[s & (a==1),]), newX=as.data.frame(x[s,]), SL.library=sl.lib, family=binomial)
    piT.mat[-s,1] <-pifit.t$SL.predict

    Sys.sleep(0.1)
    setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

    # E(Y|A=a1,T=1,X)
    xtrain <- as.data.frame(x[a==1 & t==1 & s,])
    xtest <- as.data.frame(x[-s,])
    names(xtrain)<- names(xtest)
    mufit <- SuperLearner(y[a==1 & t==1 & s],xtrain,newX=xtest, SL.library=sl.lib)
    muhat[-s,1] <- mufit$SL.predict

    Sys.sleep(0.1)
    setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

    # a = 0 ----
    # P(A=0|X)
    piA.mat[-s,2] <- 1 - piA.mat[-s,1]

    # P(T=1|A=0,X)
    pifit.t <- SuperLearner(as.numeric(t==1)[s & (a==0)],as.data.frame(x[s & (a==0),]), newX=as.data.frame(x[-s,]), SL.library=sl.lib, family=binomial)
    piT.mat[-s,2] <-pifit.t$SL.predict

    # E(Y|A=0,T=1,X)
    xtrain <- as.data.frame(x[a==0 & t==1 & s,])
    xtest <- as.data.frame(x[-s,])
    names(xtrain)<- names(xtest)
    mufit <- SuperLearner(y[a==0 & t==1 & s],xtrain,newX=xtest, SL.library=sl.lib)
    muhat[-s,2] <- mufit$SL.predict
  }


  # get IF ----
  alevel <- matrix(c(rep(1,n),rep(0,n)),nrow=n,byrow=F)
  amat <- matrix(rep(a,2),nrow=n,byrow=F)
  ymat <- matrix(rep(y,2),nrow=n,byrow=F)

  ifvals <- as.matrix( (amat==alevel & t==1)*(ymat-muhat)/(piA.mat*piT.mat) + muhat )

  # get estimates ----
  est <- apply(ifvals,2,mean, na.rm = T)
  se <- apply(ifvals,2,sd,na.rm = T)/sqrt(n)
  ci.ll <- est-1.96*se; ci.ul <- est+1.96*se
  pval <- round(2*(1-pnorm(abs(est/se))),3)
  paste("E{Y(",c(1,0),")}")
  res1 <- data.frame(parameter=paste("E{Y(",c(1,0),")}",sep=""), est,se,ci.ll,ci.ul,pval)

  ifvals2 <- ifvals[,1] - ifvals[,2]
  contrasts <- apply(cbind(1,0),1,paste,collapse=")-Y(")
  contrasts <- paste("E{Y(",contrasts,")}",sep="")

  est2 <- mean(ifvals2, na.rm = T)
  se2 <- sd(ifvals2, na.rm = T)/sqrt(n)
  ci.ll2 <- est2-1.96*se2; ci.ul2 <- est2+1.96*se2
  pval2 <- round(2*(1-pnorm(abs(est2/se2))),3)
  res2 <- data.frame(parameter=contrasts,est=est2,se=se2,ci.ll=ci.ll2,ci.ul=ci.ul2,pval=pval2)

  res <- rbind(res1,res2); rownames(res) <- NULL

  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount)
  close(pb)

  nuis <- as.data.frame(cbind(piA.mat, piT.mat,muhat))
  colnames(nuis) <- paste(rep(c("piA", "piT","mu"), rep(2,3)), colnames(nuis), sep="_")

  print(res)
  return(invisible(list(res=res, nuis=nuis, ifvals=as.data.frame(ifvals) )))


}

surv.saved <- function(y,a,t,x, nsplits=2,
                  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.rpart")){

  if(nsplits != 2){cat('Can only do 2 splits'); stop()}
  if(length(unique(a)) != 2){cat('Can only do dichotomous A'); stop()}

  n = dim(x)[1]
  s0 = sample(1:n, n/2, replace = FALSE)
  s1 = c(1:n)[-s0]
  samples = cbind(s0,s1)

  avals <- names(table(a)); n.avals = 2

  piA.mat <- piT.mat <- muhat <- matrix(rep(NA,2*n),ncol = 2)

  for(i in 1:2){
    s = samples[,i]

    # a = a1 ----
    # P(A=a1|X)
    pifit.a <- SuperLearner(as.numeric(a==avals[1])[s],as.data.frame(x[s,]), newX=as.data.frame(x[-s,]), SL.library=sl.lib, family=binomial)
    piA.mat[-s,1] <-pifit.a$SL.predict

    # P(T=1|A=a1,X)
    pifit.t <- SuperLearner(as.numeric(t==1)[s & (a==avals[1])],as.data.frame(x[s & (a==avals[1]),]), newX=as.data.frame(x[s,]), SL.library=sl.lib, family=binomial)
    piT.mat[-s,1] <-pifit.t$SL.predict

    # E(Y|A=a1,T=1,X)
    xtrain <- as.data.frame(x[a==avals[1] & t==1 & s,])
    xtest <- as.data.frame(x[-s,])
    names(xtrain)<- names(xtest)
    mufit <- SuperLearner(y[a==avals[1] & t==1 & s],xtrain,newX=xtest, SL.library=sl.lib)
    muhat[-s,1] <- mufit$SL.predict

    # a = a2 ----
    # P(A=a2|X)
    piA.mat[-s,2] <- 1 - piA.mat[-s,1]

    # P(T=1|A=a2,X)
    pifit.t <- SuperLearner(as.numeric(t==1)[s & (a==avals[2])],as.data.frame(x[s & (a==avals[2]),]), newX=as.data.frame(x[-s,]), SL.library=sl.lib, family=binomial)
    piT.mat[-s,2] <-pifit.t$SL.predict

    # E(Y|A=0,T=1,X)
    xtrain <- as.data.frame(x[a==avals[2] & t==1 & s,])
    xtest <- as.data.frame(x[-s,])
    names(xtrain)<- names(xtest)
    mufit <- SuperLearner(y[a==avals[2] & t==1 & s],xtrain,newX=xtest, SL.library=sl.lib)
    muhat[-s,2] <- mufit$SL.predict
  }


  # get IF ----
  alevel <- matrix(c(rep(avals[1],n),rep(avals[2],n)),nrow=n,byrow=F)
  amat <- matrix(rep(a,n.avals),nrow=n,byrow=F)
  ymat <- matrix(rep(y,n.avals),nrow=n,byrow=F)

  ifvals <- as.matrix( (amat==alevel & t==1)*(ymat-muhat)/(piA.mat*piT.mat) + muhat )

  # get estimates ----
  est <- apply(ifvals,2,mean, na.rm = T)
  se <- apply(ifvals,2,sd,na.rm = T)/sqrt(n)
  ci.ll <- est-1.96*se; ci.ul <- est+1.96*se
  pval <- round(2*(1-pnorm(abs(est/se))),3)
  paste("E{Y(",avals,")}")
  res1 <- data.frame(parameter=paste("E{Y(",avals,")}",sep=""), est,se,ci.ll,ci.ul,pval)

  # signdist <- function(x){ c(as.dist(outer(x,x,'-'))) }
  # ifvals2 <- t(apply(ifvals,1,signdist))
  # ifvals2 <- t(ifvals2)


  tmp <- expand.grid(1:n.avals,1:n.avals)
  tmp2 <- tmp[tmp[,1]>tmp[,2],]
  ifvals2 <- ifvals[,tmp2[,1]] - ifvals[,tmp2[,2]]
  contrasts <- apply(cbind(avals[tmp2[,1]],avals[tmp2[,2]]),1,paste,collapse=")-Y(")
  contrasts <- paste("E{Y(",contrasts,")}",sep="")

  est2 <- mean(ifvals2, na.rm = T)
  se2 <- sd(ifvals2, na.rm = T)/sqrt(n)
  ci.ll2 <- est2-1.96*se2; ci.ul2 <- est2+1.96*se2
  pval2 <- round(2*(1-pnorm(abs(est2/se2))),3)
  res2 <- data.frame(parameter=contrasts,est=est2,se=se2,ci.ll=ci.ll2,ci.ul=ci.ul2,pval=pval2)

  res <- rbind(res1,res2); rownames(res) <- NULL

  nuis <- as.data.frame(cbind(piA.mat, piT.mat,muhat))
  colnames(nuis) <- paste(rep(c("piA", "piT","mu"), rep(n.avals,3)), colnames(nuis), sep="_")

  print(res)
  return(invisible(list(res=res, nuis=nuis, ifvals=as.data.frame(ifvals) )))


}

