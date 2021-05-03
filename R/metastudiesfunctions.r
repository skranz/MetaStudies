
#' Perform Andrews and Kasy (2019) estimation
#'
#' @param X vector of reported coefficients
#' @param sigma vector of reported stanard errors
#' @param cutoffs significance thresholds that define intervals of different publication probabilities, e.g. \code{c(1.645, 1.96, 2.576)}.
#' @param symmetric if TRUE assume that positive and negative z-statistics are symmetrically distributed.
#' @param model either "normal" or "t". The assumed functional form of the distribution absent publication bias.
#' @param eval.max,iter.max,abs.tol Control parameters for [stats::nlminb].
metastudies_estimation=function(X,sigma, cutoffs=c(1.96),symmetric=FALSE, model="normal", eval.max=10^5,iter.max=10^5,abs.tol=10^(-8),stepsize=10^(-6)){

  restore.point("metastudies_estimation")
  nn=length(X)
  #%regressors for step function p
  TT=X/sigma;
  Tpowers=Tpowers_fun(TT,cutoffs, symmetric)

  if (model=="normal"){
    LLH <- function (Psi) VariationVarianceLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)],  1),cutoffs,symmetric, X, sigma, Tpowers)
    psi_hat0 = c(0,1, rep(1,length(cutoffs))) #starting values
  } else if (model=="t"){
    LLH <- function (Psi) VariationVarianceLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2,3)],  1),cutoffs,symmetric, X, sigma, Tpowers, df=Psi[3])
    psi_hat0 = c(0,1,10, rep(1,length(cutoffs)))    #starting values
  }


  LLH_only<-function (Psi){
      A<-LLH(Psi);
      return(A$LLH)
  }

  lower.b = c(-Inf,rep(0,length(psi_hat0)-1))
  upper.b=rep(Inf,length(psi_hat0))

  findmin<-nlminb(objective=LLH_only, start=psi_hat0,lower=lower.b,upper=upper.b,control = list(eval.max = eval.max, iter.max = iter.max, abs.tol = abs.tol));
  psi_hat<-findmin$par
  LLHmax<-findmin$objective

  Var_robust<-RobustVariance(stepsize, nn, psi_hat, LLH,1:nn);
  se_robust<-sqrt(diag(Var_robust));

  k = length(cutoffs)
  prob.df = data.frame(
    z.min = c(ifelse(symmetric,0,-Inf), cutoffs),
    z.max = c(cutoffs,Inf),
    pub.prob = c(psi_hat[(length(psi_hat)-k+1):(length(psi_hat))],1)
  )

  dat = data.frame(
    X = X,
    sigma = sigma,
    z = X / sigma
  )

  if (symmetric) {
    dat$interval.ind = findInterval(abs(dat$z),c(0,cutoffs, Inf))
    dat$pub.prob = prob.df$pub.prob[dat$interval.ind]
  } else {
    dat$interval.ind = findInterval(dat$z,c(-Inf,cutoffs, Inf))
    dat$pub.prob = prob.df$pub.prob[dat$interval.ind]
  }

  res = list(psi_hat=psi_hat, SE=se_robust, X=X, sigma=sigma, cutoffs=cutoffs, model=model, symmetric=symmetric, psi_vcov=Var_robust, prob.df = prob.df, est_tab = NA, dat=dat)
  res$est_tab = estimatestable(res)
  class(res) = c("MetaStudy","list")
  res
}

print.MetaStudy = function(ms) {
  cat("\nEstimated Metastudy\n")
  str(ms)
}



Tpowers_fun=function(TT,cutoffs,symmetric){
  n=length(TT)
  Tpowers=matrix (0,n,length(cutoffs)+1)
  if (symmetric) TT=abs(TT)
  Tpowers[,1]=TT<cutoffs[1]
  if (length(cutoffs)>1) {
    for (m in 2:length(cutoffs)) {
      Tpowers[,m]=(TT<cutoffs[m])*(TT>=cutoffs[m-1]);
    }
  }
  Tpowers[,length(cutoffs)+1]=(TT)>=cutoffs[length(cutoffs)]
  Tpowers
}



VariationVarianceLogLikelihood <-function(lambdabar, tauhat, betap,
                                          cutoffs, symmetric, X, sigma, Tpowers,
                                          df=Inf) { #if df argument is provided, switch to t-dist
  n=length(X);
  betap=as.matrix(betap, length(betap),1);

  #  %vector of estimated publication probabilities
  phat=Tpowers%*%betap;

  #  %%%%%%%%%%%%%%%%%%%%%%%%
  #@  vector of un-truncated likelihoods
  sigmatilde=sqrt(sigma^2 + tauhat^2)
  fX=dt((X-lambdabar)/sigmatilde, df)/sigmatilde


  #  normalizingconstant
  normalizedcutoffs=(sigma/sigmatilde)%*%t(cutoffs) - (lambdabar/sigmatilde)
  if (symmetric){
    normalizednegativecutoffs=(sigma/sigmatilde)%*%t(-cutoffs) - (lambdabar/sigmatilde)
    cdfs= pt(normalizedcutoffs, df)-pt(normalizednegativecutoffs, df)
  } else {
    cdfs= pt(normalizedcutoffs, df)
  }
  cdfs=cbind(rep(0,n),cdfs,rep(1,n))
  cellprobas =cdfs[,-1] - cdfs[,-(length(cutoffs)+2)]
  normalizingconst=cellprobas%*%betap;

  #likelihood for each observation
  L<-phat*fX/normalizingconst;
  logL<-log(L);
  # objective function; note the sign flip, since we are doing minimization
  LLH<--sum(logL);

  if (is.nan(LLH)){
      show(lambdabar)
      show(tauhat)
      show(betap)
  }

  return(list(LLH=LLH,logL=logL))

}


estimatestable=function(x = NULL, psi_hat=x$psi_hat, SE=x$SE, cutoffs=x$cutoffs, symmetric=x$symmetric, model=x$model) {
  l=length(psi_hat)
  estimates=matrix(0,2,l)
  estimates[1,]=psi_hat
  estimates[2,]=SE
  rownames(estimates)=c("estimate", "standard error")
  colnames(estimates)=rep(" ",l)
  colnames(estimates)[1]=intToUtf8(956) #mu
  colnames(estimates)[2]=intToUtf8(964) #tau
  if (model=="t"){
      colnames(estimates)[3]="df"
      shift=1
  } else {shift = 0}

  if (symmetric){
      colnames(estimates)[3+shift]=paste("[0,", cutoffs[1],"]")
      for (i in seq(2, length(cutoffs), length=max(0, length(cutoffs) - 1))) {
          colnames(estimates)[2+i+shift]=paste("(", cutoffs[i-1], ",", cutoffs[i],"]")
      }
  } else {
      colnames(estimates)[3+shift]=paste("(-", intToUtf8(8734), ",", cutoffs[1],"]")
      for (i in seq(2, length(cutoffs), length=max(0, length(cutoffs) - 1))) {
          colnames(estimates)[2+i+shift]=paste("(", cutoffs[i-1], ",", cutoffs[i],"]")
      }
  }
  estimates
}

