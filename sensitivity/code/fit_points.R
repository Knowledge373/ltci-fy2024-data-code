# Fixed author fitting only; raw X/Y stay in memory. One formal GLM call per eligible cell.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==7L)
out <- args[1]; mode <- args[2]; expectedX <- args[3]; expectedY <- args[4]; expectedC <- args[5]
do_ols <- args[6]=='TRUE'; do_frac <- args[7]=='TRUE'
stopifnot(as.character(getRversion())=='4.6.1',as.character(packageVersion('stats'))=='4.6.1')
exactcsv <- function(z,path) {
  for(k in names(z)) if(is.numeric(z[[k]])) z[[k]] <- ifelse(is.na(z[[k]]),NA_character_,sprintf('%.17g',z[[k]]))
  write.table(z,path,sep=',',row.names=FALSE,na='',quote=TRUE)
}
sha <- function(x,name,size=8L) {
 p <- file.path(out,paste0(name,'.sha256')); h <- pipe(paste('sha256sum >',shQuote(p)),open='wb')
 writeBin(x,h,size=size,endian='little');close(h)
 strsplit(readLines(p,warn=FALSE)[1],' ')[[1]][1]
}
rangeinfo <- function(z) c(min=if(all(is.finite(z)))min(z) else NA_real_,max=if(all(is.finite(z)))max(z) else NA_real_,nonfinite=sum(!is.finite(z)),below0=sum(z<0,na.rm=TRUE),above1=sum(z>1,na.rm=TRUE),numeric_boundary=sum(z<=10*.Machine$double.eps | z>=1-10*.Machine$double.eps,na.rm=TRUE))
rk <- function(x) {s <- svd(x,nu=0,nv=0)$d; tol <- if(length(s))max(dim(x))*.Machine$double.eps*max(s) else 0; list(rank=sum(s>tol),tolerance=tol,s=s)}
d <- read.csv(file('stdin'),check.names=FALSE)
X <- cbind(intercept=1,as.matrix(d[setdiff(names(d),c('y','pref'))]));Y<-d$y; p<-ncol(X);n<-nrow(X)
stopifnot(sha(as.double(t(X)),'R_X_LE_F64_C')==expectedX,sha(as.double(Y),'R_Y_LE_F64')==expectedY,sha(as.integer(d$pref),'R_CLUSTER_LE_I32',4L)==expectedC)
stopifnot(all(is.finite(X)),all(is.finite(Y)),all(Y>=0 & Y<=1),all(d$D %in% c(0,1)),colnames(X)[2]=='D')
X0<-X;X1<-X;X0[,2]<-0;X1[,2]<-1
if(do_ols) {
  warns<-character();err<-''
  fit <- tryCatch(withCallingHandlers(stats::lm.fit(X,Y,tol=1e-7,singular.ok=FALSE),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL})
  writeLines(warns,file.path(out,'OLS_WARNINGS.txt'));writeLines(err,file.path(out,'OLS_ERROR.txt'))
  if(!is.null(fit)) {
   b<-fit$coefficients;pred<-drop(X%*%b);pred0<-drop(X0%*%b);pred1<-drop(X1%*%b)
   exactcsv(data.frame(term=colnames(X),estimate_fraction=unname(b),estimate_pp=100*unname(b)),file.path(out,'OLS_COEFFICIENTS.csv'))
   exactcsv(data.frame(n=n,G=length(unique(d$pref)),p=p,rank=fit$rank,beta_fraction=unname(b[2]),beta_pp=100*unname(b[2]),max_abs_Xt_residual=max(abs(crossprod(X,Y-pred))),normal_equation_normalized=max(abs(crossprod(X,Y-pred))/(1+colSums(abs(X)))),SSE=sum((Y-pred)^2),residual_df=n-fit$rank,coefficients_finite=all(is.finite(b)),predictions_finite=all(is.finite(c(pred,pred0,pred1)))),file.path(out,'OLS_FIT_DIAGNOSTICS.csv'))
   exactcsv(data.frame(set=c('observed','D0','D1'),rbind(rangeinfo(pred),rangeinfo(pred0),rangeinfo(pred1))),file.path(out,'OLS_PREDICTION_DIAGNOSTICS.csv'))
   loo<-list()
   for(g in sort(unique(d$pref))) {
    keep<-d$pref!=g;cl<-d$pref[keep];pr<-sort(unique(cl));xx<-cbind(X[keep,1:6,drop=FALSE],sapply(pr[-1],function(k)as.numeric(cl==k)))
    rx<-rk(xx)$rank;rz<-rk(xx[,-2,drop=FALSE])$rank;identified<-rx-rz==1;beta<-NA_real_;fitrank<-NA_integer_;loerr<-'';calls<-0L
    if(rx==ncol(xx) && identified) {calls<-1L;lf<-tryCatch(stats::lm.fit(xx,Y[keep],tol=1e-7,singular.ok=FALSE),error=function(e){loerr<<-conditionMessage(e);NULL});if(!is.null(lf)){fitrank<-lf$rank;if(fitrank==ncol(xx)&&all(is.finite(lf$coefficients)))beta<-lf$coefficients[2]}}
    loo[[length(loo)+1L]]<-data.frame(removed_prefecture_code=g,n=sum(keep),rank_X=rx,rank_Z=rz,p=ncol(xx),R_fit_rank=fitrank,D_identified=identified,D_estimable=is.finite(beta),beta_fraction=beta,beta_pp=100*beta,beta_change_pp=100*(beta-b[2]),auxiliary_OLS_fit_calls=calls,error=loerr)
   }
   exactcsv(do.call(rbind,loo),file.path(out,'LEAVE_ONE_PREFECTURE_OUT.csv'))
  }
}
if(do_frac) {
  warns<-character();err<-''
  fit<-tryCatch(withCallingHandlers(stats::glm.fit(x=X,y=Y,weights=rep(1,n),start=rep(0,p),etastart=NULL,mustart=NULL,offset=rep(0,n),family=stats::quasibinomial('logit'),control=stats::glm.control(epsilon=1e-10,maxit=100,trace=FALSE),intercept=TRUE,singular.ok=FALSE),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL})
  writeLines(warns,file.path(out,'FRACTIONAL_WARNINGS.txt'));writeLines(err,file.path(out,'FRACTIONAL_ERROR.txt'))
  if(!is.null(fit)) {
   b<-fit$coefficients;f<-stats::quasibinomial('logit');eta<-drop(X%*%b);eta0<-drop(X0%*%b);eta1<-drop(X1%*%b);mu<-f$linkinv(eta);mu0<-f$linkinv(eta0);mu1<-f$linkinv(eta1)
   objective<-sum(vapply(seq_along(Y),function(i){(if(Y[i]==0)0 else Y[i]*log(mu[i]))+(if(Y[i]==1)0 else (1-Y[i])*log1p(-mu[i]))},numeric(1)))
   score<-drop(crossprod(X,Y-mu));denom<-1+colSums(abs(X)); normscore<-max(abs(score)/denom);fw<-mu*(1-mu)
   exactcsv(data.frame(term=colnames(X),estimate=unname(b)),file.path(out,'FRACTIONAL_COEFFICIENTS.csv'));bhash<-sha(as.double(b),'FRACTIONAL_COEFFICIENTS_LE_F64')
   exactcsv(data.frame(row_index=seq_len(n),eta_observed=eta,eta_D0=eta0,eta_D1=eta1,p_observed=mu,p_D0=mu0,p_D1=mu1,working_weight_last_IWLS=fit$weights,final_weight=fw),file.path(out,'FRACTIONAL_DIAGNOSTIC_ARRAYS.csv'))
   exactcsv(data.frame(term=colnames(X),raw_score=score,denominator=denom,normalized_abs_score=abs(score)/denom),file.path(out,'FRACTIONAL_SCORE.csv'))
   exactcsv(data.frame(set=c('observed','D0','D1'),rbind(rangeinfo(mu),rangeinfo(mu0),rangeinfo(mu1))),file.path(out,'FRACTIONAL_PREDICTION_DIAGNOSTICS.csv'))
   exactcsv(data.frame(formal_fit_calls=1,n=n,p=p,rank=fit$rank,iter=fit$iter,converged=fit$converged,boundary=fit$boundary,deviance=fit$deviance,objective=objective,normalized_score=normscore,delta_pp=100*mean(mu1-mu0),coefficients_finite=all(is.finite(b)),eta_finite=all(is.finite(c(eta,eta0,eta1))),predictions_finite=all(is.finite(c(mu,mu0,mu1))),working_weights_positive_finite=all(is.finite(fit$weights)&fit$weights>0),final_weights_positive_finite=all(is.finite(fw)&fw>0),coefficient_sha256=bhash),file.path(out,'FRACTIONAL_FIT_DIAGNOSTICS.csv'))
   # Persist only derived diagnostic arrays, coefficients and scalar results, never fit object or X/Y.
   saveRDS(list(coefficients=b,eta=list(observed=eta,D0=eta0,D1=eta1),probability=list(observed=mu,D0=mu0,D1=mu1),working_weights=fit$weights,final_weights=fw,raw_score=score,score_denominator=denom,normalized_score=normscore,objective=objective,deviance=fit$deviance,delta_pp=100*mean(mu1-mu0),rank=fit$rank,iter=fit$iter,converged=fit$converged,boundary=fit$boundary),file.path(out,'FRACTIONAL_DERIVED_DIAGNOSTICS.rds'),version=3)
  }
}
cat('FORMAL_POINT_STAGE_COMPLETE\n')
