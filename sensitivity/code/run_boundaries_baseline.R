# Official production boundary runner. Raw X/Y and full bootstrap arrays stay in memory.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
out <- args[1]; contrast <- args[2]; max_new <- as.integer(args[3])
script_file <- local({ z <- commandArgs(FALSE); k <- match('-f', z); if (!is.na(k)) z[k+1L] else sub('^--file=', '', grep('^--file=', z, value=TRUE)[1]) })
module_root <- normalizePath(file.path(dirname(script_file), '..'))
root <- file.path(module_root, 'environment')
.libPaths(c(file.path(root,'r-library'),.libPaths()))
library(fwildclusterboot)
stopifnot(packageVersion('fwildclusterboot')=='0.14.3', as.character(getRversion())=='4.6.1')
source(file.path(module_root, 'code', 'boundary_diagnostics.R'))
# Check installed function bodies against pinned upstream source.
e <- new.env(parent=asNamespace('fwildclusterboot'))
for (f in c('arg_checks.R','get_weights.R','p_val_null.R')) sys.source(file.path(root,'verification-source',f),envir=e)
for (n in c('check_set_full_enumeration','get_weights','p_val_null2')) {
  stopifnot(identical(body(get(n,e)),body(getFromNamespace(n,'fwildclusterboot'))),
            identical(formals(get(n,e)),formals(getFromNamespace(n,'fwildclusterboot'))))
}
write_exact_csv <- function(x,path,append=FALSE) {
  for(n in names(x)) if(is.numeric(x[[n]])) x[[n]] <- ifelse(is.na(x[[n]]),NA_character_,sprintf('%.17g',x[[n]]))
  write.table(x,path,sep=',',row.names=FALSE,col.names=!append,append=append,na='',quote=TRUE)
}
dir.create(file.path(out,'boundary_rows'),showWarnings=FALSE)
d <- read.csv(file('stdin'),check.names=FALSE)
stopifnot(all(is.finite(as.matrix(d))), all(d$D %in% c(0,1)))
G <- length(unique(d$pref)); expected <- expected_rademacher_draws(G,99999L)
nuisance <- setdiff(names(d),c('y','pref'))
f <- reformulate(nuisance,response='y')
fit <- lm(f,data=d,na.action=na.fail)
X <- model.matrix(fit)
stopifnot(nobs(fit)==nrow(d),fit$rank==ncol(X),!anyNA(coef(fit)), all(is.finite(coef(fit))))
write_exact_csv(data.frame(term=names(coef(fit)),estimate_fraction=unname(coef(fit)),estimate_pp=unname(coef(fit))*100),file.path(out,'OLS_COEFFICIENTS.csv'))
# Authorized model-dependence and leave-one-prefecture-out influence diagnostics.
pred <- fitted(fit); cf0 <- pred-coef(fit)['D']*d$D; cf1 <- cf0+coef(fit)['D']
fd <- data.frame(n=nrow(d),G=G,p=ncol(X),rank=fit$rank,beta_fraction=unname(coef(fit)['D']),beta_pp=unname(coef(fit)['D'])*100,
 max_abs_Xt_residual=max(abs(crossprod(X,resid(fit)))),fitted_min=min(pred),fitted_max=max(pred),
 fitted_below0=sum(pred<0),fitted_above1=sum(pred>1),counterfactual0_below0=sum(cf0<0),counterfactual0_above1=sum(cf0>1),
 counterfactual1_below0=sum(cf1<0),counterfactual1_above1=sum(cf1>1),CI_computed=FALSE)
write_exact_csv(fd,file.path(out,'FIT_DIAGNOSTICS.csv'))
loo <- do.call(rbind,lapply(sort(unique(d$pref)),function(g) {
 z <- d[d$pref!=g,]; zcols <- nuisance[!vapply(z[nuisance],function(v)length(unique(v))==1L,logical(1))]
 lf <- lm(reformulate(zcols,response='y'),data=z,na.action=na.fail)
 beta <- unname(coef(lf)['D']); estimable <- 'D' %in% names(coef(lf)) && is.finite(beta)
 data.frame(removed_prefecture_code=g,n=nrow(z),rank=lf$rank,columns=length(coef(lf)),
 D_estimable=estimable,beta_fraction=if(estimable)beta else NA_real_,beta_pp=if(estimable)100*beta else NA_real_,
 beta_change_pp=if(estimable)100*(beta-coef(fit)['D']) else NA_real_)
}))
write_exact_csv(loo,file.path(out,'LEAVE_ONE_PREFECTURE_OUT.csv'))
writeLines(c(paste('R',R.version.string),paste('fwildclusterboot',packageVersion('fwildclusterboot')),
 paste('RNGkind',paste(RNGkind(),collapse=';')),paste('dqrng',packageVersion('dqrng')),
 'official_source_commit=336bb574eba169ac0183317f01d0564791d8122f',
 'seeds: set.seed(20261015); dqrng::dqset.seed(20261015), reset per boundary',
 'nthreads=1; engine=R; bootstrap_type=fnw11; type=rademacher; sampling=dqrng; impose_null=TRUE',
 'ssc: adj=TRUE; fixef.K=none; cluster.adj=TRUE; cluster.df=conventional',
 'sign_draw_reproducible_reference: official get_weights pinned source, G/expected_draws, both seeds reset, numeric pref ascending',
 'full bootstrap arrays and raw X/Y are not persisted'),file.path(out,'RUNTIME.txt'))
path <- file.path(out,'BOUNDARIES.csv')
old <- if(file.exists(path)) read.csv(path,stringsAsFactors=FALSE) else data.frame()
if(nrow(old)) stopifnot(!anyDuplicated(old$id),all(old$G==G),all(old$requested_B==99999L))
new_count <- 0L
for(j in 1:200) for(side in c('lower','upper')) {
 id <- sprintf('%s_j%03d_%s',contrast,j,side)
 if(nrow(old) && id %in% old$id) next
 if(new_count>=max_new || file.exists(file.path(out,'STOP_AFTER_BOUNDARY'))) break
 direction <- if(side=='lower') '>' else '<'; delta_pp <- j/100; r <- if(side=='lower') -j/10000 else j/10000
 set.seed(20261015L); dqrng::dqset.seed(20261015L)
 warns <- character(); err <- ''
 start <- Sys.time()
 b <- tryCatch(withCallingHandlers(boottest(fit,param='D',clustid='pref',bootcluster='pref',B=99999L,r=r,
    impose_null=TRUE,engine='R',bootstrap_type='fnw11',type='rademacher',sampling='dqrng',nthreads=1L,
    p_val_type=direction,sign_level=1/60,conf_int=FALSE,
    ssc=boot_ssc(adj=TRUE,fixef.K='none',cluster.adj=TRUE,cluster.df='conventional')),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e){err<<-conditionMessage(e);NULL})
 elapsed <- as.numeric(difftime(Sys.time(),start,units='secs'))
 if(is.null(b)) b <- list(t_boot=numeric(),t_stat=NA_real_,p_val=NA_real_,boot_iter=NA_real_)
 row <- diagnose_boundary(b,id,'OFFICIAL_PRODUCTION_FIXED_SUPPORT',G,99999L,direction,delta_pp,
    observed_se=NULL,official_division_contract=TRUE,mc_independent_signs_confirmed=TRUE)
 row$r <- r; row$j <- j; row$side <- side; row$contrast <- contrast
 row$started_utc <- format(start,'%Y-%m-%dT%H:%M:%OS6Z',tz='UTC'); row$elapsed_seconds <- elapsed
 row$seed <- 20261015L; row$nthreads <- 1L; row$conf_int <- FALSE
 row$warnings <- paste(unique(warns),collapse=' | '); row$error <- err
 row$observed_t_finite <- is.finite(as.numeric(b$t_stat))
 row$bootstrap_SE_zero_status <- if(row$returned_nonfinite==0L && row$actual_M==row$expected_draws) 'ZERO_EXCLUDED_BY_OFFICIAL_DIVISION_FINITE_T_AND_EXPECTED_M' else 'UNVERIFIED_LOCAL_HOLD'
 row$observed_generated_columns <- NA_integer_
 row$generated_columns_evidence <- 'PINNED_SOURCE_CONTRACT_INFERENCE_NOT_DIRECTLY_MEASURED'
 row$sign_draw_reproducible_reference <- sprintf('336bb574|get_weights|G=%d|M=%d|Rademacher|dqrng|seed=20261015|pref_ascending',G,expected$expected_draws)
 row$production_eligibility_validated <- !row$diagnostic_hold
 # Official p/count consistency is diagnostic; never replace the official p.
 row$official_p_matches_strict_count <- is.finite(row$tail_count) && is.finite(row$official_p_unchanged) && identical(as.numeric(row$official_p_unchanged),as.numeric(row$tail_count/row$actual_M))
 if(!row$official_p_matches_strict_count) {
   row$diagnostic_hold <- TRUE; row$production_eligibility_validated <- FALSE; row$local_status <- 'HOLD'
   row$hold_reasons <- paste(c(row$hold_reasons,'OFFICIAL_P_COUNT_MISMATCH'),collapse=';')
 }
 if(nzchar(err)) row$hold_reasons <- paste(c(row$hold_reasons,'OFFICIAL_CALL_ERROR'),collapse=';')
 # Hash the returned bootstrap statistics through a pipe; no raw array file.
 hashfile <- file.path(out,'boundary_rows',paste0(id,'_t_boot.sha256'))
 hc <- pipe(paste('sha256sum >',shQuote(hashfile)),open='wb')
 writeBin(as.double(b$t_boot),hc,size=8,endian='little'); close(hc)
 row$t_boot_sha256_le_f64 <- strsplit(readLines(hashfile,warn=FALSE)[1],' ')[[1]][1]
 saveRDS(row,file.path(out,'boundary_rows',paste0(id,'.rds')),version=3)
 write_exact_csv(row,path,append=file.exists(path))
 new_count <- new_count+1L
 cat(sprintf('BOUNDARY_DONE contrast=%s id=%s seconds=%.6f M=%d p=%.17g status=%s\n',contrast,id,elapsed,row$actual_M,row$official_p_unchanged,row$local_status)); flush.console()
 rm(b); if(new_count%%20L==0L) gc(FALSE)
}
cat(sprintf('RUN_COMPLETE contrast=%s new_boundaries=%d\n',contrast,new_count)); flush.console()
