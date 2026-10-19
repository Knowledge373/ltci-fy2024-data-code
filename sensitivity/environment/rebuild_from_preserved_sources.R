# Run with ./run-r-local.sh -f rebuild_from_preserved_sources.R from runtime directory.
# Installs exact preserved versions in dependency order, without querying latest versions.
order <- c('BH','lattice','Formula','generics','gtools','JuliaConnectoR','Rcpp','rlang','sitmo','stringmagic','collapse','Matrix','RcppArmadillo','RcppEigen','dqrng','dreamerr','MASS','cli','summclust','fwildclusterboot')
stopifnot(getRversion()=='4.6.1')
for (p in order) {
 f <- list.files('packages',pattern=paste0('^',p,'_[^/]+[.]tar[.]gz$'),full.names=TRUE)
 stopifnot(length(f)==1L)
 install.packages(f,repos=NULL,type='source',Ncpus=1)
 stopifnot(requireNamespace(p,quietly=TRUE))
}
stopifnot(packageVersion('fwildclusterboot')=='0.14.3',packageVersion('dqrng')=='0.4.1')
write.csv(installed.packages()[,c('Package','Version','Built')], 'INSTALLED_PACKAGES.csv', row.names=FALSE)
capture.output(sessionInfo(),file='SESSION_INFO.txt')
