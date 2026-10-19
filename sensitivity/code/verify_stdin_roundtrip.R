# No model or bootstrap: verify only CSV float parsing and memory byte identity.
args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==1L)
out <- args[1]
d <- read.csv(file('stdin'),check.names=FALSE)
X <- cbind(1,as.matrix(d[setdiff(names(d),c('y','pref'))]))
write_hash <- function(v,name,size) {
 p <- file.path(out,paste0('R_STDIN_',name,'.sha256'))
 h <- pipe(paste('sha256sum >',shQuote(p)),open='wb')
 writeBin(v,h,size=size,endian='little'); close(h)
}
write_hash(as.double(t(X)),'X_LE_F64_C',8L)
write_hash(as.double(d$y),'Y_LE_F64',8L)
write_hash(as.integer(d$pref),'CLUSTER_LE_I32',4L)
cat('STDIN_HASH_COMPLETE\n')
