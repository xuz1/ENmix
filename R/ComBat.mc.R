ComBat.mc <-
function(dat,batch,nCores=1,...)
{
    if(nCores>detectCores()){
        nCores <- detectCores();
        cat("Only ",detectCores()," Cores avialable, nCores was reset to ",detectCores(),"\n")
    }
    cat("Analysis is running, please wait...!","\n")

    rname <- rownames(dat)
    x <- rep(1:nCores,ceiling(nrow(dat)/nCores))[1:nrow(dat)]
    x <- sample(x,length(x),replace=FALSE)
    c1 <- makeCluster(nCores)
    registerDoParallel(c1)
    datc <- foreach (s = 1:nCores,.combine=rbind,.export=c("ComBat"))  %dopar%{
        dat.c=ComBat(dat=dat, batch=batch, ...)
        return(dat.c)
        }
    stopCluster(c1)
    datc[rname,]
}
