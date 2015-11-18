ComBat.mc <-
function(dat,batch,nCores=1,...)
{
    if(nCores>detectCores()){
        nCores <- detectCores();
        cat("Only ",detectCores()," Cores avialable, nCores was reset to ",detectCores(),"\n")
    }
    cat("Analysis is running, please wait...!","\n")

    rname <- rownames(dat)
    nparts1=min(ceiling(ncol(dat)/500), floor(nrow(dat)/40000))
    parts1 <- rep(1:nparts1,ceiling(nrow(dat)/nparts1))[1:nrow(dat)]
    nCores=min(10, nCores, ceiling(nrow(dat)/(nparts1*20000)))
    dat.o=NULL
    c1 <- makeCluster(nCores)
    registerDoParallel(c1)
    for(i in 1:nparts1){
        dat1=dat[which(parts1==i),]
        parts <- rep(1:nCores,ceiling(nrow(dat1)/nCores))[1:nrow(dat1)]
        parts=sample(parts)
        dat.o1 <- foreach (s = 1:nCores,.combine=rbind,.export=c("ComBat"))  %dopar%{
            s=s;idx=which(parts==s) 
            ComBat(dat=dat1[idx,], batch=batch,...)
        }
        dat.o=rbind(dat.o,dat.o1)
    }
    stopCluster(c1)
    dat.o[rname,]
}


