bmiq.mc <-function(mdat,nCores=1,...)
{
    if(!is(mdat, "MethylSet")){stop("object needs to be of class 'MethylSet'")}
    if(nCores>detectCores()){
        nCores <- detectCores();
        cat("Only ",nCores," cores are avialable in your computer,", 
           "argument nCores was reset to nCores=",nCores,"\n")
    }

    anno <- getAnnotation(mdat)
    cat("Analysis is running, please wait...!","\n")
    beta.b <- getBeta(mdat, "Illumina")
    rm(mdat)
    beta.b[beta.b <= 0] <- 1e-06
    design.v <- as.vector(anno$Type);
    design.v[design.v == "I"]=1
    design.v[design.v == "II"]=2 
    design.v <- as.numeric(design.v)
    coln=colnames(beta.b)

    N=ceiling(ncol(beta.b)/(nCores*10))
    parts=rep(1:N,each = ceiling(ncol(beta.b)/N))[1:ncol(beta.b)]
    c1 <- makeCluster(nCores)
    registerDoParallel(c1)
    for(i in 1:N){
        id=which(parts==i)
        beta.b1=beta.b[,id]

        beta.b1 <- foreach (s = 1:ncol(beta.b1),.combine=cbind,.export=c("BMIQ"))  %dopar%{
            s=s;out <- BMIQ(beta.b1[, s], design.v=design.v, plots = FALSE,...)
            out$nbeta
        }
        beta.b[,id]=beta.b1
    }
    stopCluster(c1)
    if(is.matrix(beta.b)){if(sum(is.na(beta.b))>0){stop("BMIQ estimates encountered error, try
        to run it again")}}else{
         stop("BMIQ estimates encountered error, try to run it again")}
    colnames(beta.b) <- coln
    beta.b
}

