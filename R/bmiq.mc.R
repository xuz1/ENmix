bmiq.mc <-function(mdat,nCores=1,...)
{
    if(!is(mdat, "MethylSet")){stop("object needs to be of class 'MethylSet'")}
    if(nCores>detectCores()){
        nCores <- detectCores();
        cat("Only ",nCores," cores are avialable in your computer,", 
           "argument nCores was reset to nCores=",nCores,"\n")
    }
    cat("Analysis is running, please wait...!","\n")

    anno <- getAnnotation(mdat)
    beta.b <- getBeta(mdat, "Illumina")
    beta.b[beta.b <= 0] <- 1e-06
    design.v <- as.vector(anno$Type);
    design.v[design.v == "I"]=1
    design.v[design.v == "II"]=2 
    design.v <- as.numeric(design.v)

    c1 <- makeCluster(nCores)
    registerDoParallel(c1)
    beta.o <- foreach (s = 1:ncol(beta.b),.combine=cbind,.export=c("BMIQ"))  %dopar%{
    s=s;out <- BMIQ(beta.b[, s], design.v, plots = FALSE,sampleID = colnames(beta.b)[s],...)
        return(out$nbeta)
    }
    stopCluster(c1)
    colnames(beta.o) <- colnames(beta.b)
    beta.o
}
