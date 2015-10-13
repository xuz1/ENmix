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

    if(.Platform$OS.type == "unix"){
        BMIQ.b <-function(ID,beta.b,design.v,...)
        {
            out<-BMIQ(beta.b[,ID], design.v, plots = FALSE,...)
            out$nbeta
        }
        beta.o=mclapply(1:ncol(beta.b),BMIQ.b,beta.b=beta.b,design.v=design.v,...,
               mc.silent=TRUE,mc.cores=nCores)
        beta.o <- do.call(cbind, lapply(beta.o, unlist))
    }else{
        c1 <- makeCluster(nCores)
        registerDoParallel(c1)
        beta.o <- foreach (s = 1:ncol(beta.b),.combine=cbind,.export=c("BMIQ"))  %dopar%{
        s=s;out <- BMIQ(beta.b[, s], design.v=design.v, plots = FALSE,...)
            return(out$nbeta)
        }
        stopCluster(c1)
    }
    colnames(beta.o) <- colnames(beta.b)
    beta.o
}

