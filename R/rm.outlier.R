rm.outlier <-function(mat,byrow=TRUE,qcscore=NULL,detPthre=0.000001,nbthre=3,
     rmcr=FALSE,rthre=0.05,cthre=0.05,impute=FALSE,imputebyrow=TRUE,...)
{
    if(!(is.numeric(mat) & is.matrix(mat))){stop("Input data must be a
     numeric matrix")}

    #If qcsocre is not NULL, low quality data points will be filtered
    if(is.null(qcscore)){}else if((sum(!(rownames(mat) %in%
     rownames(qcscore$detP))) + 
    sum(!(colnames(mat) %in% colnames(qcscore$detP))))>0){
    stop("Wrong qcscore matrix, please check...\n")}else{
    temp <- qcscore$nbead<nbthre | qcscore$detP>detPthre
    temp=temp[rownames(mat),]
    temp=temp[,colnames(mat)]
    mat[temp]=NA
    }
    #remove outliers
    if(!byrow){mat=t(mat)}
    q2575 <- apply(mat,1,function(x) quantile(x, probs=c(0.25,0.75),
    na.rm=TRUE))
    qr <- q2575["75%",]-q2575["25%",]
    flag=mat < (q2575["25%",] - 3*qr) | mat > (q2575["75%",] + 3*qr)
    mat[flag] <- NA
    if(!byrow){mat=t(mat)}

    #remove rows and columns that have too many missing values
    if(rmcr){
    if(nrow(mat)>ncol(mat)){
    thre=min(2*cthre,0.3)
    cpercna=apply(is.na(mat),2,sum)/nrow(mat)
    tmat=mat[,cpercna<thre]
    rpercna=apply(is.na(tmat),1,sum)/ncol(tmat)
    rthre=max(rthre,min(3,ncol(tmat))/ncol(tmat))
    tmat=mat[rpercna<rthre,]
    nrout=sum(rpercna>=rthre)
    cpercna=apply(is.na(tmat),2,sum)/nrow(tmat)
    cthre=max(cthre,min(3,nrow(tmat))/nrow(tmat))
    mat=tmat[,cpercna<cthre]
    ncout=sum(cpercna>=cthre)
    cat(nrout," rows with percentage of missing data greater than ",rthre,
    " were excluded\n")
    cat(ncout," columns with percentage of missing data greater than ",cthre,
    " were excluded\n")
    }else{
    thre=min(2*rthre,0.3)
    rpercna=apply(is.na(mat),1,sum)/ncol(mat)
    tmat=mat[rpercna<thre,]
    cpercna=apply(is.na(tmat),2,sum)/nrow(tmat)
    cthre=max(cthre,min(3,nrow(tmat))/nrow(tmat))
    tmat=mat[,cpercna<cthre]
    ncout=sum(cpercna>=cthre)
    rpercna=apply(is.na(tmat),1,sum)/ncol(tmat)
    rthre=max(rthre,min(3,ncol(tmat))/ncol(tmat))
    mat=tmat[rpercna<rthre,]
    nrout=sum(rpercna>=rthre)
    cat(nrout," rows with percentage of missing data greater than ",rthre,
    " were excluded\n")
    cat(ncout," columns with percentage of missing data greater than ",cthre,
    " were excluded\n")
    }
    }

    #impute missing data using knn method
    if(impute)
    {
    if(imputebyrow){mat=t(mat)}
    resu=impute.knn(mat,...)

    #error checking and imperfect fix
    rg=apply(mat,2,function(x) range(x,na.rm=TRUE))
    tmat=t(resu$data);gc()
    idx=which(t(tmat<rg[1,] | tmat>rg[2,]),arr.ind = TRUE)
    rm(tmat)
    if(nrow(idx)>0){
       m=apply(mat,2,function(x) mean(x,na.rm=TRUE))
       for(i in 1:nrow(idx)){resu$data[idx[i,1],idx[i,2]]=m[idx[i,2]]}
    }
    mat=resu$data
    if(imputebyrow){mat=t(mat)}
    }
    mat
}

