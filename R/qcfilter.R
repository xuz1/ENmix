qcfilter <-function(mat,qcscore=NULL,rmoutlier=TRUE,byrow=TRUE,detPthre=0.000001,nbthre=3,
     rmcr=FALSE,rthre=0.05,cthre=0.05,impute=FALSE,imputebyrow=TRUE,fastimpute=FALSE,...)
{
    if(impute){requireNamespace("impute")}
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
    if(rmoutlier){
    if(!byrow){mat=t(mat)}
    q2575 <- apply(mat,1,function(x) quantile(x, probs=c(0.25,0.75),
    na.rm=TRUE))
    qr <- q2575["75%",]-q2575["25%",]
    flag=mat < (q2575["25%",] - 3*qr) | mat > (q2575["75%",] + 3*qr)
    mat[flag] <- NA
    if(!byrow){mat=t(mat)}
    }
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

if(fastimpute){
    mm=apply(mat,2,function(x)median(x,na.rm=TRUE))
    idx=which(is.na(mat),arr.ind = TRUE)
    for(i in 1:nrow(idx)){mat[idx[i,1],idx[i,2]]=mm[idx[i,2]]}
}else{
psize=100000
if(ncol(mat)>psize){
    npart=floor(ncol(mat)/psize)
    for(i in 1:npart){
       n1=psize*(i-1)+1
       n2=psize*i
       if(i==npart){n2=ncol(mat)}
       mat1=mat[,n1:n2]
       resu=impute::impute.knn(mat1,...)
       mat[,n1:n2]=resu$data
    }
}else{
      resu=impute::impute.knn(mat,...)
      mat=resu$data
}}

    if(imputebyrow){mat=t(mat)}
    }
    mat
}

