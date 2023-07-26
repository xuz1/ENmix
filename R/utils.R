#convert Beta value to M value
B2M<-function (x)
{
    x[x == 0] <- min(x[x != 0])
    x[x == 1] <- max(x[x != 1])
    log2(x) - log2(1 - x)
}

#convert M value to Beta value
M2B<-function(x)
{2^x/(1+2^x)}

#extract Beta value
getB <- function(mdat,type="Illumina",offset=100)
{
    if(!is(mdat, "methDataSet") & !is(mdat, "MethylSet"))
        {stop("The input is not an object of methDataSet or MethylSet\n")}
    if(type=="Illumina"){offset=100}
    beta<-assays(mdat)$Meth/(assays(mdat)$Meth+assays(mdat)$Unmeth+offset)
    beta
}

#bed example file simulation
simubed <-function(nprobe=1000)
{
   chrset=c(1:22,"X","Y")
   resu=NULL
   for(i in 1:length(chrset)){
     start=sample(1:(90*nprobe),nprobe)
     end=start+1
     p=runif(nprobe)
     resu=rbind(resu,data.frame(chr=chrset[i],start=start,end=end,p=p))
  }
  resu$probe=paste("cg",1:nrow(resu),sep="")
  return(resu)
}

#remove suffix from CpG names, and combine values for duplicated probes
rm.cgsuffix<-function(datMeth){
    cgid=sapply(strsplit(rownames(datMeth),split="_"),unlist)[1,]
    dupcg=unique(cgid[duplicated(cgid)])
    datMeth2=datMeth[cgid %in% dupcg,]
    cid=sapply(strsplit(rownames(datMeth2),split="_"),unlist)[1,]
    datMeth2=aggregate(datMeth2,by=list(cid),FUN=function(x)mean(x,na.rm=TRUE))
    rownames(datMeth2)=datMeth2[,1];datMeth2=as.matrix(datMeth2[,-1])
    datMeth=datMeth[!(cgid %in% dupcg),]
    rownames(datMeth)=sapply(strsplit(rownames(datMeth),split="_"),unlist)[1,]
    rbind(datMeth,datMeth2)
}

#Modified RCP method to normalize datMeth value to a reference
rcp2 <- function(datMeth,reference,quantile.grid=seq(0.001,0.999,by=0.001)){
    #if datMeth is a vector, convert it to a matrix
    if(is.vector(datMeth)){datMeth=as.matrix(datMeth)}

#    quantile.grid=seq(0.001,0.999,by=0.001)
    datMeth<-B2M(datMeth)

    #common CpG set between reference and datMeth
    cg=intersect(as.vector(reference$cg),rownames(datMeth))
    rownames(reference)=reference$cg
    reference=B2M(reference[cg,]$meth_mean)
    datMeth2=datMeth[cg,]
    if(is.vector(datMeth2)){datMeth2=as.matrix(datMeth2)}

#linear regression
    qtl<-function(x) quantile(x, quantile.grid, na.rm=TRUE)
    M.I=qtl(reference)
    M.II=apply(datMeth2,2,qtl)

#    datMeth.est<-mat.or.vec(2,ncol(datMeth))
    datMeth.est<-matrix(NA,2,ncol(datMeth))

    for (i in 1:ncol(datMeth)){
    index<-(M.II[,i]!=Inf & M.II[,i]!=-Inf & M.I!=Inf & M.I!=-Inf)
    X<-cbind(rep(1,sum(index)),M.II[index,i]); Y<-M.I[index]
    datMeth.est[,i]<-solve(t(X)%*%X)%*%t(X)%*%Y
    }

    datMeth.new<-matrix(NA,nrow(datMeth),ncol(datMeth))
    for (i in 1:ncol(datMeth)){
    datMeth.new[,i]<-datMeth.est[1,i]+datMeth.est[2,i]*datMeth[,i]
    }
    datMeth.new[datMeth==Inf]<-Inf; datMeth.new[datMeth==-Inf]<-(-Inf)
    dimnames(datMeth.new)=dimnames(datMeth)
    M2B(datMeth.new)
}




