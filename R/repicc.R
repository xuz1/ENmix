repicc<-function(dat,repid,mvalue=FALSE,nCores=2,qcflag=FALSE,qc=NULL,
       detPthre=0.05,nbthre=3)
{
  if(nCores>detectCores()){nCores=detectCores()}
  if(mvalue){dat=B2M(dat)}
  if(qcflag){
    if(is.null(qc)){stop("Please provide ENmix QCinfo object for quality control")} 
    dat=dat[,!(colnames(dat) %in% qc$badsample)]
    dat=dat[!(rownames(dat) %in% qc$badCpG),]
    qc$detP=qc$detP[rownames(dat),]
    qc$detP=qc$detP[,colnames(dat)]
    qc$nbead=qc$nbead[,colnames(dat)]
    qc$nbead=qc$nbead[rownames(dat),]
    dat[qc$detP>detPthre | qc$nbead<nbthre]<-NA
}      

if(!identical(repid$id,colnames(dat))){
  rownames(repid)=repid$id
  id=intersect(repid$id,colnames(dat))
  repid=repid[id,]
  dat=dat[,id]
#remove singletons
x.char=as.character(repid$idx)
flag=x.char %in% names(table(x.char)[table(x.char)<2])
repid=repid[!flag,]
dat=dat[,!flag]
}

icc<-function(probe,meth,x.char){
  y=meth[probe,]
#remove missing data and singletons
  flag=is.na(y);if(sum(flag)>0){y=y[!flag];x.char=x.char[!flag]
flag=x.char %in% names(table(x.char)[table(x.char)<2]);y=y[!flag];x.char=x.char[!flag]
}
  temp<-aggregate(y~x.char,FUN=mean)
  ga<-temp[,2]
  gs<-table(x.char)
  sst<-sum((ga-mean(y))^2*gs[temp$"x.char"])
  sse<-(sum((y-mean(y))^2)-sst)
  n<-length(y); t<-length(ga)
  mse<-sse/(n-t)
  mst<-sst/(t-1)
  c<-(n^2-sum(gs^2))/(n*(t-1))
  sigma.a.sq<-(mst-mse)/c
  c(probe,ifelse(sigma.a.sq<=0,0,sigma.a.sq/(sigma.a.sq+mse)))
}

resu1=mclapply(dimnames(dat)[[1]],icc,meth=dat,x.char=as.character(repid$idx), mc.cores=nCores)
resu <- do.call(rbind, lapply(resu1, unlist))
resu=as.data.frame(resu)
names(resu)=c("probe","icc")
resu$icc=as.numeric(as.vector(resu$icc))
resu
}


