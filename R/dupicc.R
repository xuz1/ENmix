dupicc<-function(dat,dupid,mvalue=FALSE,center=FALSE,ncores=2,qcflag=FALSE,qc=NULL,detPthre=0.05,nbthre=3,skipicc=FALSE)
{
dupid=data.frame(dupid)

if(mvalue){dat=M2B(dat)}

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

flag=dupid$id1 %in% colnames(dat) & dupid$id2 %in% colnames(dat) 
dupid=dupid[flag,]

icc_all <-function(probe,mat1,mat2)
{
        m1=mat1[probe,]
        m2=mat2[probe,]
        mm=cbind(m1,m2)
	mm=mm[!(is.na(m1) | is.na(m2)),]
        fit=try(icc(mm,type="consistency",model = c("twoway")))
        if(class(fit)[1] == "try-error"){return(c(probe,NA,NA,NA))}else{
        icc=fit$value
        p=fit$p.value
        return(c(probe,icc,p))
        }
}

dcor<-function(b1,b2)
{
cc=array()
dist=array()
for(i in 1:ncol(b1))
{
	x1=b1[,i];x2=b2[,i]
	flag=!(is.na(x1) | is.na(x2))
	x1=x1[flag];x2=x2[flag]
        cc[i]=cor(x1,x2,method="pearson",use="pairwise.complete.obs")
        dist[i]=sum(abs(x1-x2))/length(x1)
}
resu=data.frame(id=colnames(b1),cor=cc,diff=dist)
resu
}

if(!skipicc){
dupid=dupid[dupid$id1 %in% colnames(dat) & dupid$id2 %in% colnames(dat),]
dat1=dat[,as.vector(dupid$id1)]
dat2=dat[,as.vector(dupid$id2)]
resu1=mclapply(dimnames(dat1)[[1]],icc_all,mat1=dat1,mat2=dat2, mc.cores=ncores)
resu <- do.call(rbind, lapply(resu1, unlist))
resu=as.data.frame(resu)
names(resu)=c("probe","icc","p")
resu$icc=as.numeric(as.vector(resu$icc))
resu$p=as.numeric(as.vector(resu$p))
}

if(center){dat=t(scale(t(dat),center=T,scale=F))
	dat1=dat[,as.vector(dupid$id1)]
	dat2=dat[,as.vector(dupid$id2)]
}
	dupcor=dcor(dat1,dat2)

if(!skipicc){list(icc=resu,dupcor=dupcor)}else{dupcor}
}


