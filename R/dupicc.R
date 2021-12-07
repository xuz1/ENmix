dupicc<-function(dat,dupid,mvalue=FALSE,center=TRUE,nCores=2,qcflag=FALSE,qc=NULL,
       detPthre=0.05,nbthre=3,skipicc=FALSE,corfig=FALSE,model="oneway")
{
if(nCores>detectCores()){nCores=detectCores()}
dupid=data.frame(dupid)

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

flag=dupid$id1 %in% colnames(dat) & dupid$id2 %in% colnames(dat) 
dupid=dupid[flag,]

#Shrout & Fleiss (1979), ICC(3,1)
#modified from irr package
#twoway, consistency
icc <-function(meth,model) {
  r0 = 0; conf.level = .95
  meth <- as.matrix(na.omit(meth))
  alpha <- 1-conf.level
  ns <- nrow(meth)
  nr <- ncol(meth)
  SStotal <- var(as.numeric(meth))*(ns*nr-1)
  MSr <- var(apply(meth,1,mean))*nr
  MSw <- sum(apply(meth,1,var)/ns)
  MSc <- var(apply(meth,2,mean))*ns
  MSe <- (SStotal-MSr*(ns-1)-MSc*(nr-1))/((ns-1)*(nr-1))

  if (model=="oneway") {
#Asendorpf & Wallbott, S. 245, ICu; #Bartko (1966)
  coeff  <- (MSr-MSw)/(MSr+(nr-1)*MSw)
  Fvalue <- MSr/MSw*((1-r0)/(1+(nr-1)*r0))
  df1    <- ns-1
  df2    <- ns*(nr-1)
  p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)
 }else if (model=="twoway") {
  coeff  <- (MSr-MSe)/(MSr+(nr-1)*MSe)
  Fvalue <- MSr/MSe*((1-r0)/(1+(nr-1)*r0))
  df1    <- ns-1
  df2    <- (ns-1)*(nr-1)
  p.value <- pf(Fvalue, df1, df2, lower.tail=FALSE)}
  list(value = coeff,Fvalue = Fvalue,p.value = p.value)
}

icc_all <-function(probe,mat1,mat2)
{
        fit=try(icc(meth=cbind(mat1[probe,],mat2[probe,]),model=model))
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
resu=data.frame(id1=colnames(b1),id2=colnames(b2),cor=cc,diff=dist)
resu
}

if(!skipicc){
dupid=dupid[dupid$id1 %in% colnames(dat) & dupid$id2 %in% colnames(dat),]
dat1=dat[,as.vector(dupid$id1)]
dat2=dat[,as.vector(dupid$id2)]
resu1=mclapply(dimnames(dat1)[[1]],icc_all,mat1=dat1,mat2=dat2, mc.cores=nCores)
resu <- do.call(rbind, lapply(resu1, unlist))
resu=as.data.frame(resu)
names(resu)=c("probe","icc","p")
resu$icc=as.numeric(as.vector(resu$icc))
resu$p=as.numeric(as.vector(resu$p))
}

if(center){dat=t(scale(t(dat),center=T,scale=F))}
    dat1=dat[,as.vector(dupid$id1)]
    dat2=dat[,as.vector(dupid$id2)]
    dupcor=dcor(dat1,dat2)
    if(corfig){
    nd1=cor(dat1,method="pearson",use="pairwise.complete.obs");nd1=nd1[lower.tri(nd1)]
    nd2=cor(dat2,method="pearson",use="pairwise.complete.obs");nd2=nd2[lower.tri(nd2)]
    nondupcor=c(nd1,nd2)
      jpeg("dupcorfig.jpg",height=800,width=500,quality=100)
      par(mfrow=c(2,1))
      plot(dupcor$cor,xlab="Duplicate pairs",ylab="Correlation",main="",ylim=c(-1,1))
      plot(nondupcor,xlab="Non-duplicate pairs",ylab="Correlation",main="",ylim=c(-1,1))
      dev.off()
    }

if(!skipicc){list(icc=resu,dupcor=dupcor)}else{dupcor}
}


