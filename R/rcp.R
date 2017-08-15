rcp <- function(mdat,dist=25,quantile.grid=seq(0.001,0.999,by=0.001),
    qcscore=NULL,nbthre=3, detPthre=0.000001)
{
    if(!is(mdat, "MethylSet")){stop("the input needs to be of class 
    'MethylSet'")}

    beta<-getB(mdat,type="Illumina")
    raw.M<-logit2(beta)

# find therby pairs of type I probes and type II probes
    annotation<-getAnnotation(mdat)
    annotation=annotation[intersect(rownames(beta),rownames(annotation)),]

    probe.II.Name=annotation$Name[annotation$Type=="II"]
    annotation=annotation[order(annotation$chr,annotation$pos),]
    anno1=annotation[1:(nrow(annotation)-1),]
    anno2=annotation[2:nrow(annotation),]
    flag=(abs(anno1$pos-anno2$pos)<dist & anno1$chr==anno2$chr & 
     anno1$Relation_to_Island==anno2$Relation_to_Island & anno1$Type !=
     anno2$Type)
    anno1=anno1[flag,]
    anno2=anno2[flag,]
    probe.I=anno1$Name
    probe.II=anno2$Name
    probe.I[anno2$Type=="I"]=anno2$Name[anno2$Type=="I"]
    probe.II[anno1$Type=="II"]=anno1$Name[anno1$Type=="II"]

    raw.M.t=raw.M[c(probe.I,probe.II),]

#remove low quality data
    if(is.null(qcscore)){}else if((sum(!(rownames(raw.M.t) %in% 
    rownames(qcscore$detP))) +
    sum(!(colnames(raw.M.t) %in% colnames(qcscore$detP))))>0){
    stop("Wrong qcscore matrix, please check...\n")}else{
    temp <- qcscore$nbead<nbthre | qcscore$detP>detPthre
    temp=temp[rownames(raw.M.t),]
    temp=temp[,colnames(raw.M.t)]
    raw.M.t[temp]=NA
    }

#linear regression
    M.II<-raw.M.t[probe.II,]
    M.I<-raw.M.t[probe.I,]
    
    qtl<-function(x) quantile(x, quantile.grid, na.rm=TRUE)
    M.I=apply(M.I,2,qtl)
    M.II=apply(M.II,2,qtl)

    beta.est<-mat.or.vec(2,ncol(beta))

    for (i in 1:ncol(beta)){
    index<-(M.II[,i]!=Inf & M.II[,i]!=-Inf & M.I[,i]!=Inf & M.I[,i]!=-Inf)
    X<-cbind(rep(1,sum(index)),M.II[index,i]); Y<-M.I[index,i]
    beta.est[,i]<-solve(t(X)%*%X)%*%t(X)%*%Y
    }

    M.II.all<-raw.M[probe.II.Name,]
    M.II.new<-mat.or.vec(nrow(M.II.all),ncol(M.II.all))
    for (i in 1:ncol(M.II.all)){
    M.II.new[,i]<-beta.est[1,i]+beta.est[2,i]*M.II.all[,i]
    }
    M.II.new[M.II.all==Inf]<-Inf; M.II.new[M.II.all==-Inf]<-(-Inf)

    beta[probe.II.Name,]<-ilogit2(M.II.new)
    beta
}


