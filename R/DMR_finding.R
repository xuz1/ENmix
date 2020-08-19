##### a function to get a table of p-values for estimating acf
#####loc should be increasing; 
acf.table<-function(x,loc,dist.cutoff){
  flag=TRUE; lag=1; result=NULL
  while(flag){
    x1=utils::head(x,-lag); x2=utils::tail(x,-lag); dist=diff(loc,lag=lag)
    index=(dist<dist.cutoff)  
    if(all(!index)){flag=FALSE}else{
      result=rbind(result,data.frame(x1=x1[index],x2=x2[index],dist=dist[index]))
    lag=lag+1
    }
  }
  return(result)  
}

##### a function to estimate acf
get.acf<-function(data,dist.cutoff,bin.size){
  temp<-NULL
  for (chr in unique(data$V1)){
    y<-data[data$V1==chr,]; y<-y[order(y$V3),]
    temp<-rbind(temp,acf.table(y$V4,y$V3,dist.cutoff))
  }
  bin.label<-findInterval(temp$dist,seq(bin.size,dist.cutoff,bin.size))
  temp.stouffer<-by(temp,bin.label,FUN=function(x){cor.test(qnorm(x$x1),
               qnorm(x$x2),alternative="greater")},simplify=FALSE)

  cor.stouffer<-sapply(temp.stouffer,function(x){x$estimate})
  p.stouffer<-sapply(temp.stouffer,function(x){x$p.value})

  if (any(p.stouffer>0.05)){
    index=min(which(p.stouffer>0.05))
    cor.stouffer[index:length(cor.stouffer)]=0
  }
  return(cor.stouffer)
}

#regional P value plot
regplot<-function(ref,sig,extend=2000,outf="region_plot.pdf"){
  sig=sig[order(sig[,1],sig[,2]),]
  ref=ref[order(ref[,1],ref[,2]),]

  pdf(outf)
  for(i in 1:nrow(sig)){
    chr=sig$V1[i]
    pos1=sig$V2[i]
    pos2=sig$V3[i]
    subset=ref[ref$V1==chr & ref$V2>=(pos1-extend) & ref$V2<=(pos2+extend),]
    subset$cor="black"
    subset$cor[subset$V2>=pos1 &subset$V2<=pos2]="red"
    plot(subset$V2,-log10(subset$V4),col=subset$cor,pch=20,xlim=c(pos1-extend,
          pos2+extend),xlab="Chromosome position",ylab=expression(paste(
          -log[10],"(", italic(P), " value)")))
  }
  dev.off()
}

##### interval method
##### 1. get interval p-values; 2. select all intervals with fdr<seed; 
###3. combine nearby sig. intervals and sig. probes (probes with fdr<seed)
### into regions; 4. find region p-values
ipdmr<-function(data,include.all.sig.sites=TRUE,dist.cutoff=1000,bin.size=50,
                seed=0.01,region_plot=TRUE,mht_plot=TRUE,verbose=TRUE){
  data=as.data.frame(data)
  acf<-get.acf(data,dist.cutoff,bin.size)
  if(verbose){
    cat("P value correlations:\n")
    print(data.frame(bin=seq(bin.size,dist.cutoff,bin.size),acf=acf))
  }
  result<-NULL
  for (chr in unique(data$V1)){
    y=data[data$V1==chr,]; y=y[order(y$V3),]
    pos=y$V3; p=qnorm(y$V4)
    index=which(diff(pos)<dist.cutoff)
    int<-findInterval(diff(pos)[index],seq(bin.size,dist.cutoff,bin.size))
    sd<-sqrt(acf[int+1]*2+2)
    int.p<-pnorm(p[index]+p[index+1],mean=0,sd=sd)
    start=pos[index]; end=pos[index+1];
    result=rbind(result,data.frame(chr,start,end,int.p))
  }

  if (include.all.sig.sites){
    temp=data[p.adjust(data$V4,method="fdr")<seed,]
  }else{
    int.sites=unique(c(paste(result$chr,result$start),paste(result$chr,result$end)))
    data.1=data[!(paste(data$V1,data$V3) %in% int.sites),]
    temp=data.1[p.adjust(data.1$V4,method="fdr")<seed,]
  }

  result=result[p.adjust(result$int.p,method="fdr")<seed,]
  result=rbind(result,data.frame(chr=temp$V1,start=temp$V3,end=temp$V3,int.p=temp$V4))

  result.fdr=NULL
  if (nrow(result)>0){
    for (chr in unique(result$"chr")){
      y=data[data$V1==chr,]; y=y[order(y$V3),]
      pos=y$V3; p=qnorm(y$V4)

      result.chr=result[result$"chr"==chr,]
      a=IRanges::IRanges(start=result.chr$start,end=result.chr$end)
      b=IRanges::reduce(a,min.gapwidth=dist.cutoff)

      start=start(b); end=end(b)
      region.max<-max(width(b))
      temp=sapply(1:length(b),function(i){
         index.i=(pos>=start[i] & pos<=end[i]);
         if (sum(index.i)>1){  
           int<-findInterval(c(dist(pos[index.i])),seq(bin.size,region.max+bin.size,bin.size))
           sd<-sqrt(sum(ifelse(int<length(acf),acf[int+1],0))*2+sum(index.i))
           return(pnorm(sum(p[index.i]),mean=0,sd=sd))
         }else{
            return(y$V4[index.i])
         }
      })
      result.fdr=rbind(result.fdr,data.frame(chr,start,end,p=temp))
    }
  
    ##### BH correction
    result.fdr$fdr=p.adjust(result.fdr$p,method="fdr")
    result.fdr<-result.fdr[order(result.fdr$p),]

    ##### use 0-coordinate
    result.fdr$start=(result.fdr$start-1)
  }

  if(is.null(result.fdr)){cat("Number of identified DMR:  0\n")}else{
    ndmr=nrow(result.fdr)
    cat("Number of DMRs identified:  ",ndmr, "\n")
    if(region_plot){
      cat("Drawing regional plot: region_plot.pdf ...\n")
      sig=result.fdr
      colnames(sig)=c("V1","V2","V3","V4","V5")
      regplot(ref=data,sig)
    }
    if(mht_plot){
      cat("Drawing manhattan plot: mht.jpg ...\n")
      set2=c()
      for(i in 1:ndmr){
        set2=c(set2,as.vector(data$V5[data$V1==result.fdr$chr[i] 
              & data$V2>=result.fdr$start[i] & data$V2<=result.fdr$end[i]]))
      }
      mhtplot(data,color="gray",markprobe=set2)
  }
  write.table(result.fdr,"resu_ipdmr.csv",row.names=FALSE,sep=",")
  } 
}

##### comb_p-like method
##### 1. get smoothed p-values; 2. select all probes with smoothed p-values with fdr<seed; 3. combine nearby sig. probes into regions; 4. find region p-values
combp<-function(data,dist.cutoff=1000,bin.size=310,seed=0.01,
               region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE){
  if(nCores>detectCores()){nCores=detectCores()}
  data=as.data.frame(data)
  acf<-get.acf(data,dist.cutoff,bin.size)
  if(verbose){
    cat("P value correlations:\n")
    print(data.frame(bin=seq(bin.size,dist.cutoff,bin.size),acf=acf))
  }

  result<-mclapply(unique(data$V1), function(chr){
    y=data[data$V1==chr,]; y=y[order(y$V3),]
    pos=y$V3; p=qnorm(y$V4)

    temp=sapply(pos,function(i){
      index.i=(abs(pos-i)<bin.size);
      if (sum(index.i)>1){  
        int<-findInterval(c(dist(pos[index.i])),c(bin.size,2*bin.size))
        sd<-sqrt(sum(acf[int+1])*2+sum(index.i))
        return(pnorm(sum(p[index.i]),mean=0,sd=sd))
      }else{return(y$V4[index.i])}
    })

    return(data.frame(chr,start=pos,end=pos,s.p=temp))
  },mc.cores=nCores)

  result <- do.call("rbind", result)
  names(result)=c("chr","start","end","s.p")

  result=result[p.adjust(result$s.p,method="fdr")<seed,]

  result.fdr=NULL
  if (nrow(result)>0){
    for (chr in unique(result$"chr")){
      y=data[data$V1==chr,]; y=y[order(y$V3),]
      pos=y$V3; p=qnorm(y$V4)

      result.chr=result[result$"chr"==chr,]
      a=IRanges::IRanges(start=result.chr$start,end=result.chr$end)
      b=IRanges::reduce(a,min.gapwidth=dist.cutoff)

      start=start(b); end=end(b)
      region.max<-max(width(b))
      temp=sapply(1:length(b),function(i){
        index.i=(pos>=start[i] & pos<=end[i]);
        if (sum(index.i)>1){  
          int<-findInterval(c(dist(pos[index.i])),seq(bin.size,region.max+bin.size,bin.size))
          sd<-sqrt(sum(ifelse(int<length(acf),acf[int+1],0))*2+sum(index.i))
          return(pnorm(sum(p[index.i]),mean=0,sd=sd))
        }else{return(y$V4[index.i])}
      })
      result.fdr=rbind(result.fdr,data.frame(chr,start,end,p=temp))
    }

    ##### BH FDR correction and Sidak correction
    result.fdr$fdr=p.adjust(result.fdr$p,method="fdr")
    result.fdr$sidak=(1-(1-result.fdr$p)^(nrow(data)/(result.fdr$end-result.fdr$start+1)))
    result.fdr<-result.fdr[order(result.fdr$p),]

    ##### use 0-coordinate
    result.fdr$start=(result.fdr$start-1)
  }

  if(is.null(result.fdr)){cat("Number of identified DMR:  0\n")}else{
    ndmr=nrow(result.fdr)
    cat("Number of DMRs identified:  ",ndmr, "\n")
    if(region_plot){
      cat("Drawing regional plot: region_plot.pdf ...\n")
      sig=result.fdr
      colnames(sig)=c("V1","V2","V3","V4","V5","V6")
      regplot(ref=data,sig)
    }
  if(mht_plot){
    cat("Drawing manhattan plot: mht.jpg ...\n")
    set2=c()
    for(i in 1:ndmr){
        set2=c(set2,as.vector(data$V5[data$V1==result.fdr$chr[i]
           & data$V2>=result.fdr$start[i] & data$V2<=result.fdr$end[i]]))
    }
  mhtplot(data,color="gray",markprobe=set2)
  }

  write.table(result.fdr,"resu_combp.csv",row.names=FALSE,sep=",")
  }
}


