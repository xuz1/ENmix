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
  temp.stouffer<-by(temp,bin.label,FUN=function(x){stats::cor.test(qnorm(x$x1),
               qnorm(x$x2),alternative="greater")},simplify=FALSE)

  cor.stouffer<-base::sapply(temp.stouffer,function(x){x$estimate})
  p.stouffer<-base::sapply(temp.stouffer,function(x){x$p.value})

  if (any(p.stouffer>0.05)){
    index=min(which(p.stouffer>0.05))
    cor.stouffer[index:length(cor.stouffer)]=0
  }
  return(cor.stouffer)
}

#regional P value plot
regplot<-function(ref,sig,extend=2000,outf="region_plot.pdf"){
  sig=sig[base::order(sig[,1],sig[,2]),]
  ref=ref[base::order(ref[,1],ref[,2]),]

  grDevices::pdf(outf)
  for(i in 1:base::nrow(sig)){
    chr=sig$V1[i]
    pos1=sig$V2[i]
    pos2=sig$V3[i]
    subset=ref[ref$V1==chr & ref$V2>=(pos1-extend) & ref$V2<=(pos2+extend),]
    subset$cor="black"
    subset$cor[subset$V2>=pos1 &subset$V2<=pos2]="red"
    graphics::plot(subset$V2,-log10(subset$V4),col=subset$cor,pch=20,xlim=c(pos1-extend,
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
    bin=seq(bin.size,dist.cutoff,bin.size)
    if(!(dist.cutoff%%bin.size==0)){bin=base::c(bin,dist.cutoff)}
    print(data.frame(bin=bin,acf=acf))
  }
  result<-NULL
  for (chr in base::unique(data$V1)){
    y=data[data$V1==chr,]; y=y[base::order(y$V3),]
    pos=y$V3; p=stats::qnorm(y$V4)
    index=which(base::diff(pos)<dist.cutoff)
    int<-base::findInterval(base::diff(pos)[index],base::seq(bin.size,dist.cutoff,bin.size))
    sd<-base::sqrt(acf[int+1]*2+2)
    int.p<-stats::pnorm(p[index]+p[index+1],mean=0,sd=sd)
    start=pos[index]; end=pos[index+1];
    result=base::rbind(result,base::data.frame(chr,start,end,int.p))
  }

  if (include.all.sig.sites){
    temp=data[stats::p.adjust(data$V4,method="fdr")<seed,]
  }else{
    int.sites=base::unique(base::c(base::paste(result$chr,result$start),
                  base::paste(result$chr,result$end)))
    data.1=utils::data[!(base::paste(data$V1,data$V3) %in% int.sites),]
    temp=data.1[stats::p.adjust(data.1$V4,method="fdr")<seed,]
  }

  result=result[stats::p.adjust(result$int.p,method="fdr")<seed,]
  result=base::rbind(result,base::data.frame(chr=temp$V1,start=temp$V3,end=temp$V3,int.p=temp$V4))

  result.fdr=NULL
  if (base::nrow(result)>0){
    for (chr in base::unique(result$"chr")){
      y=data[data$V1==chr,]; y=y[base::order(y$V3),]
      pos=y$V3; p=stats::qnorm(y$V4)

      result.chr=result[result$"chr"==chr,]
      a=IRanges::IRanges(start=result.chr$start,end=result.chr$end)
      b=IRanges::reduce(a,min.gapwidth=dist.cutoff)

      start=IRanges::start(b); end=IRanges::end(b)
      region.max<-base::max(IRanges::width(b))
      temp=base::sapply(1:base::length(b),function(i){
         index.i=(pos>=start[i] & pos<=end[i]);
         if (base::sum(index.i)>1){  
           int<-base::findInterval(base::c(stats::dist(pos[index.i])),
                  base::seq(bin.size,region.max+bin.size,bin.size))
           sd<-base::sqrt(base::sum(base::ifelse(int<base::length(acf),acf[int+1],0))*2+base::sum(index.i))
           base::return(stats::pnorm(base::sum(p[index.i]),mean=0,sd=sd))
         }else{
            base::return(y$V4[index.i])
         }
      })
      result.fdr=base::rbind(result.fdr,base::data.frame(chr,start,end,p=temp))
    }
  
    ##### BH correction
    result.fdr$fdr=stats::p.adjust(result.fdr$p,method="fdr")
    result.fdr<-result.fdr[base::order(result.fdr$p),]

    ##### use 0-coordinate
    result.fdr$start=(result.fdr$start-1)
  }

  if(base::is.null(result.fdr)){base::cat("Number of identified DMR:  0\n")}else{
    ndmr=base::nrow(result.fdr)
    base::cat("Number of DMRs identified:  ",ndmr, "\n")
    if(region_plot){
      base::cat("Drawing regional plot: region_plot.pdf ...\n")
      sig=result.fdr
      base::colnames(sig)=c("V1","V2","V3","V4","V5")
      regplot(ref=data,sig)
    }
    if(mht_plot){
      base::cat("Drawing manhattan plot: mht.jpg ...\n")
      set2=base::c()
      for(i in 1:ndmr){
        set2=base::c(set2,base::as.vector(data$V5[data$V1==result.fdr$chr[i] 
              & data$V2>=result.fdr$start[i] & data$V2<=result.fdr$end[i]]))
      }
      mhtplot(probe=data$V5,chr=data$V1,pos=data$V2,p=data$V4,color="gray",markprobe=set2)
  }
  utils::write.table(result.fdr,"resu_ipdmr.csv",row.names=FALSE,sep=",")
  } 
}

##### comb_p-like method
##### 1. get smoothed p-values; 2. select all probes with smoothed p-values with fdr<seed; 3. combine nearby sig. probes into regions; 4. find region p-values
combp<-function(data,dist.cutoff=1000,bin.size=310,seed=0.01,
               region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE){
  if(nCores>parallel::detectCores()){nCores=parallel::detectCores()}
  data=base::as.data.frame(data)
  acf<-get.acf(data,dist.cutoff,bin.size)
  if(verbose){
    base::cat("P value correlations:\n")
    bin=base::seq(bin.size,dist.cutoff,bin.size)
    if(!(dist.cutoff%%bin.size==0)){bin=base::c(bin,dist.cutoff)}
    print(base::data.frame(bin=bin,acf=acf))
  }

  result<-parallel::mclapply(base::unique(data$V1), function(chr){
    y=data[data$V1==chr,]; y=y[base::order(y$V3),]
    pos=y$V3; p=stats::qnorm(y$V4)

    temp=base::sapply(pos,function(i){
      index.i=(base::abs(pos-i)<bin.size);
      if (base::sum(index.i)>1){  
        int<-base::findInterval(base::c(stats::dist(pos[index.i])),base::c(bin.size,2*bin.size))
        sd<-base::sqrt(base::sum(acf[int+1])*2+base::sum(index.i))
        base::return(stats::pnorm(base::sum(p[index.i]),mean=0,sd=sd))
      }else{base::return(y$V4[index.i])}
    })

    base::return(base::data.frame(chr,start=pos,end=pos,s.p=temp))
  },mc.cores=nCores)

  result <- base::do.call("rbind", result)
  names(result)=base::c("chr","start","end","s.p")

  result=result[stats::p.adjust(result$s.p,method="fdr")<seed,]

  result.fdr=NULL
  if (base::nrow(result)>0){
    for (chr in base::unique(result$"chr")){
      y=data[data$V1==chr,]; y=y[base::order(y$V3),]
      pos=y$V3; p=stats::qnorm(y$V4)

      result.chr=result[result$"chr"==chr,]
      a=IRanges::IRanges(start=result.chr$start,end=result.chr$end)
      b=IRanges::reduce(a,min.gapwidth=dist.cutoff)

      start=IRanges::start(b); end=IRanges::end(b)
      region.max<-base::max(IRanges::width(b))
      temp=base::sapply(1:base::length(b),function(i){
        index.i=(pos>=start[i] & pos<=end[i]);
        if (base::sum(index.i)>1){  
          int<-base::findInterval(base::c(stats::dist(pos[index.i])),
              base::seq(bin.size,region.max+bin.size,bin.size))
          sd<-base::sqrt(base::sum(base::ifelse(int<base::length(acf),
              acf[int+1],0))*2+sum(index.i))
          base::return(stats::pnorm(base::sum(p[index.i]),mean=0,sd=sd))
        }else{base::return(y$V4[index.i])}
      })
      result.fdr=base::rbind(result.fdr,base::data.frame(chr,start,end,p=temp))
    }

    ##### BH FDR correction and Sidak correction
    result.fdr$fdr=stats::p.adjust(result.fdr$p,method="fdr")
    result.fdr$sidak=(1-(1-result.fdr$p)^(base::nrow(data)/(result.fdr$end-result.fdr$start+1)))
    result.fdr<-result.fdr[base::order(result.fdr$p),]

    ##### use 0-coordinate
    result.fdr$start=(result.fdr$start-1)
  }

  if(is.null(result.fdr)){cat("Number of identified DMR:  0\n")}else{
    ndmr=base::nrow(result.fdr)
    base::cat("Number of DMRs identified:  ",ndmr, "\n")
    if(region_plot){
      cat("Drawing regional plot: region_plot.pdf ...\n")
      sig=result.fdr
      colnames(sig)=c("V1","V2","V3","V4","V5","V6")
      regplot(ref=data,sig)
    }
  if(mht_plot){
    base::cat("Drawing manhattan plot: mht.jpg ...\n")
    set2=base::c()
    for(i in 1:ndmr){
        set2=base::c(set2,base::as.vector(data$V5[data$V1==result.fdr$chr[i]
           & data$V2>=result.fdr$start[i] & data$V2<=result.fdr$end[i]]))
    }
  mhtplot(probe=data$V5,chr=data$V1,pos=data$V2,p=data$V4,color="gray",markprobe=set2)
  }

  utils::write.table(result.fdr,"resu_combp.csv",row.names=FALSE,sep=",")
  }
}


