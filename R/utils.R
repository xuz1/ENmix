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

getBeta <- function(mdat,type="Illumina",offset=100)
{
    if(!is(mdat, "MethylSet"))
        {stop("The input must be an object of MethylSet\n")}
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


