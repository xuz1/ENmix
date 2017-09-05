#convert Beta value to M value
B2M<-function (x)
{
    x[x == 0] <- min(x[x != 0])
    x[x == 1] <- max(x[x != 1])
    log(x/(1 - x))
}

#convert M value to Beta value
M2B<-function(x)
{exp(x)/(1+exp(x))}

#extract Beta value
getB <- function(mdat,type="Illumina",offset=100)
{
    if(!is(mdat, "MethylSet"))
        {stop("The input must be an object of MethylSet\n")}
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


