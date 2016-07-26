
bincount <- function(x,breaks)
{
    x <- x[!is.na(x)]
    bc <- table(.bincode(x, breaks, TRUE, TRUE))
    temp=data.frame(id=c(1: (length(breaks)-1)))
    bc=data.frame(id=as.numeric(names(bc)),counts=as.numeric(bc))
    resu=merge(temp,bc,by="id",all.x=TRUE)
    resu$counts[is.na(resu$counts)]=0
    resu$counts[order(resu$id)]
}

multifreqpoly <- function(mat, nbreaks=100, col=1:ncol(mat), xlab="", 
    ylab="Frequency", legend = list(x = "top", fill=col,
    legend = if(is.null(colnames(mat))) paste(1:ncol(mat)) else 
    colnames(mat)),...)
{
    if(!is.matrix(mat)) stop("Warning: input data is not a numeric matrix\n")
    if(is.null(col)) col="black"
    col=rep(col,ceiling(ncol(mat)/length(col)))
    if(nbreaks > nrow(mat)) nbreaks=min(15,round(nrow(mat)/2))

    breaks <- seq(min(mat,na.rm=TRUE), max(mat,na.rm=TRUE), 
      diff(range(mat,na.rm=TRUE))/nbreaks)
    mids <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
    counts <- sapply(data.frame(mat),bincount,breaks=breaks)
    plot(range(mids),c(0,max(counts)),type="n",xlab=xlab,ylab=ylab,...)
    for(i in 1:ncol(counts)){lines(mids,counts[,i],col=col[i],...)}
    if(is.list(legend)) do.call(graphics::legend, legend)
}

freqpoly <- function(mat, nbreaks=15, col="black", xlab="", ylab="Frequency",
     type="l",append=FALSE,...)
{
    if(!is.numeric(mat)) stop("Warning: input data is not a numeric vector\n")
    if(nbreaks > length(mat)) nbreaks=min(15,round(length(mat)/2))

    breaks <- seq(min(mat,na.rm=TRUE), max(mat,na.rm=TRUE),
      diff(range(mat,na.rm=TRUE))/nbreaks)
    mids <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
    counts <- bincount(mat,breaks=breaks)
    if(!append){plot(range(mids),c(0,max(counts)),type="n",xlab=xlab,
    ylab=ylab,...)}
    lines(mids,counts,col=col,type=type,...)
}


