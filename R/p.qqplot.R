
p.qqplot<-function(pvalues,outf="qq",outfmt="jpg",draw.conf=TRUE,
        conf.col="lightgray",conf.alpha=.95,pch=20,col="black",reducesize=0,...)
{

    conf.alpha=1-conf.alpha
    ## Calculate lambda
    lambda <- qchisq(median(pvalues,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
    cat("\nThe genomic inflation factor lambda=",lambda,"\n")

    pvalues=pvalues[!is.na(pvalues)]
    if(!is.numeric(pvalues)){stop("pvalues is not a numeric vector")}
    if(min(pvalues)<=0){stop("Some pvalues is less or equal to 0")}
    if(max(pvalues)>1){stop("Some pvalues is greater than 1")}

    pvalues=sort(pvalues)
    n=length(pvalues)+1
    exp.x <- -log10((1:(n-1))/n)
    pvalues <- -log10(pvalues)

    if(draw.conf){
        conf.points = n-1;
        cpts<-matrix(nrow=conf.points*2, ncol=2)
        for(i in seq(from=1, to=conf.points)) {
           cpts[i,1]<- -log10((i)/n)
           cpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
           cpts[conf.points*2+1-i,1]<- -log10((i)/n)
           cpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
        }
    }

#    xlab=expression(Expected~~-log[10](italic(P))~~"under NULL hypothesis")
#    ylab=expression(Observed~~-log[10](italic(P)))

     xlab=bquote('Expected  '*'-log'['10']*'(P)'*'  under NULL hypothesis')
     ylab=bquote('Observed  '*'-log'['10']*'(P)')

if(outfmt == "jpg"){
jpeg(paste(outf,".",outfmt,sep=""),height=600,width=600,quality=100)
par(mar=c(5.1, 5.1, 4.1, 2.1))
}else if(outfmt == "eps"){
grDevices::postscript(paste(outf,".",outfmt,sep=""),width=3.5,height=3.5,
    paper = "special",pointsize=5,onefile=FALSE)
par(mar=c(5.1, 5.1, 4.1, 2.1))
while(reducesize>0){
dat=data.frame(ex=exp.x,p=pvalues)
dat=dat[order(dat$ex),]
dat$id=1:nrow(dat)
dat$v1=dat$ex/max(dat$ex)
dat$v2=dat$p/max(dat$p)
dat$v3=c(2,abs(dat$v1[2:nrow(dat)]-dat$v1[1:(nrow(dat)-1)])+
  abs(dat$v2[2:nrow(dat)]-dat$v2[1:(nrow(dat)-1)]))
dat$V3[nrow(dat)]=2
dat=dat[order(dat$v3),]
exid=dat$id[1:round(nrow(dat)/2)]
exid=exid[exid%%2==0]
dat=dat[!(dat$id %in% exid),]
reducesize=reducesize-1
exp.x=dat$ex
pvalues=dat$p
}
}else{stop("Outfmt should be jpg or eps\n")}

    plot(c(0,0),type="n",xlim=c(0,max(exp.x)), ylim=c(0,max(pvalues)), xlab=xlab,ylab=ylab,cex.lab=1.3,...)
    if(draw.conf){polygon(x=cpts[,1],y=cpts[,2], col="lightgray", lty=0)}
    points(exp.x,pvalues,pch=pch,col=col)
    abline(0,1,lty=2)
    dev.off()
}



