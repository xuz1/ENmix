lmmatrix <- function(y,cov)
{
    p <- matrix(NA,ncol=ncol(cov),nrow=ncol(y))
    for(j in 1:ncol(cov))
    {
    x <- cov[,j]
    for(i in 1:ncol(y))
    {
    fit <- summary(lm(y[,i]~x,na.action=na.omit))
    f <- fit$fstatistic
    p[i,j] <- pf(f["value"],f["numdf"],f["dendf"],lower.tail=FALSE)
    }
    }
    colnames(p) <- names(cov)
    p
}

plotp<-function(jpgfile,p,yaxis,xmax,title)
{
    jpeg(filename=jpgfile,width=1000,height=700,quality = 100)
    par(mar=c(5, 4.2, 4, 2) + 0.1)
    plot(1,xlim=c(0,xmax),ylim=c(0,length(yaxis)+1),type="n",bty="n",
    axes=FALSE,xlab="Principal Component",ylab="",main=title)
    axis(1,at=c(1:xmax),pos=0.5,las=1,lwd=3)
    for(i in 1:length(yaxis)){text(0.3,i, yaxis[i],xpd=TRUE,adj=1)}
    for(i in 1:ncol(p)){
    for(j in 1:nrow(p)){
    pp <- p[j,i]
    colcode <- "white"
    if(pp <= 10e-10){colcode="darkred"}
    else if (pp <= 10e-5){colcode="red"}
    else if (pp <= 0.01){colcode="orange"}
    else if(pp <= 0.05){colcode="pink"}
    polygon(c(j-0.5,j-0.5,j+0.5,j+0.5), c(i-0.5,i+0.5,i+0.5,i-0.5),col=colcode,
    border=NA)
    }}
    legend("topright",c("<0.05","<0.01","<10E-5","<10E-10"),
    col=c("pink","orange","red","darkred"),
    pch=15,pt.cex=2,bty="o",horiz=TRUE,xpd=TRUE)
    dev.off()
}

pcrplot<-function(beta,cov,npc=50)
{
    if(!is.matrix(beta)){stop("beta is not a data matirx")}
    if(!is.data.frame(cov)){stop("cov is not a data frame")}
    if(ncol(beta)!=nrow(cov))
    {stop("number of columns in beta is not equal to number of rows in cov")}
    cat("Analysis is running, please wait...!","\n")
    npc <- min(ncol(beta),npc)
    svd <- prcomp(t(beta),center=TRUE,scale=TRUE,retx=TRUE,na.action="na.omit")
    jpeg(filename="svdscreeplot.jpg",width=1000,height=500,quality = 100)
    screeplot(svd,npc,type="barplot")
    dev.off()
    cat("svdscreeplot.jpg was plotted","\n")
    eigenvalue <- svd[["sdev"]]**2
    prop <- (sum(eigenvalue[1:npc])/sum(eigenvalue))*100
    cat("Top ",npc," principal components can explain ", prop, "% of data 
    variation","\n")
    p <- lmmatrix(svd$x[,1:npc],cov)
    yaxis <- colnames(p)
    plotp("pcr_diag.jpg",p,yaxis,npc,
    title="Principal Component Regression Analysis")
    cat("pcr_diag.jpg was plotted","\n")
}

