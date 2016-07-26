plotCtrl <- function(rgSet,IDorder=NULL)
{
    if(!is(rgSet, "RGChannelSet")){stop("object needs to be of class 
    'RGChannelSet'")}
    if(!is.null(IDorder))
    {
    if(sum(!(IDorder %in% colnames(rgSet)))>0)
    {
       idmissing=IDorder[!(IDorder %in% colnames(rgSet))]
       cat("The IDs were not found in the input data: ",idmissing,"\n")
       stop("Wrong ids in IDorder, please check")
    }else{
    rgSet=rgSet[,IDorder]
    }
    }
    ctrls<-getProbeInfo(rgSet,type="Control")
    ctrls <- ctrls[ctrls$Address %in% featureNames(rgSet),]
    ctrl_r <- getRed(rgSet)[ctrls$Address,]
    ctrl_g <- getGreen(rgSet)[ctrls$Address,]
    contlid <- c("STAINING","EXTENSION","HYBRIDIZATION","TARGET REMOVAL",
     "BISULFITE CONVERSION I", "BISULFITE CONVERSION II","SPECIFICITY I",
     "SPECIFICITY II","NON-POLYMORPHIC","NEGATIVE",
     "NORM_A","NORM_C","NORM_G","NORM_T","NORM_ACGT")
    col <- as.vector(ctrls$Color)
    col[col == "-99"] <- NA
    col[col == "Aqua"] <- "aquamarine2"
    col[col == "Crimson"] <- "firebrick2"
    col[col == "Fuchsia"] <- "deeppink1"
    col[col == "Indigo"] <- "darkviolet"
    col[col == "Lime"] <- "yellowgreen"
    col[col == "Olive"] <- "darkolivegreen"
    col[col == "Silver"] <- "azure4"
    col[col == "Teal"] <- "cyan4"
    ctrls$Color <- col

for(ctype in contlid)
{
    fn <- ctype;fn <- gsub(" ","_",fn)
    cat("Plotting ",fn,".jpg","\n")
    jpeg(paste(fn,".jpg",sep=""),width=1100,height=500,quality=100)
    par(mfrow=c(1,2))
    if(ctype == "NORM_ACGT"){
    cc <- ctrls[ctrls$Type %in% c("NORM_A","NORM_C","NORM_G","NORM_T"),]}
    else{cc <- ctrls[ctrls$Type %in% ctype,]}
    red <- ctrl_r[cc$Address,]
    grn <- ctrl_g[cc$Address,]
    ymax <- max(red,grn)*1.01
if(ctype == "NEGATIVE")
{
    par(mar=c(5, 4, 4, 2))
    colnames(red) <- 1:ncol(red)
    colnames(grn) <- 1:ncol(grn)
    boxplot(grn,ylim=c(0,ymax),main=paste(ctype," Green",sep=""),bty="o",
    xlab="Sample",ylab="Intensity",cex.lab=1.2)
    boxplot(red,ylim=c(0,ymax),main=paste(ctype," Red",sep=""),bty="o",
    xlab="Sample",ylab="Intensity",cex.lab=1.2)
}
else if(ctype %in% c("NORM_A","NORM_C","NORM_G","NORM_T"))
{
    par(mar=c(5, 4, 4, 1))
    colnames(red)<-1:ncol(red)
    colnames(grn)<-1:ncol(grn)
    boxplot(grn,ylim=c(0,ymax),col=unique(as.vector(cc$Color)),main=
    paste(ctype," Green",sep=""),bty="o",xlab="Sample",ylab="Intensity",
    cex.lab=1.2)
    boxplot(red,ylim=c(0,ymax),col=unique(as.vector(cc$Color)),main=
    paste(ctype," Red",sep=""),bty="o",xlab="Sample",ylab="Intensity",
    cex.lab=1.2)
}
else if(ctype == "NORM_ACGT")
{
    par(mar=c(5, 4, 4, 8.2))
    label <- c("NORM_A","NORM_C","NORM_G","NORM_T")
    colcode <- c("Red","Green","Blue","Purple")
    idx <- t(replicate(nrow(red),1:ncol(red)));
    loc <- cc$Color;loc[loc == "Red"] <- -0.2
    loc[loc == "Green"] <- -0.1
    loc[loc == "Blue"] <- 0.1
    loc[loc == "Purple"]=0.2;loc=as.numeric(loc)
    loc1 <- matrix(rep(loc,ncol(red)),ncol=ncol(red));idx=idx+loc1
    plot(idx,grn,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
    " Green",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
    legend(x=(ncol(grn)+0.2)*1.05,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",
    legend=label, col=colcode,xpd=TRUE,pch=15,cex=0.8)
    plot(idx,red,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
    " Red",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
    legend(x=(ncol(red)+0.2)*1.05,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",
    legend=label, col=colcode,xpd=TRUE,pch=15,cex=0.8)
}
else
{
    par(mar=c(5, 4, 4, 8.2))
    idx <- t(replicate(nrow(red),1:ncol(red)))
    plot(idx,grn,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
    " Green",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
    legend(x=ncol(grn)*1.05,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",legend=
    as.vector(cc$ExtendedType),col=as.vector(cc$Color),xpd=TRUE,pch=15,cex=0.8)
    plot(idx,red,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
    " Red",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
    legend(x=ncol(red)*1.045,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",legend=
    as.vector(cc$ExtendedType),col=as.vector(cc$Color),xpd=TRUE,pch=15,cex=0.8)
}
    dev.off()
}
}
