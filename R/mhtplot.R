

mhtplot <-function(dat,color="bg",sigthre=NULL,markprobe=NULL,markcolor="red",
          outf="mht.jpg")
{
tmp=toupper(as.vector(dat$V1))
tmp[tmp %in% "X"]="23"
tmp[tmp %in% "Y"]="24"
dat$V1=as.numeric(tmp)
dat$V2=as.numeric(as.vector(dat$V2))
dat=dat[order(dat$V1,dat$V2),]

chrlen=aggregate(dat$V2,by=list(dat$V1),FUN=max)
chrlen$cumx=cumsum(chrlen$x)
chrlen$mid=(chrlen$cumx-c(0,chrlen$cumx[1:(nrow(chrlen)-1)]))/2+
            c(0,chrlen$cumx[1:(nrow(chrlen)-1)])
for(i in 2:nrow(chrlen))
{dat$V2[dat$V1 %in% chrlen$Group.1[i]]=dat$V2[dat$V1 %in% chrlen$Group.1[i]]+
            chrlen$cumx[i-1]}
chrlen$Group.1[chrlen$Group.1 == 23]="X"
chrlen$Group.1[chrlen$Group.1 == 24]="Y"

if(color=="bg"){
dat$col="gray"
dat$col[dat$V1 %in% chrlen$Group.1[(1:nrow(chrlen) %% 2)==0]]="black"
}else{dat$col=color}

dsig=dat[dat$V5 %in% markprobe,]

jpeg(outf,height=600,width=900,quality=100)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(dat$V2,-log10(dat$V4),col=dat$col,pch=20,lwd=1,xaxt="n",
    xlab="CpG physical location (by chromosome)",
    ylab=expression(paste(-log[10],"(", italic(P), " value)")),
    main="",cex.lab=1.5)

axis(1,at=c(0,chrlen$cumx),las=1,lwd=1,labels=FALSE,cex=0.5,font=2)
graphics::mtext(chrlen$Group.1,side=1,at=chrlen$mid,las=2,adj=0.5,line=1)

if(!is.null(markprobe)){points(dsig$V2,-log10(dsig$V4),col=markcolor,pch=20)}

if(!is.null(sigthre)){abline(h=-log10(sigthre),lty=2)}
dev.off()
}







