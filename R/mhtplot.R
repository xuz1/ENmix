
mhtplot <-function(probe=NULL,chr=NULL, pos=NULL, p=NULL, color="bg",sigthre=NULL,markprobe=NULL,markcolor="red",
          outf="mht.jpg")
{
if(is.null(probe)|is.null(chr)|is.null(pos)|is.null(p)){stop("Please provide probe,chr, pos and p\n")}
if(!(length(probe)==length(chr) & length(chr)==length(pos) & length(chr)==length(p))){
     stop("The length of probe, chr, pos and p are different")}
flag=!is.na(probe) & !is.na(chr) & !is.na(pos) & !is.na(p)

dat=data.frame(probe=probe,chr=chr,pos=pos,p=p)
dat=dat[flag,]
tmp=toupper(as.vector(dat$chr))
tmp[tmp %in% "X"]="23"
tmp[tmp %in% "Y"]="24"
dat$chr=as.numeric(tmp)
dat$pos=as.numeric(as.vector(dat$pos))
dat=dat[order(dat$chr,dat$pos),]

chrlen=aggregate(dat$pos,by=list(dat$chr),FUN=max)
chrlen$cumx=cumsum(chrlen$x)
chrlen$mid=(chrlen$cumx-c(0,chrlen$cumx[1:(nrow(chrlen)-1)]))/2+
            c(0,chrlen$cumx[1:(nrow(chrlen)-1)])
for(i in 2:nrow(chrlen))
{dat$pos[dat$chr %in% chrlen$Group.1[i]]=dat$pos[dat$chr %in% chrlen$Group.1[i]]+
            chrlen$cumx[i-1]}

if(color=="bg"){
dat$col="gray"
dat$col[dat$chr %in% chrlen$Group.1[(1:nrow(chrlen) %% 2)==0]]="black"
}else{dat$col=color}

chrlen$Group.1[chrlen$Group.1 == 23]="X"
chrlen$Group.1[chrlen$Group.1 == 24]="Y"

dsig=dat[dat$probe %in% markprobe,]

jpeg(outf,height=600,width=900,quality=100)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(dat$pos,-log10(dat$p),col=dat$col,pch=20,lwd=1,xaxt="n",
    xlab="Chromosome location",
    ylab=expression(paste(-log[10],"(", italic(P), " value)")),
    main="",cex.lab=1.5)

axis(1,at=c(0,chrlen$cumx),las=1,lwd=1,labels=FALSE,cex=0.5,font=2)
graphics::mtext(chrlen$Group.1,side=1,at=chrlen$mid,las=2,adj=0.5,line=1)

if(!is.null(markprobe)){points(dsig$pos,-log10(dsig$p),col=markcolor,pch=20)}

if(!is.null(sigthre)){abline(h=-log10(sigthre),lty=2)}
dev.off()
}







