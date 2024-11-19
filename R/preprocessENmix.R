

log.erfc  <- function(x){
    log(2)+pnorm(x*sqrt(2),lower.tail=FALSE,log.p=TRUE)
}

loglik.x  <- function(x,background,lambda,mu2,sigma2,p1,p2,mu,sigma){
    if (background){
    temp.1 <- lambda*(2*mu+lambda*sigma^2-2*x)/2+log.erfc((mu+lambda*
    sigma^2-x)/(sqrt(2)*sigma))
    temp.2 <- dnorm(x,mu+mu2,sqrt(sigma^2+sigma2^2),log=TRUE)
    return(sum(ifelse(temp.1>temp.2,temp.1+log(p1*lambda/2+p2*exp(temp.2-
    temp.1)),temp.2+log(p1*lambda/2*exp(temp.1-temp.2)+p2))))
    }
    else{
    temp.1 <- dexp(x,lambda,log=TRUE)
    temp.2 <- dnorm(x,mu2,sigma2,log=TRUE)
    return(sum(ifelse(temp.1>temp.2,temp.1+log(p1+p2*exp(temp.2-temp.1)),
    temp.2+log(p1*exp(temp.1-temp.2)+p2))))
    }
}

EM_estimate  <- 
function(x,start=c(max(density(x)$y),mean(range(x)),diff(range(x))/6,0.5),
    epsilon=c(0.0001,0.001,0.001,0.001)){
    lambda <- start[1]
    mu <- start[2]
    sigma.sq <- start[3]^2
    p1 <- start[4]
    n <- length(x)

    diff <- TRUE
    while (diff){
    z <- p1/(p1+(1-p1)*exp(dnorm(x,mu,sqrt(sigma.sq),log=TRUE)-dexp(x,
    lambda,log=TRUE)))

    lambda.new <- sum(z)/sum(x*z)
    mu.new <- sum((1-z)*x)/sum(1-z)
    sigma.sq.new <- sum(((x-mu.new)^2)*(1-z))/sum(1-z)
    p1.new <- sum(z)/n

    diff <- !(all(c(abs(lambda.new-lambda),abs(mu.new-mu),
    abs(sigma.sq.new-sigma.sq),abs(p1.new-p1))<epsilon))

    lambda <- lambda.new; mu <- mu.new; sigma.sq <- sigma.sq.new; p1 <- p1.new
    }
    list(lambda,mu,sqrt(sigma.sq),p1)
}

new_cm  <- 
function(lambda,mu2,sigma2,p1,p2,mu,sigma,s){
    a <- mu+lambda*sigma^2
    C <- (sigma2^2*(s-mu)+sigma^2*mu2)/(sigma2^2+sigma^2)
    D <- sigma2*sigma/sqrt(sigma2^2+sigma^2)

    temp_1 <- p1*lambda*exp(lambda^2*sigma^2/2-lambda*(s-mu))*(pnorm(s,a,sigma)-
     pnorm(0,a,sigma))
    temp_2 <- p2*dnorm(s,mu2+mu,sqrt(sigma2^2+sigma^2))/(1-pnorm(0,mu2,sigma2))

    num_1 <- temp_1
    num_2 <- temp_2*(pnorm(s,C,D)-pnorm(0,C,D))

    denom_1 <- temp_1*(s-a+sigma*(dnorm((s-a)/sigma)-dnorm(a/sigma))/(
     pnorm(s,a,sigma)-pnorm(0,a,sigma)))
    denom_2 <- temp_2*(C*(pnorm(s,C,D)-pnorm(0,C,D))+D*(exp(-C^2/(2*D^2))
     -exp(-(s-C)^2/(2*D^2)))/sqrt(2*pi))

    (denom_1+denom_2)/(num_1+num_2)
}

firstpeak <- function(x,y,sn,dat)
{
    n <- length(y)
##sn to avoid using small peak position
    nn <- sn*2+1
    v <- matrix(NA,ncol=nn,nrow=n-nn+1)
    for(i in 1:nn){v[,i]=y[i:(n-nn+i)]}
    ix <- sn+which(apply(v<v[,(sn+1)],1,sum) == (nn-1))
    if(length(ix)>0){mu=x[ix[1]]}
    if(length(ix) == 0 | sum(dat<mu)/length(dat)>=0.15)
    {
    dat <- dat[dat<quantile(dat,p=0.10)]
    temp <- density(dat)
    flag <- temp$x>=min(dat) & temp$x<=max(dat)
    temp$x <- temp$x[flag];temp$y=temp$y[flag]
    mu <- temp$x[which.max(temp$y)]
    }
    mu
}


enmix_adj <- function(meth_i=NULL,bg_i=NULL,bgParaEst)
{
    if(sum(is.na(meth_i))>0)
    {stop("ENmix background correction does not allow missing value")}
    meth_i[meth_i <= 0]=1
    mu <- bg_i$mu[1]
    sigma <- bg_i$sigma[1]
if(bgParaEst == "est" | bgParaEst == "neg" | bgParaEst == "oob")
{
    x <- (meth_i[meth_i>=mu]-mu)
    temp <- EM_estimate(x)
    lambda <- temp[[1]]
    mu2 <- temp[[2]]
    sigma2<-ifelse(temp[[3]] <= sigma, 0.1, sqrt(temp[[3]]^2-sigma^2))

    p1 <- (sum(meth_i<mu)+temp[[4]]*length(x))/length(meth_i)
    p2 <- 1-p1
    meth_adj <- new_cm(lambda,mu2,sigma2,p1,p2,mu,sigma,meth_i)
    meth_adj[meth_adj <= 0]=0.01 #restrict to positive values, only a few
}else if (bgParaEst == "subtract_neg" | bgParaEst == "subtract_estBG"
    | bgParaEst == "subtract_q5neg" | bgParaEst == "subtract_oob"){
    meth_adj=meth_i-mu
    meth_adj[meth_adj <= 0]=0.01 #restrict to positive values
}
    meth_adj
}

enmix <- function(meth,bg,bgParaEst,nCores)
{
    colnm <- colnames(meth)
    meth.o <- foreach(i=1:ncol(meth),.combine=cbind,.export=c("EM_estimate",
     "new_cm","enmix_adj")) %dopar% {
      i=i;enmix_adj(meth[,i],bg[i,],bgParaEst)}
    if(is.matrix(meth.o)){if(sum(is.na(meth.o))>0){stop("Computation ran out
    of memory, try to set nCores with a smaller value")}}else{
     stop("Computation ran out of memory, try to set nCores with a smaller
      value")}
    colnames(meth.o)=colnm
    gc(); meth.o
}

#copy from MASS package
huber<-function (y, k = 1.5, tol = 1e-06)
{
    y <- y[!is.na(y)]
    n <- length(y)
    mu <- median(y)
    s <- mad(y)
    if (s == 0)
        stop("cannot estimate scale: MAD is zero for this sample")
    repeat {
        yy <- pmin(pmax(mu - k * s, y), mu + k * s)
        mu1 <- sum(yy)/n
        if (abs(mu - mu1) < tol * s)
            break
        mu <- mu1
    }
    list(mu = mu, s = s)
}

huber_mus <- function(x){ests <- try(huber(x)); if(class(ests)[1]=="try-error"){
    cat("Warning:Check negtive control data, or do quality control before ENmix\n");
    c(mu=median(x,na.rm=TRUE),s=sd(x,na.rm=TRUE))
    }else{c(mu=ests$mu,s=ests$s)}}
huber_mu <- function(x){ests <- try(huber(x));if(class(ests)[1]=="try-error"){
    cat("Warning: Check NORM control data, or do quality control before ENmix\n");
    median(x,na.rm=TRUE)
    }else{ests$mu}}

estBG  <- function(meth_i)
{
    meth_i[meth_i<=0]=1e-06
    temp <- density(meth_i)
    temp <- density(meth_i[meth_i<temp$x[which.max(temp$y)]])
    flag <- temp$x>=min(meth_i) & temp$x<=max(meth_i)
    temp$x <- temp$x[flag];temp$y=temp$y[flag]
    mu <- temp$x[which.max(temp$y)]
    ##first mode
    if((sum(meth_i<mu)/length(meth_i))>=0.15){mu=firstpeak(temp$x,temp$y,
      sn=5,meth_i)}
    perc <- sum(meth_i<mu)/length(meth_i)
    sigma <- sqrt(sum((meth_i[meth_i<mu]-mu)^2)/sum(meth_i<mu))
    c(mu,sigma,perc)
}

##background correction
preprocessENmix  <- function(rgSet, bgParaEst="oob", dyeCorr="RELIC",
    QCinfo=NULL, exQCsample=TRUE,
    exQCcpg=TRUE, exSample=NULL, exCpG=NULL, nCores=2)
{
    if(is(rgSet, "rgDataSet")){
    if(!is.null(QCinfo)){exSample=unique(c(QCinfo$badsample, exSample))}
    exSample=exSample[exSample %in% colnames(rgSet)]
    if(length(exSample)>0){
    rgSet=rgSet[,!(colnames(rgSet) %in% exSample)]
    cat(length(exSample), " samples were excluded before ENmix correction\n")
    }
    mdat <- getmeth(rgSet)
    }else if(is(rgSet, "RGChannelSet")){
    if(!is.null(QCinfo)){exSample=unique(c(QCinfo$badsample, exSample))}
    exSample=exSample[exSample %in% colnames(rgSet)]
    if(length(exSample)>0){
    rgSet=rgSet[,!(colnames(rgSet) %in% exSample)]
    cat(length(exSample), " samples were excluded before ENmix correction\n")
    }
    mdat <- preprocessRaw(rgSet)
    }else if(is(rgSet, "methDataSet") | is(rgSet, "MethylSet")){
    if(!is.null(QCinfo) & exQCsample){exSample=unique(c(QCinfo$badsample,
     exSample))}
    exSample=exSample[exSample %in% colnames(rgSet)]
    if(length(exSample)>0){
    rgSet=rgSet[,!(colnames(rgSet) %in% exSample)]
    cat(length(exSample), " samples were excluded before ENmix correction\n")
    }
    mdat=rgSet; bgParaEst="est"; 
    if(!(dyeCorr=="none")){cat("Warning: Input data need to be a rgDataSet or 
        RGChannelSet to perform dye-bias correction\n");
       cat("Warning: dye-bias correction will not be performed\n")}
    dyeCorr="none"
    }else{stop("Error: object needs to be of class 'RGChannelSet' or 
      'MethylSet'")}
    if(nCores>detectCores()){
    nCores=detectCores(); 
    cat("Only ",detectCores(), " cores avialable, nCores was reset to ",
     detectCores(),"\n")
    }
    if(!is.null(QCinfo) & exQCcpg) {exCpG=unique(c(exCpG, QCinfo$badCpG))}
    exCpG=exCpG[exCpG %in% rownames(mdat)]
    if(length(exCpG)>0){
    mdat=mdat[!(rownames(mdat) %in% exCpG),]
    cat(length(exCpG), " CpGs were excluded before ENmix correction\n")
    }
    rm(QCinfo)
    
    if(is(mdat, "methDataSet")){
      probe_type=rowData(mdat)$Infinium_Design_Type
      col=rowData(mdat)$Color_Channel
      probe_type[probe_type %in% c("I","snpI") & col=="Grn"]="IGrn"
      probe_type[probe_type %in% c("I","snpI") & col=="Red"]="IRed"
      probe_type[probe_type %in% c("snpII")]="II"
    }else if(is(mdat, "MethylSet")){
    probe_type <- getProbeType(mdat, withColor=TRUE)}
    cat("Analysis is running, please wait...!\n")
##estimate background parameters
    if(bgParaEst == "neg" | bgParaEst == "subtract_neg")
    {
    if(is(rgSet,"rgDataSet")){
        ctrls <- getCGinfo(rgSet,type="ctrl")
    }else if(is(rgSet,"RGChannelSet")){
        ctrls <- getProbeInfo(rgSet,type="Control")}

    ctrls <- ctrls[ctrls$Address %in% rownames(rgSet),]
    ctrl_address <- as.vector(ctrls$Address[ctrls$Type %in% "NEGATIVE"])

    ctrl_r <- assays(rgSet)$Red[ctrl_address,]
    ctrl_g <- assays(rgSet)$Green[ctrl_address,]
    ctrl_r[ctrl_r<=0]=1e-06;ctrl_g[ctrl_g<=0]=1e-06
    temp <- apply(ctrl_r,2,huber_mus)
    mu <- temp["mu",];sigma <- temp["s",]
    bgRI <- as.data.frame(cbind(mu,sigma))
    temp <- apply(ctrl_g,2,huber_mus);
    mu <- temp["mu",];sigma <- temp["s",]
    bgGI <- as.data.frame(cbind(mu,sigma))
    bgRII <- bgRI;bgGII <- bgGI
    }else if (bgParaEst == "oob" | bgParaEst == "subtract_oob")
    {
    if(is(rgSet,"rgDataSet")){
        I_probe <- getCGinfo(rgSet,type="I")
        I_probe=I_probe[I_probe$Color_Channel=="Red",]
    }else if(is(rgSet,"RGChannelSet")){
        I_probe <- getProbeInfo(rgSet, type="I-Red")}
    I_green_bg_M <-  assays(rgSet)$Green[I_probe$AddressB,]
    I_green_bg_U <-  assays(rgSet)$Green[I_probe$AddressA,]
    ctrl_g <- rbind(I_green_bg_M,I_green_bg_U)
    if(is(rgSet,"rgDataSet")){
        I_probe <- getCGinfo(rgSet,type="I")
        I_probe=I_probe[I_probe$Color_Channel=="Grn",]
    }else if(is(rgSet,"RGChannelSet")){
        I_probe <- getProbeInfo(rgSet, type="I-Green")}
    I_red_bg_M <-  assays(rgSet)$Red[I_probe$AddressB,]
    I_red_bg_U <-  assays(rgSet)$Red[I_probe$AddressA,]
    ctrl_r <- rbind(I_red_bg_M,I_red_bg_U)
    ctrl_r[ctrl_r<=0]=1e-06;ctrl_g[ctrl_g<=0]=1e-06
    temp <- apply(ctrl_r,2,huber_mus)
    mu <- temp["mu",];sigma=temp["s",]
    bgRI <- as.data.frame(cbind(mu,sigma))
    temp <- apply(ctrl_g,2,huber_mus);
    mu <- temp["mu",];sigma=temp["s",]
    bgGI <- as.data.frame(cbind(mu,sigma))
    bgRII <- bgRI;bgGII=bgGI
    rm(list=c("I_green_bg_M","I_green_bg_U","ctrl_g","I_red_bg_M","I_red_bg_U",
     "ctrl_r"))

    }else if (bgParaEst == "subtract_q5neg"){
    if(is(rgSet,"rgDataSet")){
        ctrls <- getCGinfo(rgSet,type="ctrl")
    }else if(is(rgSet,"RGChannelSet")){
        ctrls <- getProbeInfo(rgSet,type="Control")}

    ctrls <- ctrls[ctrls$Address %in% rownames(rgSet),]
    ctrl_address <- as.vector(ctrls$Address[ctrls$Type %in% "NEGATIVE"])
    ctrl_r <- assays(rgSet)$Red[ctrl_address,]
    ctrl_g <- assays(rgSet)$Green[ctrl_address,]
    ctrl_r[ctrl_r<=0]=1e-06;ctrl_g[ctrl_g<=0]=1e-06  ## may not need this
    mu <- apply(ctrl_r,2,function(x) quantile(x,probs=0.05,na.rm=TRUE));
    sigma <- apply(ctrl_r,2,function(x)sd(x,na.rm=TRUE));
    bgRI <- as.data.frame(cbind(mu,sigma))
    mu <- apply(ctrl_g,2,function(x) quantile(x,probs=0.05,na.rm=TRUE));
    sigma <- apply(ctrl_g,2,function(x)sd(x,na.rm=TRUE));
    bgGI <- as.data.frame(cbind(mu,sigma))
    bgRII <- bgRI;bgGII=bgGI
    }else if(bgParaEst == "est" | bgParaEst == "subtract_estBG"){
    
    mdat_subset <- mdat[probe_type == "IRed",]
    m_I_red <- rbind(assays(mdat_subset)$Meth,assays(mdat_subset)$Unmeth)
    mdat_subset <- mdat[probe_type == "IGrn",]
    m_I_grn <- rbind(assays(mdat_subset)$Meth,assays(mdat_subset)$Unmeth)
    mdat_subset <- mdat[probe_type == "II",]
    mII <- rbind(assays(mdat_subset)$Meth,assays(mdat_subset)$Unmeth)
    rm(mdat_subset)
    bgRI <- as.data.frame(t(apply(m_I_red,2,estBG)));names(bgRI) <- c("mu",
     "sigma","perc")
    bgGI <- as.data.frame(t(apply(m_I_grn,2,estBG)));names(bgGI) <- c("mu",
     "sigma","perc")
    bgII <- as.data.frame(t(apply(mII,2,estBG)));names(bgII) <- c("mu",
     "sigma","perc")
    ##empirically adjusting the estimates
    pp <- apply(cbind(bgRI$perc,bgGI$perc,bgII$perc),1,max)-
     apply(cbind(bgRI$perc,bgGI$perc,bgII$perc),1,min)
    avgp <- apply(cbind(bgRI$perc,bgGI$perc,bgII$perc),1,mean)
    for(i in 1:nrow(bgGI)){if(pp[i]>=0.04){
    bgRI$mu[i] <- quantile(m_I_red[,i],probs=avgp[i])
    bgGI$mu[i] <- quantile(m_I_grn[,i],probs=avgp[i])
    bgII$mu[i] <- quantile(mII[,i],probs=avgp[i])
    bgRI$perc[i] <- avgp[i];bgGI$perc[i]=avgp[i];bgII$perc[i]=avgp[i]
    }}
    bgRI <- bgRI[,c("mu","sigma")]
    bgGI <- bgGI[,c("mu","sigma")]
    bgII <- bgII[,c("mu","sigma")]
    A1=sum(bgII$mu)*2/(sum(bgRI$mu)+sum(bgGI$mu))
    A2=sum(bgII$sigma)*2/(sum(bgRI$sigma)+sum(bgGI$sigma))
    bgGII=bgGI;bgGII$mu=bgGII$mu*A1;bgGII$sigma=bgGII$sigma*A2
    bgRII=bgRI;bgRII$mu=bgRII$mu*A1;bgRII$sigma=bgRII$sigma*A2
    rm(list=c("m_I_red","m_I_grn","mII"))
    }
    c1 <- makeCluster(nCores)
    registerDoParallel(c1)

    if (dyeCorr == "mean"){
    if(is(rgSet,"rgDataSet")){
        ctrls <- getCGinfo(rgSet,type="ctrl")
    }else if(is(rgSet,"RGChannelSet")){
        ctrls <- getProbeInfo(rgSet,type="Control")}
      ctrls <- ctrls[ctrls$Address %in% rownames(rgSet),]
      ctrl_r <- assays(rgSet)$Red[ctrls$Address,]
      ctrl_g <- assays(rgSet)$Green[ctrls$Address,]
      CG.controls <- ctrls$Type %in% c("NORM_C", "NORM_G")
      AT.controls <- ctrls$Type %in% c("NORM_A", "NORM_T")

      cg_grn=ctrl_g[CG.controls,]
      at_red=ctrl_r[AT.controls,]
      cg_grn <- enmix(cg_grn,bgGI,bgParaEst,nCores)
      at_red <- enmix(at_red,bgRI,bgParaEst,nCores)

      Green.avg <- apply(cg_grn,2,huber_mu)
      Red.avg <- apply(at_red,2,huber_mu)
      ref <- mean(c(Red.avg,Green.avg))
      Grn.factor <- ref/Green.avg
      Red.factor <- ref/Red.avg
      }else if(dyeCorr =="RELIC"){
    if(is(rgSet,"rgDataSet")){
        ctrls <- getCGinfo(rgSet,type="ctrl")
    }else if(is(rgSet,"RGChannelSet")){
        ctrls <- getProbeInfo(rgSet,type="Control")}

    ctrls<-ctrls[ctrls$Address %in% rownames(rgSet),]
    ctrl_r<-assays(rgSet)$Red[ctrls$Address,]
    ctrl_g<-assays(rgSet)$Green[ctrls$Address,]
    CG.controls<-ctrls$Type %in% c("NORM_C","NORM_G")
    AT.controls<-ctrls$Type %in% c("NORM_A","NORM_T")
    cg_grn<-ctrl_g[CG.controls,];rownames(cg_grn)=
     ctrls$ExtendedType[CG.controls]
    at_red<-ctrl_r[AT.controls,];rownames(at_red)=
     ctrls$ExtendedType[AT.controls]
      cg_grn <- enmix(cg_grn,bgGI,bgParaEst,nCores)
      at_red <- enmix(at_red,bgRI,bgParaEst,nCores)
      }
    rm(rgSet)

    methData <- assays(mdat)$Meth
    N=ceiling(ncol(methData)/(nCores*10))
    parts=rep(1:N,each = ceiling(ncol(methData)/N))[1:ncol(methData)]
    for(i in 1:N){
    id=which(parts==i)
    methD=methData[,id]
    methD[probe_type == "IGrn",] <- enmix(methD[probe_type == "IGrn",],
     bgGI[id,], bgParaEst, nCores)
    methD[probe_type == "IRed",] <- enmix(methD[probe_type == "IRed",],
     bgRI[id,], bgParaEst, nCores)
    methD[probe_type == "II",] <- enmix(methD[probe_type == "II",],
     bgGII[id,], bgParaEst, nCores)
    methData[,id]=methD;
    }
    if (dyeCorr == "mean"){
      methData[probe_type == "IGrn",] <- sweep(methData[probe_type == "IGrn",],
       2, FUN="*", Grn.factor)
      methData[probe_type == "II",] <- sweep(methData[probe_type == "II",], 2,
       FUN="*", Grn.factor)
      methData[probe_type == "IRed",] <- sweep(methData[probe_type == "IRed",],
       2, FUN="*", Red.factor)
    }
    assays(mdat)$Meth <- methData
    rm(methData)

    unmethData <- assays(mdat)$Unmeth
    for(i in 1:N){
    id=which(parts==i)
    unmethD=unmethData[,id]
    unmethD[probe_type == "IGrn",] <- enmix(unmethD[probe_type == "IGrn",],
     bgGI[id,], bgParaEst, nCores)
    unmethD[probe_type == "IRed",] <- enmix(unmethD[probe_type == "IRed",],
     bgRI[id,], bgParaEst, nCores)
    unmethD[probe_type == "II",] <- enmix(unmethD[probe_type == "II",],
     bgRII[id,], bgParaEst, nCores)
    unmethData[,id]=unmethD;
    }
    if (dyeCorr == "mean"){
      unmethData[probe_type == "IGrn",] <- sweep(unmethData[probe_type ==
       "IGrn",], 2, FUN="*", Grn.factor)
      unmethData[probe_type == "IRed",] <- sweep(unmethData[probe_type ==
       "IRed",], 2, FUN="*", Red.factor)
      unmethData[probe_type == "II",] <- sweep(unmethData[probe_type ==
       "II",], 2, FUN="*", Red.factor)
    }
    assays(mdat)$Unmeth <- unmethData
    rm(unmethData)
    stopCluster(c1)
    if(dyeCorr =="RELIC"){mdat=relic(mdat,at_red,cg_grn)}

    annotation=c(paste("Backgroud_corr: ENmix,",bgParaEst,sep=""),
         paste("dyeBiasCorrection: ",dyeCorr,sep=""))
    if(is(mdat, "methDataSet")){
      metadata(mdat)$preprocessMethod=annotation
    }else if(is(mdat, "MethylSet")){
    mdat@preprocessMethod <- annotation
    }
    mdat
}

