nmode_est <-
function(x,minN=3,modedist=0.2)
{
    x <- x[!is.na(x)]
    if(sum(!is.na(x))>3){
    minN <- max(minN,ceiling(length(x)*0.05))
    y<-density(x)
    sn <- 3;nn <- sn*2+1
    n <- length(y$y)
    v <- matrix(NA,ncol=nn,nrow=n-nn+1)
    for(i in 1:nn){v[,i] <- y$y[i:(n-nn+i)]}
    ix <- sn+which(apply(v<v[,(sn+1)],1,sum)==(nn-1))
#number of samples under each peak must greater than minN
    if(length(ix)>1){
    valley <- array()
    for(i in 1:(length(ix)-1))
    {valley[i] <- y$x[ix[i]:ix[i+1]][which.min(y$y[ix[i]:ix[i+1]])]}
    v1 <- c(min(x),valley,max(x)+0.1)
    nn <- array();
    for(i in 1:(length(v1)-1)){nn[i] <- length(x[x>=v1[i] & x<v1[i+1]])}
    ix=ix[nn>=minN]
    nmode <- max(1,length(ix))
    }else{nmode=1}
    }else{nmode=1}
#nmode <- min(5,max(1,length(ix)))
    if(nmode>=2){
    flag <- 1;
    xx <- sort(y$x[ix][order(y$y[ix],decreasing=TRUE)[1:nmode]])
    }else{flag <- 0}
#distance between peaks must greater than modedist
    while(flag)
    {
    pdist <- xx[2:length(xx)]-xx[1:(length(xx)-1)]

    valley <- array()
    pidx <- which(y$x %in% xx)
    for(i in 1:(length(pidx)-1))
    {valley[i] <- y$x[pidx[i]:pidx[i+1]][which.min(y$y[pidx[i]:pidx[i+1]])]}
    v1 <- c(min(x),valley,max(x)+0.1)
    nn <- array();
    for(i in 1:(length(v1)-1)){nn[i] <- length(x[x>=v1[i] & x<v1[i+1]])}

    if(min(pdist)<modedist)
    {
    id_rm <- which.min(pdist)
    if(nn[id_rm]<nn[id_rm+1]){id_rm <- id_rm+1}
    xx <- xx[-id_rm]
    }else if(min(nn)<minN){
    id_rm <- which.min(nn)
    xx <- xx[-id_rm]
    }else{flag <- 0}
    nmode <- length(xx)
    if(nmode==1){flag <- 0}
    }
    nmode
}

nmode.mc <-function(x,minN=3,modedist=0.2,nCores=1)
{
    if(nCores>detectCores()){
    nCores <- detectCores();
    cat("Only ",detectCores()," Cores avialable, nCores was reset to ",
    detectCores(),"\n")
    }
    cat("Analysis is running, please wait...!","\n")

    CpG <- rownames(x)
    if(!(.Platform$OS.type == "unix")){
    c1 <- makeCluster(nCores)
    registerDoParallel(c1)
    }

    nmode=NULL
    N=ceiling(nrow(x)/(nCores*1200))
    parts=rep(1:N,each = ceiling(nrow(x)/N))[1:nrow(x)]
    for(i in 1:N){
    id=which(parts==i)
    x.1=x[id,]

    if(.Platform$OS.type == "unix"){
    nmode_est1 <-function(id,beta,minN,modedist)
    {nmode_est(beta[id,],minN=minN,modedist=modedist)}
    nmode.1=mclapply(1:nrow(x.1),nmode_est1,beta=x.1,minN=minN,modedist=modedist,
    mc.cores=nCores)
    nmode.1=unlist(nmode.1)
    }else{

    nmode.1 <- foreach(i = 1:nrow(x.1),.combine=c,.export=c("nmode_est")) %dopar% {
     i=i; nmode_est(x.1[i,], minN=minN,modedist=modedist)}
    }
    nmode=c(nmode,nmode.1)
    }
    if(!(.Platform$OS.type == "unix")){stopCluster(c1)}
    names(nmode) <- CpG
    nmode
}

