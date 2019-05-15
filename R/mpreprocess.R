mpreprocess <-function(rgSet,nCores=1,bgParaEst="oob",dyeCorr="RELIC", 
        qc=FALSE,qnorm=TRUE,qmethod="quantile1",
        foutlier=TRUE,rmcr=FALSE,impute=FALSE)
{
    if(!is(rgSet, "RGChannelSet") & !is(rgSet, "MethylSet") & !is(rgSet, "rgDataSet")
      & !is(rgSet, "methDataSet")){
    stop("ERROR: rgSet is not an object of class 'rgDataSet','methDataSet', 'RGChannelSet' or 'MethylSet'")}
    if(nCores>detectCores()){
    nCores <- detectCores();
    cat("Only ",nCores," cores are avialable in your computer,", 
       "argument nCores was reset to nCores=",nCores,"\n")
    }
    if(qc){
      if(is(rgSet, "rgDataSet") | is(rgSet, "RGChannelSetExtended")){qc <- QCinfo(rgSet)}else{
          qc=NULL; 
          cat("rgSet is not an object of class 'rgDataSet' or 'RGChannelSetExtended';\n",
              "QC will not be performed")}
        }else{qc=NULL}
    mdat <- preprocessENmix(rgSet, bgParaEst=bgParaEst, QCinfo=qc, 
        dyeCorr=dyeCorr,nCores=nCores)
    if(qnorm){mdat=norm.quantile(mdat,method=qmethod)}
    beta=rcp(mdat,qcscore=qc)
    if(foutlier){beta <- rm.outlier(beta,qcscore=qc,impute=impute,rmcr=rmcr)}
    beta
}

