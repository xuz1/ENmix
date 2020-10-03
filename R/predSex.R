predSex <-function(mdat,cutoff=2){
  if(!is(mdat, "rgDataSet") & !is(mdat, "methDataSet")){
    stop("[predSex] The input must be a rgDataSet or methDataSet")}

  if(is(mdat, "rgDataSet")){mdat=getmeth(mdat)}
  xmdat=mdat[rowData(mdat)$chr %in% "chrX",]
  ymdat=mdat[rowData(mdat)$chr %in% "chrY",]
  xCN <- log2(assays(xmdat)$Meth + assays(xmdat)$Unmeth)
  yCN <- log2(assays(ymdat)$Meth + assays(ymdat)$Unmeth)
  xMed=apply(xCN,2,median,na.rm=TRUE)    
  yMed=apply(yCN,2,median,na.rm=TRUE)
  diff=xMed -yMed
  sex <- ifelse(diff > cutoff, "F", "M")
  return(data.frame(id=colnames(mdat),sex=sex))

#add in prediction based on beta values
#kmeans clustering, in minfi
}



