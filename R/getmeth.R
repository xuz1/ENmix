#create methDataSet from rgDataSet
getmeth <-function(rgSet){
  if(!is(rgSet, "rgDataSet") & !is(rgSet, "RGChannelSet")){
    stop("[getmeth] The input must be a rgDataSet or RGChannelSet")}

if(is(rgSet, "RGChannelSet")){return(preprocessRaw(rgSet))}else{

cginfo=getCGinfo(rgSet)

probeI=cginfo[cginfo$Infinium_Design_Type %in% c("I","snpI"),]
probeIred=probeI[probeI$Color_Channel=="Red",]
probeIgrn=probeI[probeI$Color_Channel=="Grn",]
probeII=cginfo[cginfo$Infinium_Design_Type %in% c("II","snpII"),]

MIRed <- assays(rgSet)$Red[probeIred$AddressB,]
UIRed <- assays(rgSet)$Red[probeIred$AddressA,]
MIGrn <- assays(rgSet)$Green[probeIgrn$AddressB,]
UIGrn <- assays(rgSet)$Green[probeIgrn$AddressA,]

dimnames(MIRed)[[1]]=probeIred$Name
dimnames(UIRed)[[1]]=probeIred$Name
dimnames(MIGrn)[[1]]=probeIgrn$Name
dimnames(UIGrn)[[1]]=probeIgrn$Name

MII <- assays(rgSet)$Green[probeII$AddressA,]
UII <- assays(rgSet)$Red[probeII$AddressA,]

dimnames(MII)[[1]]=probeII$Name
dimnames(UII)[[1]]=probeII$Name

M=rbind(MIRed,MIGrn,MII)
U=rbind(UIRed,UIGrn,UII)

rowData=rbind(probeI,probeII)
dimnames(rowData)[[1]]=rowData$Name
rowData=rowData[dimnames(M)[[1]],]

    out <- methDataSet(
        Meth=M,
        Unmeth = U,
        rowData=as(rowData,"DataFrame")
    )
    out
}}


