
#detection P values
calc_detP <- function(rgSet, detPtype = "negative") {
    if(!is(rgSet, "rgDataSet") & !is(rgSet, "RGChannelSetExtended"))
      stop("[QCinfo] The input should be an object of 'rgDataSet' or 'RGChannelSetExtended'")

    if(is(rgSet, "RGChannelSetExtended")){
    locusNames <- getManifestInfo(rgSet, "locusNames")
    detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                   dimnames = list(locusNames, colnames(rgSet)))
    if(detPtype == "negative"){
    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
    rBg <- assays(rgSet)$Red[controlIdx,,drop=FALSE]
    gBg <- assays(rgSet)$Green[controlIdx,,drop=FALSE]
    }else if(detPtype == "oob"){
       Igrn=getProbeInfo(rgSet,type="I-Green")
       Igrn=c(Igrn$AddressA,Igrn$AddressB)
       Ired=getProbeInfo(rgSet,type="I-Red")
       Ired=c(Ired$AddressA,Ired$AddressB)
    rBg <- assays(rgSet)$Red[Igrn,,drop=FALSE]
    gBg <- assays(rgSet)$Green[Ired,,drop=FALSE]
    }
    rMu <- matrixStats::colMedians(rBg)
    rSd <- matrixStats::colMads(rBg)
    gMu <- matrixStats::colMedians(gBg)
    gSd <- matrixStats::colMads(gBg)

    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")

    Ired_int=assays(rgSet)$Red[TypeI.Red$AddressA,] + assays(rgSet)$Red[TypeI.Red$AddressB,]
    Igrn_int=assays(rgSet)$Green[TypeI.Green$AddressA,] + assays(rgSet)$Green[TypeI.Green$AddressB,]
    II_int=assays(rgSet)$Red[TypeII$AddressA, ] + assays(rgSet)$Green[TypeII$AddressA, ]
    for (i in 1:ncol(rgSet)) {
        detP[TypeI.Red$Name, i] <- 1-pnorm(Ired_int[,i], mean=rMu[i]*2, sd=rSd[i]*sqrt(2))
        detP[TypeI.Green$Name, i] <- 1-pnorm(Igrn_int[,i], mean=gMu[i]*2, sd=gSd[i]*sqrt(2))
        detP[TypeII$Name, i] <- 1-pnorm(II_int[,i], mean=rMu[i]+gMu[i], sd=sqrt(rSd[i]^2+gSd[i]^2))
    }
}else if(is(rgSet, "rgDataSet")){
    if(detPtype == "negative"){
        ctrl=metadata(rgSet)$ictrl
        controlIdx <- ctrl$Address[ctrl$Type %in% "NEGATIVE"]
        rBg <- assays(rgSet)$Red[controlIdx,,drop=FALSE]
        gBg <- assays(rgSet)$Green[controlIdx,,drop=FALSE]
    }else if(detPtype=="oob"){
        anno=rowData(rgSet)
        rBg <- assays(rgSet)$Red[anno$Address[anno$Color_Channel %in% "Grn"],]
        gBg <- assays(rgSet)$Green[anno$Address[anno$Color_Channel %in% "Red"],]
    }

    rMu <- matrixStats::colMedians(rBg)
    rSd <- matrixStats::colMads(rBg)
    gMu <- matrixStats::colMedians(gBg)
    gSd <- matrixStats::colMads(gBg)

    cginfo=getCGinfo(rgSet)
    probeI=cginfo[cginfo$Infinium_Design_Type %in% c("I","snpI"),]
    probeIred=probeI[probeI$Color_Channel=="Red",]
    probeIgrn=probeI[probeI$Color_Channel=="Grn",]
    probeII=cginfo[cginfo$Infinium_Design_Type %in% c("II","snpII"),]

    Ired_int=assays(rgSet)$Red[probeIred$AddressA,]+assays(rgSet)$Red[probeIred$AddressB,]
    Igrn_int=assays(rgSet)$Green[probeIgrn$AddressA,]+assays(rgSet)$Green[probeIgrn$AddressB,]
    II_int=assays(rgSet)$Red[probeII$AddressA,]+assays(rgSet)$Green[probeII$AddressA,]
detP_Ired=matrix(NA_real_, ncol = ncol(Ired_int), nrow = nrow(Ired_int),dimnames =list(probeIred$Name,colnames(Ired_int)))
detP_Igrn=matrix(NA_real_, ncol = ncol(Igrn_int), nrow = nrow(Igrn_int),dimnames =list(probeIgrn$Name,colnames(Igrn_int)))
detP_II=matrix(NA_real_, ncol = ncol(II_int), nrow = nrow(II_int),dimnames =list(probeII$Name,colnames(II_int)))
for (i in 1:ncol(rgSet)) {
detP_Ired[,i]=1-pnorm(Ired_int[,i], mean=rMu[i]*2, sd=rSd[i]*sqrt(2))
detP_Igrn[,i]=1-pnorm(Igrn_int[,i], mean=gMu[i]*2, sd=gSd[i]*sqrt(2))
detP_II[,i]=1-pnorm(II_int[,i], mean=rMu[i]+gMu[i], sd=sqrt(rSd[i]^2+gSd[i]^2))
}
detP=rbind(detP_Ired,detP_Igrn,detP_II)
    }
    detP
}

