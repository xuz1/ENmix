.guessArrayTypes <- function(nProbes) {
    if (nProbes >= 622000 && nProbes <= 623000) {
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylation450k",
            annotation = "ilmn12.hg19")
    } else if (nProbes >= 1032000 && nProbes <= 1053000) {
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylationEPIC",
            annotation = "ilm10b2.hg19")
    } else if (nProbes >= 54000 && nProbes <= 56000) {
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylation27k",
            annotation = "ilmn12.hg19")
    } else {
        arrayAnnotation <- c(array = "Unknown", annotation = "Unknown")
    }
    arrayAnnotation
}

readidat <- function(path = NULL,manifestfile=NULL,recursive = FALSE, verbose = FALSE) {
#    if(!require("illuminaio")){stop("Can not load in illuminaio package")}
    Grn.files <- list.files(path = path, pattern = "_Grn.idat$", recursive = recursive,
        ignore.case = FALSE, full.names = TRUE)
    Red.files <- list.files(path = path,pattern = "_Red.idat$",recursive = recursive,
        ignore.case = FALSE, full.names = TRUE)
    if (length(Grn.files) == 0 || length(Red.files) == 0) {
        stop("[readidat] No files with suffix _Grn.idat or _Red.idat were found")
    }else{message("[readidat] Found ",length(Grn.files)," files with suffix _Grn.idat")
          message("[readidat] Found ",length(Red.files)," files with suffix _Red.idat")}
    commonFiles <- intersect(
        x = sub("_Grn.idat$", "", Grn.files),
        y = sub("_Red.idat$", "", Red.files))
    if (length(commonFiles) == 0) {
        stop("[readidat] No IDAT files with both Red and Green channel were found")
    }
    basenames = commonFiles

    G.files <- paste(basenames, "_Grn.idat", sep = "")
    names(G.files) <- basename(basenames)
    R.files <- paste(basenames, "_Red.idat", sep = "")
    names(R.files) <- basename(basenames)

    G.Quants <- lapply(G.files, function(xx) {
        Quants <- readIDAT(xx)[["Quants"]]
        if (verbose) message("[readidat] Loaded ", basename(xx), ", ",nrow(Quants)," probes")
        Quants[,c("Mean","NBeads")]
    })

    allNProbes <- vapply(G.Quants, nrow, integer(1L))
    arrayTypes <- cbind(do.call(rbind, lapply(allNProbes, .guessArrayTypes)),
                        size = allNProbes)
    sameLength <- (length(unique(allNProbes)) == 1)
    sameArray <- (length(unique(arrayTypes[, "array"])) == 1)
    if (!sameLength){
        message(allNProbes)
        stop("[readidat] Array size are different.\n")
    }
    if(!sameArray){
        message(arrayTypes)
        stop("[readidat] Samples are from different array")
     }else{arraytype=as.character(arrayTypes[,"array"][1])
           annotation=as.character(arrayTypes[,"annotation"][1])
    }

    commonAddresses <- as.character(
        Reduce("intersect", lapply(G.Quants, rownames)))



    if (verbose) {message("[readidat] Loading array manifestfile...")}
    if(!is.null(manifestfile))
    {
      manifestfilepath=manifestfile 
      manifestfile=base::strsplit(manifestfilepath,split="/")[[1]]
      if(length(manifestfile)>1){
         manifestfile=manifestfile[length(manifestfile)]
         if (file.exists(manifestfile))file.remove(manifestfile)
         system(paste("wget ",manifestfilepath,sep=""))
      }
      manifest = readmanifest(manifestfile)
    }else{
      manifestpkg=paste(arraytype,"manifest",sep="")
#    manifest = getmanifest (manifestpkg)
      manifest = getmanifest (arraytype,annotation)
    }
    tmp=c(manifest$assay[,"Address"],manifest$ictrl[,"Address"])
    diff=setdiff(tmp,commonAddresses)
    if(length(diff)>5000){stop("IDAT file and manifest file are differed
      by more than 5000 probes")}else{
    commonAddresses=commonAddresses[commonAddresses %in% tmp]}

    if (verbose) {message("[readidat] Creating data matrices at green channel...")}
    GrnMean <- do.call(cbind,
        lapply(G.Quants, function(xx) xx[commonAddresses, "Mean"]))
    if (verbose) {message("[readidat] Creating data matrices for NBeads...")}
    NBeads <- do.call(cbind,
        lapply(G.Quants, function(xx) xx[commonAddresses, "NBeads"]))
    rm(G.Quants)
    
    R.Quants <- lapply(R.files, function(xx) {
        Quants <- readIDAT(xx)[["Quants"]]
        if (verbose) message("[readidat] Loaded ", basename(xx), ", ",nrow(Quants)," probes")
        Quants[,c("Mean")]
    })
    commonAddresses <- as.character(
        Reduce("intersect", lapply(R.Quants, names)))

    tmp=c(manifest$assay[,"Address"],manifest$ictrl[,"Address"])
    diff=setdiff(tmp,commonAddresses)
    if(length(diff)>5000){stop("IDAT file and manifest file are differed
      by more than 5000 probes")}else{
    commonAddresses=commonAddresses[commonAddresses %in% tmp]}

    if(!identical(rownames(GrnMean),commonAddresses)){
      commonAddresses=intersect(rownames(GrnMean),commonAddresses)
      GrnMean=GrnMean[commonAddresses,]
      NBeads=NBeads[commonAddresses,]
    }

    if (verbose) {message("[readidat] Creating data matrices at red channel...")}
    RedMean <- do.call(cbind,
        lapply(R.Quants, function(xx) xx[commonAddresses]))
    message("[readidat] Red Green Channel contains ",nrow(GrnMean)," probes and ",
        ncol(GrnMean)," Samples")

    tmp=manifest$ictrl
    tmp=tmp[tmp$Address %in% commonAddresses,]
    ictrl.anno=tmp
    tmp=manifest$assay[1:nrow(ictrl.anno),]
    tmp$Name=paste("ictrl",1:nrow(ictrl.anno),sep="")
    tmp$Address=ictrl.anno$Address
    tmp$Infinium_Design_Type="ctrl"
    tmp[,-which(names(tmp) %in% c("Name","Address","Infinium_Design_Type"))]=NA

    assay.anno=manifest$assay
    assay.anno=rbind(assay.anno,tmp)
    rownames(assay.anno)=assay.anno$Address
    assay.anno=assay.anno[rownames(GrnMean),]

    rownames(ictrl.anno)=ictrl.anno$Address

#check
    flag=!(assay.anno$Infinium_Design_Type %in% c("IA","IB","snpIA","snpIB") & 
                     !(assay.anno$Color_Channel %in% c("Red","Grn")))
    RedMean=RedMean[flag,]
    GrnMean=GrnMean[flag,]
    NBeads=NBeads[flag,]
    assay.anno=assay.anno[flag,]
    
    rgSet <- rgDataSet(
        Red = RedMean,
        Green = GrnMean,
        NBeads = NBeads,
	rowData=as(assay.anno,"DataFrame"),
	ictrl=as(ictrl.anno,"DataFrame")
    )
    metadata(rgSet)$Array=arraytype
    metadata(rgSet)$annotation=annotation
    rgSet
}


