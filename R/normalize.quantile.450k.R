normalize.q <- function(x)
{
    if(sum(is.na(x))>0)
    {stop("Multi-array quantile normalization does not allow missing value")}
    rname <- rownames(x)
    cname <- colnames(x)
    x <- normalize.quantiles(x)
    rownames(x) <- rname
    colnames(x) <- cname
    x
}


normalize.quantile.450k <- function(mdat,method="quantile1")
{
    if(!is(mdat, "MethylSet")){stop("object needs to be of class 'MethylSet'")}
    cat("Analysis is running, please wait...!","\n")

    if (method == "quantile1")
    {
        anno <- getAnnotation(mdat)
        mdat_I <- mdat[anno$Type == "I",]
        mdat_II <- mdat[anno$Type == "II",]
        assayDataElement(mdat_I, "Meth")<-normalize.q(assayData(mdat_I)$Meth)
        assayDataElement(mdat_I, "Unmeth")<-normalize.q(assayData(mdat_I)$Unmeth)
        assayDataElement(mdat_II, "Meth")<-normalize.q(assayData(mdat_II)$Meth)
        assayDataElement(mdat_II, "Unmeth")<-normalize.q(assayData(mdat_II)$Unmeth)
        methData <- rbind(assayData(mdat_I)$Meth,assayData(mdat_II)$Meth)
        unmethData <- rbind(assayData(mdat_I)$Unmeth,assayData(mdat_II)$Unmeth)
        CpGID <- rownames(mdat)
        SampleID <- colnames(mdat)
        methData <- methData[CpGID,]
        unmethData <- unmethData[CpGID,]
        assayDataElement(mdat, "Meth")<-methData
        assayDataElement(mdat, "Unmeth")<-unmethData
    }
    else if(method == "quantile2")
    {
        anno <- getAnnotation(mdat)
        mdat_I <- mdat[anno$Type == "I",]
        mdat_II <- mdat[anno$Type == "II",]
        mat <- rbind(assayData(mdat_I)$Meth,assayData(mdat_I)$Unmeth)
        mat  <-  normalize.q(mat)
        assayDataElement(mdat_I, "Meth") <- mat[1:(nrow(mat)/2),]
        assayDataElement(mdat_I, "Unmeth") <- mat[((nrow(mat)/2)+1):nrow(mat),]
        mat <- rbind(assayData(mdat_II)$Meth,assayData(mdat_II)$Unmeth)
        mat  <-  normalize.q(mat)
        assayDataElement(mdat_II, "Meth") <- mat[1:(nrow(mat)/2),]
        assayDataElement(mdat_II, "Unmeth") <- mat[((nrow(mat)/2)+1):nrow(mat),]
        methData <- rbind(assayData(mdat_I)$Meth,assayData(mdat_II)$Meth)
        unmethData <- rbind(assayData(mdat_I)$Unmeth,assayData(mdat_II)$Unmeth)
        CpGID <- rownames(mdat)
        SampleID <- colnames(mdat)
        methData <- methData[CpGID,]
        unmethData <- unmethData[CpGID,]
        assayDataElement(mdat, "Meth")<-methData
        assayDataElement(mdat, "Unmeth")<-unmethData

    }
    else if(method == "quantile3")
    {
        ##quantile normalization of combined intensity values
        ##similar with lumi lumimethyN
        mat <- rbind(assayData(mdat)$Meth,assayData(mdat)$Unmeth)
        mat  <-  normalize.q(mat)
        methData <- mat[1:(nrow(mat)/2),]
        unmethData <- mat[((nrow(mat)/2)+1):nrow(mat),]
        assayDataElement(mdat, "Meth")<-methData
        assayDataElement(mdat, "Unmeth")<-unmethData
    }
    mdat@preprocessMethod <- c(mu.norm="normalize.quantile", preprocessMethod(mdat))
    mdat
}
