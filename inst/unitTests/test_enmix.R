test_enmix <- function() {
    stopifnot(require(minfiData))
    stopifnot(require(digest))
    load(file.path(path.package("ENmix"), "unitTests", "testDigests.rda"))

    digestMatrix <- function(mat) {
        content <- sprintf("%.6f", mat)
        content[content == "-0.000000"] <- "0.000000"
        digest(c(content, rownames(mat), colnames(mat)))
    }

    ##read in raw intensity data
    sheet <- read.450k.sheet(file.path(find.package("minfiData"),
        "extdata"), pattern = "csv$")
##to save time for package build
#    rgSet <- read.450k.exp(targets = sheet,extended = TRUE)
#    qc<-QCinfo(rgSet)
#    mdat <- preprocessRaw(rgSet)
#    mdat.bg=bgcorrect(rgSet,bgParaEst="oob",nCores=6)
#    mdat.filter=QCfilter(mdat,qcinfo=qc, samplethre = 0.01, CpGthre = 0.05)
#    mdat.quantile=normalization(mdat,method="quantile1")

#    checkEquals(testDigests$qc$detP, digestMatrix(qc$detP))
#    checkEquals(testDigests$qc$nbead, digestMatrix(qc$nbead))
#    checkEquals(testDigests$enmix$Meth, digestMatrix(getMeth(mdat.bg)))
#    checkEquals(testDigests$enmix$Unmeth, digestMatrix(getUnmeth(mdat.bg)))
#    checkEquals(testDigests$filter$Meth, digestMatrix(getMeth(mdat.filter)))
#    checkEquals(testDigests$filter$Unmeth, digestMatrix(getUnmeth(mdat.filter)))
#    checkEquals(testDigests$quantile$Meth, digestMatrix(getMeth(mdat.quantile)))
#    checkEquals(testDigests$quantile$Unmeth, digestMatrix(getUnmeth(mdat.quantile)))
}

