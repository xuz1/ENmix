library(ENmix)
library(minfiData)
library(digest)

digestMatrix <- function(mat) {
    content <- sprintf("%.6f", mat)
    content[content == "-0.000000"] <- "0.000000"
    digest(c(content, rownames(mat), colnames(mat)))
}

#read in raw intensity data
sheet <- read.metharray.sheet(file.path(find.package("minfiData"),
        "extdata"), pattern = "csv$")
rgSet <- read.metharray.exp(targets = sheet,extended = TRUE)
qc<-QCinfo(rgSet)
mdat <- preprocessRaw(rgSet)
mdat.bg=preprocessENmix(rgSet,bgParaEst="oob",nCores=6)
mdat.filter=QCfilter(mdat,qcinfo=qc, samplethre = 0.05, CpGthre = 0.05)
mdat.quantile=normalize.quantile.450k(mdat,method="quantile1")

testDigests <- list(
    qc = list(detP = digestMatrix(qc$detP),
      nbead=digestMatrix(qc$nbead)),
   enmix = list(Meth = digestMatrix(getMeth(mdat.bg)),
      Unmeth = digestMatrix(getUnmeth(mdat.bg))),
    filter = list(Meth = digestMatrix(getMeth(mdat.filter)),
      Unmeth = digestMatrix(getUnmeth(mdat.filter))),
    quantile = list(Meth = digestMatrix(getMeth(mdat.quantile)),
      Unmeth = digestMatrix(getUnmeth(mdat.quantile)))
    )

save(testDigests, file = "../unitTests/testDigests.rda")

gc()
sessionInfo()

rm(list = ls())



