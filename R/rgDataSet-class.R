# rgSet class

setClass("rgDataSet", contains = "SummarizedExperiment")

# Validity methods

setValidity2("rgDataSet", function(rgSet) {
    if(!all((metadata(rgSet)$ictrl)$Address %in% rownames(rgSet))){
    stop("Some control probes do not have intensity data, please provide correct manifest file")}
    probeIA=rowData(rgSet)$IlmnID[rowData(rgSet)$Infinium_Design_Type %in% c("IA","snpIA")]
    probeIB=rowData(rgSet)$IlmnID[rowData(rgSet)$Infinium_Design_Type %in% c("IB","snpIB")]
    if(!all(probeIA %in% probeIB) | !all(probeIA %in% probeIB)){
    stop("Some type I probes missed reads from channels A or B")}
    NULL
})

# rgSet functions
rgDataSet <- function(Red = new("matrix"), Green = new("matrix"),
    NBeads = new("matrix"),rowData=new("DataFrame"),ictrl= new("DataFrame"),...){
    assays <- SimpleList(Red = Red, Green = Green, NBeads = NBeads)
    metadata <- list(ictrl=ictrl)
    new("rgDataSet",
        SummarizedExperiment(assays = assays,metadata=metadata,rowData=rowData,...)
    )
}

# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "rgDataSet"), function(object) {
    callNextMethod()
    cat("Array: ", metadata(object)$Array, "\n")
    cat("Annotation: ", metadata(object)$annotation,"\n")
})


