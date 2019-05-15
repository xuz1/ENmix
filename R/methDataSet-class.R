# Exported classes -------------------------------------------------------------

setClass(
    "methDataSet",
    contains = "SummarizedExperiment"
)

# Exported functions -----------------------------------------------------------

methDataSet <- function(Meth = new("matrix"), Unmeth = new("matrix"),
    rowData=new("DataFrame"),...) {
    assays <- SimpleList(Meth = Meth, Unmeth = Unmeth)
    metadata <- list(preprocessMethod="")
    new("methDataSet",
        SummarizedExperiment(assays = assays, metadata=metadata,
        rowData=rowData,...)
    )
}

# Exported methods -------------------------------------------------------------

setMethod("show", signature(object = "methDataSet"), function(object) {
    callNextMethod()
    cat("PreprocessMethod: ", metadata(object)$preprocessMethod,"\n")
})


