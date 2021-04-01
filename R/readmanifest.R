readmanifest <- function(file) {
    control.line <- try(system(sprintf("grep -n \\\\[Controls\\\\] %s", file), intern = TRUE),silent=TRUE)
    if(class(control.line)[1]=="try-error"){control.line <- system(sprintf("findstr -n \\[Controls\\] %s", file), intern = TRUE)}

    control.line <- as.integer(sub(":.*", "", control.line))
    stopifnot(length(control.line) == 1 && is.integer(control.line) && !is.na(control.line))

    assay.line <- try(system(sprintf("grep -n \\\\[Assay\\\\] %s", file), intern = TRUE),silent=TRUE)
    if(class(assay.line)[1]=="try-error"){assay.line <- system(sprintf("findstr -n \\[Assay\\] %s", file), intern = TRUE)}

    assay.line <- as.integer(sub(":.*", "", assay.line))
    stopifnot(length(assay.line) == 1 && is.integer(assay.line) && !is.na(assay.line))
    colNames <- readLines(file, n = assay.line + 1L)[assay.line + 1L]
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) <- colNames
    names(colClasses) <- make.names(names(colClasses))
    colClasses[c("MAPINFO")] <- "integer"

    manifest <- read.table(file, header = FALSE, col.names = names(colClasses),
                           sep = ",", comment.char = "", quote = "",
                           skip = assay.line + 1L, colClasses = colClasses,
                           nrows = control.line - assay.line - 2L)

    manifest$AddressA_ID <- gsub("^0*", "", manifest$AddressA_ID)
    manifest$AddressB_ID <- gsub("^0*", "", manifest$AddressB_ID)

if(sum(table(manifest$Name)>1)>0){manifest$Name=paste0(manifest$Name,"_",manifest$AddressA_ID)}

typeI=manifest[manifest$Infinium_Design_Type=="I",]
typeI$Infinium_Design_Type[grep("^rs",typeI$Name)]="snpI"
typeII=manifest[manifest$Infinium_Design_Type=="II",]
typeII$Infinium_Design_Type[grep("^rs",typeII$Name)]="snpII"
typeIA=typeI
typeIB=typeI
typeIA$Address=typeIA$AddressA_ID
typeIB$Address=typeIB$AddressB_ID
typeII$Address=typeII$AddressA_ID
typeIA$Infinium_Design_Type=paste(typeIA$Infinium_Design_Type,"A",sep="")
typeIB$Infinium_Design_Type=paste(typeIB$Infinium_Design_Type,"B",sep="")

var1=c("Name","Address","Infinium_Design_Type","Color_Channel","CHR","MAPINFO","Relation_to_UCSC_CpG_Island")
var2=names(manifest)[!(names(manifest) %in% c(var1,"AddressB_ID","AddressA_ID","IlmnID"))]
var12=c(var1,var2)
var12=var12[var12 %in% names(typeIA)]
typeIA=typeIA[,var12]
typeIB=typeIB[,var12]
typeII=typeII[,var12]

assay.anno=rbind(typeIA,typeIB,typeII)
names(assay.anno)[which(names(assay.anno)=="CHR")]="chr"
names(assay.anno)[which(names(assay.anno)=="MAPINFO")]="pos"
names(assay.anno)[which(names(assay.anno)=="Relation_to_UCSC_CpG_Island")]="Relation_to_Island"

    controls <- read.table(
        file = file,
        skip = control.line,
        sep = ",",
        comment.char = "",
        quote = "",
        colClasses = c(rep("character", 5)))[, 1:5]
    ictrl.anno <- controls[, 1:4]
    names(ictrl.anno) <- c("Address", "Type", "Color", "ExtendedType")

flag=(assay.anno$Address=="") | is.na(assay.anno$Address)
if(sum(flag)>0){assay.anno=assay.anno[!flag,]}
flag=(ictrl.anno$Address=="") | is.na(ictrl.anno$Address)
if(sum(flag)>0){ictrl.anno=ictrl.anno[!flag,]}

#check duplicate entries
assay.anno=assay.anno[!duplicated(assay.anno$Address),]
ictrl.anno=ictrl.anno[!duplicated(ictrl.anno$Address),]
assay.anno=assay.anno[!(assay.anno$Address %in% ictrl.anno$Address),]

return(list(assay=assay.anno,ictrl=ictrl.anno))
}

getmanifest<- function (arraytype,annotation)
{
#arraytype="IlluminaHumanMethylation450k"
#annotation="ilmn12.hg19"
manifestpkg=paste(arraytype,"manifest",sep="")
annopak=paste(arraytype,"anno.",annotation,sep="")


if(!require(manifestpkg, character.only = TRUE)){
    stop(sprintf("cannot load manifest package %s", manifestpkg))}
manifest=get(manifestpkg)

if(!require(annopak,character.only = TRUE)){
    stop(sprintf("cannot load annotation package %s", annopak))}
anno=get(annopak)
anno <- do.call(cbind, lapply(anno@defaults, get))

typeI=anno[anno$Type=="I",]
typeII=anno[anno$Type=="II",]
typeIA=typeI
typeIB=typeI
typeIA$Address=typeIA$AddressA
typeIB$Address=typeIB$AddressB
typeII$Address=typeII$AddressA

var=names(typeII)[!(names(typeII) %in% c("AddressB","AddressA","Color","Type","Name"))]
var=var[var %in% names(typeIA)]

typeIA=typeIA[,var]
typeIB=typeIB[,var]
typeII=typeII[,var]
anno=rbind(typeIA,typeIB,typeII)

#manifest
typeI=manifest@data$TypeI
typeIsnp=manifest@data$TypeSnpI
typeII=manifest@data$TypeII
typeIIsnp=manifest@data$TypeSnpII
ictrl.anno=manifest@data$TypeControl

typeI$Infinium_Design_Type="I"
typeIsnp$Infinium_Design_Type="snpI"
typeII$Infinium_Design_Type="II"
typeIIsnp$Infinium_Design_Type="snpII"
typeI=rbind(typeI,typeIsnp)
typeII=rbind(typeII,typeIIsnp)

typeIA=typeI[,c("Name","AddressA","Infinium_Design_Type","Color")]
typeIA$Infinium_Design_Type=paste(typeIA$Infinium_Design_Type,"A",sep="")
typeIB=typeI[,c("Name","AddressB","Infinium_Design_Type","Color")]
typeIB$Infinium_Design_Type=paste(typeIB$Infinium_Design_Type,"B",sep="")
typeII=typeII[,c("Name","AddressA","Infinium_Design_Type")]
typeII$Color_Channel=""
names(typeIA)=c("IlmnID","Address","Infinium_Design_Type","Color_Channel")
names(typeIB)=c("IlmnID","Address","Infinium_Design_Type","Color_Channel")
names(typeII)=c("IlmnID","Address","Infinium_Design_Type","Color_Channel")

assay.anno=rbind(typeIA,typeIB,typeII)
names(assay.anno)[which(names(assay.anno)=="IlmnID")]="Name"

assay.anno=merge(assay.anno,anno,by="Address",all.x=TRUE)

#order variables
var1=c("Name","Address","Infinium_Design_Type","Color_Channel","chr","pos","Relation_to_Island")
var2=names(assay.anno)[!(names(assay.anno) %in% c(var1))]
assay.anno=assay.anno[,c(var1,var2)]

#check duplicate entries
assay.anno=assay.anno[!duplicated(assay.anno$Address),]
ictrl.anno=ictrl.anno[!duplicated(ictrl.anno$Address),]
assay.anno=assay.anno[!(assay.anno$Address %in% ictrl.anno$Address),]

return(list(assay=assay.anno,ictrl=ictrl.anno))
}

