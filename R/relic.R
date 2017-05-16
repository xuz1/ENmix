relic <- function(mdat,at_red=NULL,cg_grn=NULL)
{
    if(!is(mdat, "MethylSet"))
    {stop("The input must be an object of MethylSet\n")}
    if(is.null(at_red) | is.null(cg_grn))
    {stop("Internal control intensity matrics at_red and cg_grn must be
     provided\n")}

    # make sure the order is right
    name_grn <- gsub("G","A",gsub("C","T",rownames(cg_grn)))
    common=intersect(name_grn,rownames(at_red))
    if(length(common)<50)
       {stop("Please check Internal control intensity matrics at_red and
     cg_grn.\n Please assign internal control ExtendedType as row name 
     for matrics at_red and cg_grn\n")}
    at_red=at_red[common,];cg_grn=cg_grn[match(common,name_grn),]

    probe_type <- getProbeType(mdat, withColor=TRUE)
    m_grn=assays(mdat)$Meth[probe_type %in% c("IGrn","II"),]
    um_grn=assays(mdat)$Unmeth[probe_type %in% c("IGrn"),]
    for(i in 1:ncol(cg_grn))
    {
    index<-(at_red[,i]>10 & cg_grn[,i]>10)
    X<-cbind(rep(1,sum(index)),log(cg_grn[index,i])); Y<-log(at_red[index,i])
    temp<-solve(t(X)%*%X)%*%t(X)%*%Y;
    m_grn[,i]<-exp(log(m_grn[,i])*temp[2]+temp[1]);
    um_grn[,i]<-exp(log(um_grn[,i])*temp[2]+temp[1]);
    }
    assays(mdat)$Meth[probe_type %in% c("IGrn","II"),] <- m_grn
    assays(mdat)$Unmeth[probe_type %in% c("IGrn"),] <- um_grn
    mdat
}

