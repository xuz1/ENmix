four.identical<-function(x1,x2,x3,x4){
    return(identical(x1,x2) & identical(x1,x3) & identical(x1,x4))
}

oxBS.MLE<-function(beta.BS,beta.oxBS,N.BS,N.oxBS){
    if (!four.identical(rownames(beta.BS),rownames(beta.oxBS),rownames(N.BS),
    rownames(N.oxBS))){stop("Row names are not consistent!")}
    else if (!four.identical(colnames(beta.BS),colnames(beta.oxBS),
    colnames(N.BS),colnames(N.oxBS))){stop("Column names are not consistent!")
    }else {
    # filter the measurements to do estimation
    index<-!(is.na(beta.BS) | is.na(beta.oxBS) | is.na(N.BS) | is.na(N.oxBS)
       | N.BS==0 | N.oxBS==0)
    weight.BS<-N.BS[index]/(N.BS[index]+N.oxBS[index])
    weight.oxBS<-(1-weight.BS);

    # oxBS-MLE
    beta.5mC<-matrix(NA,nrow(beta.BS),ncol(beta.BS))
    beta.5mC[index]<-ifelse(beta.BS[index]>=beta.oxBS[index],beta.oxBS[index],
    weight.BS*beta.BS[index]+weight.oxBS*beta.oxBS[index])
    beta.5hmC<-matrix(NA,nrow(beta.BS),ncol(beta.BS))
    beta.5hmC[index]<-ifelse(beta.BS[index]>=beta.oxBS[index],beta.BS[index]
    -beta.oxBS[index],0)

    # name the rows and the columns
    rownames(beta.5mC)<-rownames(beta.BS)
    rownames(beta.5hmC)<-rownames(beta.BS)
    colnames(beta.5mC)<-colnames(beta.BS)
    colnames(beta.5hmC)<-colnames(beta.BS);

    # return the estimation
    return(list("5mC"=beta.5mC,"5hmC"=beta.5hmC))
    }
}
