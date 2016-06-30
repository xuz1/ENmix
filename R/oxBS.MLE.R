oxBS.MLE <- function(beta.BS,beta.oxBS,N.BS,N.oxBS)
{
# filter the measurements to do estimation
    index<-!(is.na(beta.BS) | is.na(beta.oxBS) | is.na(N.BS) | is.na(N.oxBS) 
        | N.BS==0 | N.oxBS==0)
    weight.BS<-N.BS[index]/(N.BS[index]+N.oxBS[index]); 
    weight.oxBS<-(1-weight.BS)

# oxBS-MLE
    beta.5mC<-matrix(NA,nrow(beta.BS),ncol(beta.BS))
    beta.5mC[index]<-ifelse(beta.BS[index]>=beta.oxBS[index],beta.oxBS[index],
        weight.BS*beta.BS[index]+weight.oxBS*beta.oxBS[index])
    beta.5hmC<-matrix(NA,nrow(beta.BS),ncol(beta.BS))
    beta.5hmC[index]<-ifelse(beta.BS[index]>=beta.oxBS[index],beta.BS[index]
        -beta.oxBS[index],0)

# return the estimation
    return(list("5mC"=beta.5mC,"5hmC"=beta.5hmC))
}

