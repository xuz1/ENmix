
methscore<-function(datMeth,datPheno=NULL,fastImputation=FALSE,normalize=TRUE,GrimAgeComponent=NULL,UserRef=NULL,ForceUserRef=FALSE)
{
#input datMeth value matrix must be in between 0 to 1, percentage of cg methylated
if(min(datMeth,na.rm=TRUE)<0 | max(datMeth,na.rm=TRUE)>1){
    stop("Warning: Methylation datMeth value should be within [0,1]")}

if(ForceUserRef & is.null(UserRef)){ForceUserRef=FALSE; message("Warning: UserRef is not provided, set ForceUserRef=FALSE")}
if(!is.null(UserRef)){if(!all(c("cg","meth_mean") %in% colnames(UserRef)))stop("UserRef must include variables cg and meth_mean")}

#if GrimAge Components were provided
if(!is.null(GrimAgeComponent)){
	component= c("DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS", "DNAmPAI1", "DNAmTIMP1","DNAmGrimAge")
	component=component[!component %in% colnames(GrimAgeComponent)]
	if(length(component)>0){
		stop(paste0("Variables ",component," were not found in GrimAgeComponent"))
	}

	if(sum(colnames(datMeth) %in% as.vector(GrimAgeComponent$SampleID))<ncol(datMeth)){
		stop("Some samples in datMeth do not have data in GrimAgeComponent")
	}else{rownames(GrimAgeComponent)=as.vector(GrimAgeComponent$SampleID)}
}

#remove CpGs with missing values for all samples
datMeth=datMeth[!rowSums(is.na(datMeth))==ncol(datMeth),]

#if datMeth is a vector, convert it to a matrix
if(is.vector(datMeth)){datMeth=as.matrix(datMeth)}

if(is.null(datPheno)){ stop("Error: datPheno is missing.")} else {
  datPheno=as.data.frame(datPheno)
  if(!("Age" %in% colnames(datPheno))){
    stop("Error: datPheno must have a column named Age")
  }
  if(!("Female" %in% colnames(datPheno))){
    stop("Error: datPheno must have a column named Female")
  }
  if(!("SampleID" %in% colnames(datPheno))){
    stop("Error: datPheno must have a column named SampleID")
  }else if(sum(colnames(datMeth) %in% datPheno$SampleID)<ncol(datMeth)){
     stop("Error: some samples in datMeth do not have phenotype in datPheno")
  }else{
     rownames(datPheno)=datPheno$SampleID
     datPheno=datPheno[colnames(datMeth),]
  }
  if(sum(is.na(datPheno))>0){stop("Error: missing value in phenotype data is not allowed")}
}

#define a few global variables
refmeth=mPOA_Models=episcore_model=pcc=refmeth2=DNAmFitnessModels=NULL
cAge_CpG=bAge_CpG=EpiToc_CpGs=MiAgedat=methscore_dict=NULL
#load in all model and reference data
load(system.file("mage_ref.RData",package="ENmix"))
load(system.file("mPOA_Models.RData",package="ENmix"))
#load(system.file("pcclock_model.RData",package="ENmix"))
pcc1=pcc2=pcc3=pcc4=NULL
load(system.file("pcclock_model1.RData",package="ENmix"))
load(system.file("pcclock_model2.RData",package="ENmix"))
load(system.file("pcclock_model3.RData",package="ENmix"))
load(system.file("pcclock_model4.RData",package="ENmix"))
pcc=c(pcc1,pcc2,pcc3,pcc4)

#check whether datMeth probe names include suffix, such as in EPICv2 array
if(sum(as.vector(refmeth$cg) %in% rownames(datMeth))<50){datMeth=rm.cgsuffix(datMeth)}

#combine reference CpG
#normalize to Horvath refmeth before combine
refcg=as.vector(refmeth$cg)
refcg=unique(c(refcg,as.vector(mPOA_Models$gold_standard_probes)))
refcg=unique(c(refcg,as.vector(episcore_model$CpG_Site)))
refcg=unique(c(refcg,as.vector(frailty_model$cg[-1])))
refcg=unique(c(refcg,names(pcc$pcCpGs)))
refcg=unique(c(refcg,as.vector(refmeth2$cg)))
refcg=unique(c(refcg,as.vector(cAge_CpG$means$cpg)))
refcg=unique(c(refcg,as.vector(bAge_CpG$cpgs$CpG_Site)))
refcg=unique(c(refcg,as.vector(DNAmFitnessModels$AllCpGs)))

#remove unnecessary probes for estimation
datMeth=datMeth[rownames(datMeth) %in% refcg,]

#imputation to fill missing values
noMissingPerCpG=apply(is.na(datMeth),1,sum)
if(max(noMissingPerCpG)>0 & ncol(as.matrix(datMeth))>1){
message("Imputation to fill missing values...")
if(!fastImputation){
    set.seed(1)
    datMeth=t(impute::impute.knn(t(datMeth))$data)
}else{
    for (i in which(noMissingPerCpG>0)){
        idx=is.na(datMeth[i,])
        datMeth[i,][is.na(datMeth[i,])] = mean(as.numeric(datMeth[i,]),na.rm=T)
}}}


#Count of model required CpGs that are missing in user provided data
missedcg=NULL
cgcheck=data.frame(predictor=NA,nCpG_required=NA,nCpG_present=NA,nCpG_missing=NA)
np=1
cgcheck[np,1]="HorvathAge";cgcheck[np,2]=nrow(horvath)-1;cgcheck[np,3]=sum(as.vector(horvath$cg) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="PhenoAge";cgcheck[np,2]=nrow(phenoage)-1;cgcheck[np,3]=sum(as.vector(phenoage$cg) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="HannumAge";cgcheck[np,2]=nrow(hannum);cgcheck[np,3]=sum(as.vector(hannum$cg) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="PACE";cgcheck[np,2]=length(mPOA_Models$model_probes);cgcheck[np,3]=sum(as.vector(mPOA_Models$model_probes) %in% rownames(datMeth))
#frailty
np=np+1
cgcheck[np,1]="Frailty";cgcheck[np,2]=nrow(frailty_model)-1;cgcheck[np,3]=sum(as.vector(frailty_model$cg[-1]) %in% rownames(datMeth))
#other,PEDBE,EpiToc,EpiToc2,Zhang10CpG,Horvath2,MiAge,DNAmTL,PEDBE,GACPC,GARPC,GARRPC,Bohlin,Knight
np=np+1
cgcheck[np,1]="DNAmTL";cgcheck[np,2]=nrow(DNAmTL_CpGs);cgcheck[np,3]=sum(as.vector(DNAmTL_CpGs$ID) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="pcgtAge";cgcheck[np,2]=length(EpiToc_CpGs);cgcheck[np,3]=sum(EpiToc_CpGs %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="TNSC";cgcheck[np,2]=nrow(EpiToc2_CpGs);cgcheck[np,3]=sum(rownames(EpiToc2_CpGs) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="TNSC2";cgcheck[np,2]=nrow(EpiToc2_CpGs);cgcheck[np,3]=sum(rownames(EpiToc2_CpGs) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="Zhang10CpG";cgcheck[np,2]=nrow(Zhang_10_CpG);cgcheck[np,3]=sum(as.vector(Zhang_10_CpG$Marker) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="Horvath2";cgcheck[np,2]=nrow(Horvath2_CpGs);cgcheck[np,3]=sum(as.vector(Horvath2_CpGs$ID) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="MiAge";cgcheck[np,2]=nrow(MiAgedat$MiAge_CpGs);cgcheck[np,3]=sum(as.vector(MiAgedat$MiAge_CpGs$CpGs) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="PedBE";cgcheck[np,2]=nrow(PEDBE_CpGs);cgcheck[np,3]=sum(as.vector(PEDBE_CpGs$ID) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="GACPC";cgcheck[np,2]=nrow(GACPC_CpGs);cgcheck[np,3]=sum(as.vector(GACPC_CpGs$CpG) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="GARPC";cgcheck[np,2]=nrow(GARPC_CpGs);cgcheck[np,3]=sum(as.vector(GARPC_CpGs$CpG) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="GARRPC";cgcheck[np,2]=nrow(GARRPC_CpGs);cgcheck[np,3]=sum(as.vector(GARRPC_CpGs$CpG) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="BohlinGAge";cgcheck[np,2]=nrow(Bohlin_CpGs);cgcheck[np,3]=sum(as.vector(Bohlin_CpGs$CpG) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="KnightGAge";cgcheck[np,2]=nrow(Knight_CpGs);cgcheck[np,3]=sum(as.vector(Knight_CpGs$CpG) %in% rownames(datMeth))
#pcclock
pcvar1 = c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", 
                     "PCDNAmTL", "PCPACKYRS", "PCADM", "PCB2M", "PCCystatinC",
                     "PCGDF15", "PCLeptin", "PCPAI1", "PCTIMP1", "PCGrimAge")
pcvar=c(pcvar1)
for(vv in pcvar){
np=np+1
cgcheck[np,1]=vv;cgcheck[np,2]=length(pcc$pcCpGs);cgcheck[np,3]=sum(names(pcc$pcCpGs) %in% rownames(datMeth))
}
#cAge
np=np+1
cgcheck[np,1]="cAge";cgcheck[np,2]=nrow(cAge_CpG$means);cgcheck[np,3]=sum(as.vector(cAge_CpG$means$cpg) %in% rownames(datMeth))
#bAge
np=np+1
cgcheck[np,1]="bAge";cgcheck[np,2]=nrow(bAge_CpG$cpgs);cgcheck[np,3]=sum(as.vector(bAge_CpG$cpgs$CpG_Site) %in% rownames(datMeth))
np=np+1
cgcheck[np,1]="bAge_Years";cgcheck[np,2]=nrow(bAge_CpG$cpgs);cgcheck[np,3]=sum(as.vector(bAge_CpG$cpgs$CpG_Site) %in% rownames(datMeth))
#episcore
for(ms in unique(episcore_model$Predictor)){
    np=np+1
    tmp=episcore_model[episcore_model$Predictor == ms,]
    cgcheck[np,1]=ms
    cgcheck[np,2]=nrow(tmp)
    cgcheck[np,3]=sum(as.vector(tmp$CpG_Site) %in% rownames(datMeth))
}
#DNAm Fit Age
fa=c("DNAmGait_noAge","DNAmGrip_noAge","DNAmVO2max","DNAmGait_wAge","DNAmGrip_wAge","DNAmFEV1_wAge","DNAmFitAge")
for(i in 1:length(fa)){
  np=np+1
  cgcheck[np,1]=fa[i];cgcheck[np,2]=length(DNAmFitnessModels$AllCpGs);cgcheck[np,3]=sum(as.vector(DNAmFitnessModels$AllCpGs) %in% rownames(datMeth))
}

message("Number of CpGs missing for estimates\n")
cgcheck$nCpG_missing=cgcheck$nCpG_required-cgcheck$nCpG_present
print(cgcheck)
if(max(cgcheck$nCpG_missing)>0){message("Missing probes will be imputed with reference values\n")}
message("See file summary_methscore_CpG.csv for the CpG count summary and citation for each predictor")
cgcheck=merge(methscore_dict,cgcheck,by="predictor")
rownames(cgcheck)=cgcheck$predictor
cgcheck=cgcheck[as.vector(methscore_dict$predictor),]
write.csv(cgcheck,file="summary_methscore_CpG.csv",row.names=FALSE)

#function for normalization using a modified RCP method, and imputing missing CpGs
norm_imputing <-function(dat,refdat,normalize,missedcg){
   if(normalize){dat=rcp2(dat,refdat)}
    #replace missing probes with reference values
    if(length(missedcg)>0){
        refdat=refdat[as.vector(refdat$cg) %in% missedcg,]
        missedcg.dat=matrix(rep(refdat$meth_mean,ncol(datMeth)),nrow=length(missedcg))
        rownames(missedcg.dat)=as.vector(refdat$cg)
        colnames(missedcg.dat)=colnames(dat)
        dat=rbind(dat,missedcg.dat)
     }
     return(dat)    
}

#Age transformation and probe annotation functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }


###Methylation predictors
#########################################
#methylation ages
message("Calculating methylation scores ...")
mScore=data.frame(SampleID=colnames(datMeth))

modelcg=unique(c(as.vector(horvath$cg[-1]),as.vector(phenoage$cg[-1]),as.vector(hannum$cg)))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=refmeth;names(refdat)=c("cg","meth_mean")
#If user provied reference data include all modelcg, then use UserRef
if(!is.null(UserRef) & ForceUserRef){
    tmp=modelcg[!modelcg %in% as.vector(UserRef$cg)]
    if(length(tmp)==0){refdat=UserRef}else{
            message(paste0("The following model required CpGs are missing in UserRef:"))
            message(paste0(tmp,collaps=" "))
            message("System reference were used for HorvathAge,mAge_Hannum,PhenoAge")
    }
}
datMeth2=norm_imputing(datMeth,refdat=refdat,normalize,missedcg)

#HorvathAge
intercept=horvath$coef[horvath$cg=="(Intercept)"]
horvath=horvath[as.vector(horvath$cg) %in% rownames(datMeth2),]
mAge=anti.trafo(colSums(as.numeric(horvath$coef) * datMeth2[as.vector(horvath$cg),],na.rm=T)+intercept)
mScore$HorvathAge=mAge[as.vector(mScore$SampleID)]
#Hannum Age
hannum=hannum[hannum$cg %in% rownames(datMeth2),]
mAge=colSums(hannum$coef * datMeth2[as.vector(hannum$cg),],na.rm=T)
mScore$HannumAge=mAge[as.vector(mScore$SampleID)]
#PhenoAge
intercept=phenoage$coef[phenoage$cg=="intercept"]
phenoage=phenoage[phenoage$cg %in% rownames(datMeth2),]
mAge=colSums(phenoage$coef * datMeth2[as.vector(phenoage$cg),],na.rm=T)+intercept
mScore$PhenoAge=mAge[as.vector(mScore$SampleID)]

#########################################
#PACE estimation
modelcg=unique(as.vector(mPOA_Models$model_probes))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=data.frame(cg=mPOA_Models$gold_standard_probes,meth_mean=mPOA_Models$gold_standard_means)
#If user provied reference data include all modelcg, then use UserRef
if(!is.null(UserRef) & ForceUserRef){
    tmp=modelcg[!modelcg %in% as.vector(UserRef$cg)]
    if(length(tmp)==0){refdat=UserRef}else{
            message(paste0("The following model required CpGs are missing in UserRef:"))
            message(paste0(tmp,collaps=" "))
            message("System reference were used for PACE")
    }
}
datMeth2=norm_imputing(datMeth,refdat=refdat,normalize,missedcg)

PACE=DunedinPACE(betas=datMeth2, proportionOfProbesRequired = 0.8)
mScore$PACE=PACE[as.vector(mScore$SampleID)]

#########################################
#Frailty
modelcg=unique(as.vector(frailty_model$cg[-1]))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=frailty_model[-1,c("cg","meth_mean")]
#If user provied reference data include all modelcg, then use UserRef
if(!is.null(UserRef) & ForceUserRef){
    tmp=modelcg[!modelcg %in% as.vector(UserRef$cg)]
    if(length(tmp)==0){refdat=UserRef}else{
            message(paste0("The following model required CpGs are missing in UserRef:"))
            message(paste0(tmp,collaps=" "))
            message("System reference were used for eFRS")
    }
}
datMeth2=norm_imputing(datMeth,refdat=refdat,normalize,missedcg)

intercept=frailty_model$coef[frailty_model$cg=="intercept"]
frailty_model=frailty_model[frailty_model$cg %in% rownames(datMeth2),]
eFRS=colSums(frailty_model$coef * datMeth2[as.vector(frailty_model$cg),],na.rm=T)+intercept
mScore$Frailty=eFRS[as.vector(mScore$SampleID)]


#########################################
#other,PEDBE,EpiToc,EpiToc2,Zhang10CpG,Horvath2,MiAge,DNAmTL,PEDBE,GACPC,GARPC,GARRPC,Bohlin,Knight
modelcg=unique(as.vector(PEDBE_CpGs$ID))
modelcg=unique(c(modelcg,as.vector(EpiToc_CpGs)))
modelcg=unique(c(modelcg,rownames(EpiToc2_CpGs)))
modelcg=unique(c(modelcg,as.vector(Zhang_10_CpG$Marker)))
modelcg=unique(c(modelcg,as.vector(Horvath2_CpGs$ID)))
modelcg=unique(c(modelcg,as.vector(MiAgedat$MiAge_CpGs$CpGs)))
modelcg=unique(c(modelcg,as.vector(DNAmTL_CpGs$ID)))
modelcg=unique(c(modelcg,as.vector(PEDBE_CpGs$ID)))
modelcg=unique(c(modelcg,as.vector(GACPC_CpGs$CpG)))
modelcg=unique(c(modelcg,as.vector(GARPC_CpGs$CpG)))
modelcg=unique(c(modelcg,as.vector(GARRPC_CpGs$CpG)))
modelcg=unique(c(modelcg,as.vector(Bohlin_CpGs$CpG)))
modelcg=unique(c(modelcg,as.vector(Knight_CpGs$CpG)))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=refmeth2
#If user provied reference data include all modelcg, then use UserRef
if(!is.null(UserRef)){
    tmp=modelcg[!modelcg %in% as.vector(UserRef$cg)]
    if(length(tmp)==0){refdat=UserRef}else{
	    message(paste0("The following model required CpGs are missing in UserRef:"))
            message(paste0(tmp,collaps=" "))
	    message("System reference were used for PEDBE,EpiToc,EpiToc2,Zhang10CpG,Horvath2,MiAge,DNAmTL,PEDBE,GACPC,GARPC,GARRPC,Bohlin,Knight")
    }
}

names(refdat)=c("cg","meth_mean")
datMeth2=norm_imputing(datMeth,refdat=refdat,normalize,missedcg)

#DNAmTL
DNAmTL_CpGs=DNAmTL_CpGs[DNAmTL_CpGs$ID %in% rownames(datMeth2),]
DNAmTL=colSums(DNAmTL_CpGs$Coef * datMeth2[as.vector(DNAmTL_CpGs$ID),],na.rm=T)-7.924780053
mScore$DNAmTL=DNAmTL[as.vector(mScore$SampleID)]

#EpiTOC
pcgtAge <- colMeans(datMeth2[rownames(datMeth2) %in% EpiToc_CpGs,],na.rm=TRUE)
mScore$pcgtAge=pcgtAge[as.vector(mScore$SampleID)]

#EpiTOC2
EpiToc2_CpGs=EpiToc2_CpGs[rownames(EpiToc2_CpGs) %in% rownames(datMeth2),]
tmp.m=datMeth2[rownames(EpiToc2_CpGs),]
TNSC <- 2*colMeans((1/(EpiToc2_CpGs[,1]*(1-EpiToc2_CpGs[,2]))) * (tmp.m - EpiToc2_CpGs[,2]),na.rm=TRUE)
#approximated = T
TNSC2 <- 2*colMeans((1/EpiToc2_CpGs[,1]) * tmp.m,na.rm=TRUE)
mScore$TNSC=TNSC[as.vector(mScore$SampleID)]
mScore$TNSC2=TNSC2[as.vector(mScore$SampleID)]

#Zhang 10 CpG clock
Zhang_10_CpG=Zhang_10_CpG[as.vector(Zhang_10_CpG$Marker) %in% rownames(datMeth2),]
Zhang10CpG=colSums(Zhang_10_CpG$coef * datMeth2[as.vector(Zhang_10_CpG$Marker),],na.rm=T)
mScore$Zhang10CpG=Zhang10CpG[as.vector(mScore$SampleID)]

#Horvath2
Horvath2_CpGs=Horvath2_CpGs[as.vector(Horvath2_CpGs$ID) %in% rownames(datMeth2),]
Horvath2=anti.trafo(colSums(Horvath2_CpGs$Coef * datMeth2[as.vector(Horvath2_CpGs$ID),],na.rm=T)-0.447119319)
mScore$Horvath2=Horvath2[as.vector(mScore$SampleID)]

#MiAge mitotic age
#function used by the original MiAge code found <http://www.columbia.edu/~sw2206/softwares.htm>
MiAge_fr2 <- function (x, b, c, d, betaj)
{
    nj = x
    return(sum((c + b^(nj - 1) * d - betaj)^2, na.rm = T))
}
MiAge_grr2 <-function (x, b, c, d, betaj)
{
    nj = x
    return(2 * sum((c + b^(nj - 1) * d - betaj) * b^(nj - 1) *
        log(b) * d, na.rm = T))
}
mitotic.age <- function (beta, b, c, d)
{
#    library(methods)
    upperage = 10000
    lowerage = 10
    n = rep(500, ncol(beta));names(n)=colnames(beta)
    no.initial.n = 5
    for (j in 1:ncol(beta)) {
        current.value = MiAge_fr2(n[j], b, c, d, beta[, j])
        columnoptim = vector("list", no.initial.n)
        val = rep(NA, no.initial.n)
        for (jj in 1:(no.initial.n - 1)) {
            temp = try(optim(par = lowerage + jj * (upperage -
                lowerage)/no.initial.n, fn = MiAge_fr2, gr = MiAge_grr2,
                b = b, c = c, d = d, betaj = beta[, j], method = "L-BFGS-B",
                lower = lowerage, upper = upperage, control = list(factr = 1)),
                silent = T)
            if (!is(temp, "try-error")) {
                columnoptim[[jj]] = temp
                val[jj] = temp$value
            }
        }
        temp = try(optim(par = n[j], fn = MiAge_fr2, gr = MiAge_grr2,
            b = b, c = c, d = d, betaj = beta[, j], method = "L-BFGS-B",
            lower = lowerage, upper = upperage, control = list(factr = 1)),
            silent = T)
        if (!is(temp, "try-error")) {
            columnoptim[[no.initial.n]] = temp
            val[no.initial.n] = temp$value
        }
        temp = columnoptim[[which(val == min(val, na.rm = T))[1]]]
        if (!is(temp, "try-error")) {
            n[j] = temp$par
        }
        else {
            print(2)
            print(temp)
        }
    }
    return(n)
}

MiAge=mitotic.age(datMeth2[na.omit(match(MiAgedat$MiAge_CpGs$CpGs,rownames(datMeth2))),],MiAgedat$MiAge_parameters[[1]],MiAgedat$MiAge_parameters[[2]],MiAgedat$MiAge_parameters[[3]])
mScore$MiAge=MiAge[as.vector(mScore$SampleID)]

#PedBE The Pediatric-Buccal-Epigenetic (PedBE) clock
PEDBE_CpGs=PEDBE_CpGs[PEDBE_CpGs$ID %in% rownames(datMeth2),]
PedBE=anti.trafo(colSums(PEDBE_CpGs$Coef * datMeth2[as.vector(PEDBE_CpGs$ID),],na.rm=T)-2.10)
mScore$PedBE=PedBE[as.vector(mScore$SampleID)]

#Placental epigenetic clocks
#Control placental clock (CPC)
GACPC_CpGs=GACPC_CpGs[GACPC_CpGs$CpG %in% rownames(datMeth2),]
GACPC=colSums(GACPC_CpGs$coef * datMeth2[as.vector(GACPC_CpGs$CpG),],na.rm = T)+13.06182
mScore$GACPC=GACPC[as.vector(mScore$SampleID)]

#Robust placental clock (RPC)
GARPC_CpGs=GARPC_CpGs[GARPC_CpGs$CpG %in% rownames(datMeth2),]
GARPC=colSums(GARPC_CpGs$coef * datMeth2[as.vector(GARPC_CpGs$CpG),], na.rm = T)+24.99772
mScore$GARPC=GARPC[as.vector(mScore$SampleID)]

#Refined robust placental clock for uncomplicated term pregnancies
GARRPC_CpGs=GARRPC_CpGs[GARRPC_CpGs$CpG %in% rownames(datMeth2),]
GARRPC=colSums(GARRPC_CpGs$coef * datMeth2[as.vector(GARRPC_CpGs$CpG),], na.rm = T)+ 30.74966
mScore$GARRPC=GARRPC[as.vector(mScore$SampleID)]

#Bohlin Gestational age
Bohlin_CpGs=Bohlin_CpGs[Bohlin_CpGs$CpG %in% rownames(datMeth2),]
BohlinGAge=colSums(Bohlin_CpGs$coef * datMeth2[as.vector(Bohlin_CpGs$CpG),],na.rm=T)+ 277.2421
mScore$BohlinGAge=BohlinGAge[as.vector(mScore$SampleID)]

#Knight Gestational age
Knight_CpGs=Knight_CpGs[Knight_CpGs$CpG %in% rownames(datMeth2),]
KnightGAge=colSums(Knight_CpGs$coef * datMeth2[as.vector(Knight_CpGs$CpG),],na.rm=T) + 41.7
mScore$KnightGAge=KnightGAge[as.vector(mScore$SampleID)]

#######################################
#PC clocks
modelcg=names(pcc$pcCpGs)
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=data.frame(cg=names(pcc$pcCpGs),meth_mean=pcc$pcCpGs)
#If user provied reference data include all modelcg, then use UserRef
if(!is.null(UserRef) & ForceUserRef){
    tmp=modelcg[!modelcg %in% as.vector(UserRef$cg)]
    if(length(tmp)==0){refdat=UserRef}else{
            message(paste0("The following model required CpGs are missing in UserRef:"))
            message(paste0(tmp,collaps=" "))
            message("System reference were used for all PC clocks")
    }
}
datMeth2=norm_imputing(datMeth,refdat=refdat,normalize,missedcg)

if(!is.null(datPheno)){
    pcclock=calcPCClocks(methdat=datMeth2,phenodat=datPheno,pcc)
    pcclock=pcclock[,!(names(pcclock) %in% c("Age","Female"))]
    mScore=merge(mScore,pcclock,by="SampleID")
}

#########################################
#cAge and bAge
modelcg=unique(as.vector(cAge_CpG$means$cpg))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
refdat=data.frame(cg=as.vector(cAge_CpG$means$cpg),meth_mean=cAge_CpG$means$mean)
#If user provied reference data include all modelcg, then use UserRef
if(!is.null(UserRef) & ForceUserRef){
    tmp=modelcg[!modelcg %in% as.vector(UserRef$cg)]
    if(length(tmp)==0){refdat=UserRef}else{
            message(paste0("The following model required CpGs are missing in UserRef:"))
            message(paste0(tmp,collaps=" "))
            message("System reference were used for cAge, bAge and bAge_Years")
    }
}
datMeth2=norm_imputing(datMeth,refdat=refdat,normalize,missedcg)

mScore$cAge=calc_cAge(datMeth2,cAge_CpG)[as.vector(mScore$SampleID)]

#######################################
#use episcore reference data for bAge
modelcg=unique(as.vector(bAge_CpG$cpgs$CpG_Site))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
tmp0=episcore_model[!duplicated(episcore_model$CpG_Site),]
refdat=data.frame(cg=tmp0$CpG_Site,meth_mean=tmp0$Mean_Beta_Value)
datMeth2=norm_imputing(datMeth,refdat=refdat,normalize,missedcg)

tmp=calc_bAge(datMeth2,datPheno,mScore,bAge_CpG,GrimAgeComponent)[as.vector(mScore$SampleID),]
mScore$bAge=tmp$bAge
mScore$bAge_Years=tmp$bAge_Years


##########################################
#episcore
modelcg=unique(as.vector(episcore_model$CpG_Site))
missedcg=modelcg[!(modelcg %in% rownames(datMeth))]
tmp0=episcore_model[!duplicated(episcore_model$CpG_Site),]
refdat=data.frame(cg=tmp0$CpG_Site,meth_mean=tmp0$Mean_Beta_Value)
#If user provied reference data include all modelcg, then use UserRef
if(!is.null(UserRef) & ForceUserRef){
    tmp=modelcg[!modelcg %in% as.vector(UserRef$cg)]
    if(length(tmp)==0){refdat=UserRef}else{
            message(paste0("The following model required CpGs are missing in UserRef:"))
            message(paste0(tmp,collaps=" "))
            message("System reference were used for all episcores")
    }
}
datMeth2=norm_imputing(datMeth,refdat=refdat,normalize,missedcg)

out <- data.frame(SampleID=colnames(datMeth2))
for(i in unique(episcore_model$Predictor)){
  tmp_coef = episcore_model[episcore_model$Predictor %in% i, ]
  tmp_coef=tmp_coef[tmp_coef$CpG_Site %in% rownames(datMeth2),]
  if(nrow(tmp_coef) > 1) {
    out[,i]=colSums(tmp_coef$Coefficient * datMeth2[as.vector(tmp_coef$CpG_Site),],na.rm=T)
  } else {
    out[,i] = tmp_coef$Coefficient * datMeth2[as.vector(tmp_coef$CpG_Site),]
  }
}
out$'Epigenetic Age (Zhang)' <- out$'Epigenetic Age (Zhang)' + 65.79295

mScore=merge(mScore,out,by="SampleID")

#####################################
#DNAmFitAge
FitAge=calc_DNAmFitAge(datMeth,datPheno,mScore,GrimAgeComponent,DNAmFitnessModels)
mScore=merge(mScore,FitAge,by="SampleID")

#column order
mScore=mScore[,c("SampleID",as.vector(methscore_dict$predictor))]
###################Age residuals
if(!is.null(datPheno)){
    rownames(datPheno)=datPheno$SampleID
    datPheno=datPheno[as.vector(mScore$SampleID),]
    clockColumns=colnames(mScore)[!(colnames(mScore)=="SampleID")]
    mScore[,"Age"] = datPheno$Age
    for (i in clockColumns){
        mScore[,paste0(i,"Resid")] = resid(lm(mScore[,i] ~ datPheno$Age))
    }
}
return(mScore)
}

#calculation of PC Clocks
#Higgins-Chen et al, PMID 36277076,https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9586209/
calcPCClocks <- function(methdat, phenodat,pcc=NULL){

  if(!("Age" %in% colnames(phenodat))){
    stop("Error: phenodat must have a column named Age")
  }
  if(!("Female" %in% colnames(phenodat))){
    stop("Error: phenodat must have a column named Female")
  }
  if(!("SampleID" %in% colnames(phenodat))){
    stop("Error: phenodat must have a column named SampleID")
  }else if(sum(colnames(methdat) %in% phenodat$SampleID)<ncol(methdat)){
     stop("Error: some samples in methdat do not have phenotype in phenodat")
  }else{
     rownames(phenodat)=phenodat$SampleID
     phenodat=phenodat[colnames(methdat),]
  }

anti.trafo <-function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

  #check number of missing CpGs
  flag=apply(is.na(methdat),1,sum)!=ncol(methdat)
  methdat=methdat[flag,]
#  load(file = paste(path_to_PCClocks_directory,"pcclock_model.RData", sep = ""))
#  if(is.null(pcc)){load(system.file("pcclock_model.RData",package="ENmix"))}
  CpGs=names(pcc$pcCpGs)
  methdat=methdat[rownames(methdat) %in% CpGs,]
  missingCpGs=CpGs[!(CpGs %in% rownames(methdat))]
#  message(paste0("PC clocks: Missing ",length(missingCpGs)," out of ",length(CpGs)," CpGs required for PC clocks calculation"))
  #impute missing probes and values
  tmp=matrix(rep(pcc$pcCpGs[missingCpGs],ncol(methdat)),ncol=ncol(methdat))
  rownames(tmp)=missingCpGs
  methdat=rbind(methdat,tmp)
  meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
  methdat <- apply(methdat,1,meanimpute)
  methdat=methdat[,CpGs]
  phenodat=phenodat[rownames(methdat),]
  #Calculate PC Clocks
  #Initialize a data frame for PC clocks
  DNAmAge <- data.frame(phenodat)
#  message("Calculating PC Clocks now")
  DNAmAge$PCHorvath1 <- as.numeric(anti.trafo(sweep(as.matrix(methdat),2,pcc$PCHorvath1$center) %*% pcc$PCHorvath1$model + pcc$PCHorvath1$intercept))
  DNAmAge$PCHorvath2 <- as.numeric(anti.trafo(sweep(as.matrix(methdat),2,pcc$PCHorvath2$center) %*% pcc$PCHorvath2$model + pcc$PCHorvath2$intercept))
  DNAmAge$PCHannum <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCHannum$center) %*% pcc$PCHannum$model + pcc$PCHannum$intercept)
  DNAmAge$PCPhenoAge <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCPhenoAge$center) %*% pcc$PCPhenoAge$model + pcc$PCPhenoAge$intercept)
  DNAmAge$PCDNAmTL <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCDNAmTL$center) %*% pcc$PCDNAmTL$model + pcc$PCDNAmTL$intercept)
  pp=cbind(Female = DNAmAge$Female,Age = DNAmAge$Age)
  DNAmAge$PCPACKYRS <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCPACKYRS$model + pp %*% pcc$PCPACKYRS$pfa + pcc$PCPACKYRS$intercept)
  DNAmAge$PCADM <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCADM$model + pp %*% pcc$PCADM$pfa + pcc$PCADM$intercept)
  DNAmAge$PCB2M <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCB2M$model + as.matrix(pp[,2]) %*% pcc$PCB2M$pfa + pcc$PCB2M$intercept)
  DNAmAge$PCCystatinC <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCCystatinC$model + pp %*% pcc$PCCystatinC$pfa + pcc$PCCystatinC$intercept)
  DNAmAge$PCGDF15 <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCGDF15$model + as.matrix(pp[,2]) %*% pcc$PCGDF15$pfa + pcc$PCGDF15$intercept)
  DNAmAge$PCLeptin <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCLeptin$model + pp %*% pcc$PCLeptin$pfa + pcc$PCLeptin$intercept)
  DNAmAge$PCPAI1 <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCPAI1$model + as.matrix(pp[,1]) %*% pcc$PCPAI1$pfa + pcc$PCPAI1$intercept)
  DNAmAge$PCTIMP1 <- as.numeric(sweep(as.matrix(methdat),2,pcc$PCGrimAge$center) %*% pcc$PCTIMP1$model + pp %*% pcc$PCTIMP1$pfa + pcc$PCTIMP1$intercept)
  DNAmAge$PCGrimAge <- as.numeric(as.matrix(DNAmAge[,pcc$PCGrimAge$components]) %*% pcc$PCGrimAge$model + pcc$PCGrimAge$intercept)

 return(DNAmAge)
}

#cAge
calc_cAge <- function(data,cAge_CpG){
#load(system.file("cAge.RData",package="ENmix"))
coef_linear=cAge_CpG$coef_linear
coef_log=cAge_CpG$coef_log
intercept=cAge_CpG$intercept
intercept_log=cAge_CpG$intercept_log
means=cAge_CpG$means
rownames(coef_linear)=coef_linear$CpG_Site
rownames(coef_log)=coef_log$CpG_Site

data <- data[rownames(data) %in% rownames(means),]
## impute missing CpGs
if (!nrow(data)==nrow(means)) {
    missing_cpgs <- means[!rownames(means) %in% rownames(data),]
    mat=matrix(rep(missing_cpgs$mean,ncol(data)),ncol=ncol(data))
    rownames(mat)=rownames(missing_cpgs)
    colnames(mat)=colnames(data)
    data <- rbind(data,mat)
}
## impute missing values
na_to_mean <-function(x) {x[is.na(x)] <- mean(x, na.rm=T);x}
data <- t(apply(data,1,na_to_mean))

## Prep for linear predictor
cg1=as.vector(coef_linear$CpG_Site[-grep('_2', coef_linear$CpG_Site)])
cg2=as.vector(coef_linear$CpG_Site[grep('_2', coef_linear$CpG_Site)])
cg2=gsub("_2","",cg2)
tmp <- data[cg2,]^2;rownames(tmp)=paste0(rownames(tmp),"_2")
scores_linear <- rbind(data[cg1,],tmp)

## Prep for log predictor
cg1=as.vector(coef_log$CpG_Site[-grep('_2', coef_log$CpG_Site)])
cg2=as.vector(coef_log$CpG_Site[grep('_2', coef_log$CpG_Site)])
cg2=gsub("_2","",cg2)
tmp <- data[cg2,]^2;rownames(tmp)=paste0(rownames(tmp),"_2")
scores_log <- rbind(data[cg1,],tmp)

## Calculate cAge with linear model
coef_linear <- coef_linear[rownames(scores_linear),]
pred_linear_pp <- colSums(scores_linear * coef_linear$Coefficient) + intercept

## Identify any individuals predicted as under 20s, and re-run with model trained on log(age)
over20s <- names(pred_linear_pp[pred_linear_pp > 20])
pred_linear_pp1 <- pred_linear_pp[over20s]
under20s <- names(pred_linear_pp[pred_linear_pp < 20])

## Now re-run model for those
coef_log <- coef_log[rownames(scores_log),]
pred_log_pp <- colSums(scores_log * coef_log$Coefficient) + intercept_log
pred_log_pp <- exp(pred_log_pp[under20s])

c(pred_log_pp, pred_linear_pp1)
}

#bAge
calc_bAge<-function(data,phenodat,grim,bAge_CpG,GrimAgeComponent){
grim$DNAmGrimAge=grim$PCGrimAge
rownames(grim)=grim$SampleID
rownames(phenodat)=phenodat$SampleID
cpgs=bAge_CpG$cpgs
coefficients=bAge_CpG$coefficients
##standardize CpG by CpG
coef <- data[rownames(data) %in% as.vector(cpgs$CpG_Site),]
ids <- colnames(coef)
coef <- t(apply(coef, 1, scale))
colnames(coef) <- ids

##impute missing CpGs
if (!nrow(coef)==length(unique(cpgs$CpG_Site))) {
    missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)), c("CpG_Site", "Mean_Beta_Value")]
    missing_cpgs=missing_cpgs[!duplicated(as.vector(missing_cpgs$CpG_Site)),]
    mat=matrix(rep(missing_cpgs$Mean_Beta_Value,ncol(coef)),ncol=ncol(coef))
    rownames(mat)=as.vector(missing_cpgs$CpG_Site)
    colnames(mat)=colnames(coef)
    coef = rbind(coef,mat)
}
## impute NA with CpG means
na_to_mean <-function(methyl) {
  methyl[is.na(methyl)] <- mean(methyl, na.rm=T)
  return(methyl)
}
coef <- t(apply(coef,1,function(x) na_to_mean(x)))

#Calculate Episcores
loop <- unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  tmp=coef[as.vector(tmp_coef$CpG_Site),]
  if(nrow(tmp_coef) > 1) {
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp,na.rm=T)
  } else {
    out[colnames(coef),i] = tmp*tmp_coef$Coefficient
  }
}

###### Standadize GrimAge components
samples <- rownames(out)
grim_pred <- grim[samples, c("DNAmGrimAge"), drop = FALSE]
if(!is.null(GrimAgeComponent)){
  grim <- GrimAgeComponent[samples, c("DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS", "DNAmPAI1", "DNAmTIMP1")]
}else{
    grim <- grim[samples, c("PCADM", "PCB2M", "PCCystatinC", "PCGDF15", "PCLeptin", "PCPACKYRS", "PCPAI1", "PCTIMP1")]
}
colnames(grim)=c("DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS", "DNAmPAI1", "DNAmTIMP1")
grim <- scale(grim)
## Calculate bAge
scores <- cbind(Age=phenodat[samples, c("Age")], grim, out)
scores <- scores[, coefficients$Variable]
pred_pp <- colSums(coefficients[,"Coefficient"] * t(scores))

## Scale to same scale as age
scale_pred <- function(x, mean_pred, sd_pred, mean_test, sd_test) {
  scaled <- mean_test + (x - mean_pred)*(sd_test/sd_pred)
  return(scaled)
}
# Scale to same Z scale
scale_Z <- function(x, mean_pred, sd_pred) {
  scaled <- (x - mean_pred)/sd_pred
  return(scaled)
}

mean_pred <- mean(pred_pp)
mean_test <- mean(phenodat$Age) # Mean age in testing data
sd_pred <- sd(pred_pp)
sd_test <- sd(phenodat$Age) # SD age in testing data

bAge <- scale_Z(pred_pp, mean_pred, sd_pred)
bAge_Years <- scale_pred(pred_pp, mean_pred, sd_pred, mean_test, sd_test)

pred=data.frame(SampleID=names(bAge),bAge=bAge,bAge_Years=bAge_Years)
rownames(pred)=pred$SampleID
return(pred)
}


#######DNAmFitAge begin
calc_DNAmFitAge<-function(datMeth,datPheno,mScore,GrimAgeComponent,DNAmFitnessModels){

	if(!is.null(GrimAgeComponent)){grim1=GrimAgeComponent[,c("SampleID","DNAmGrimAge")]}else{
        grim1=mScore[,c("SampleID","PCGrimAge")];names(grim1)[2]="DNAmGrimAge"}
	pheno=merge(datPheno,grim1,by="SampleID")
	pheno=pheno[,c("SampleID","Age","Female","DNAmGrimAge")]
	rownames(pheno)=pheno$SampleID
	pheno=pheno[colnames(datMeth),]

data_prep <- function(dataset,pheno){
  dataset <- dataset[rownames(dataset) %in% DNAmFitnessModels$AllCpGs,]
  if(nrow(dataset) != length(DNAmFitnessModels$AllCpGs)){
      cpgs_toadd <- DNAmFitnessModels$AllCpGs[!DNAmFitnessModels$AllCpGs %in% rownames(dataset)]

      #separate by sex to impute missing CpGs
      output=NULL
      rownames(pheno)=pheno$SampleID
      cid=intersect(colnames(dataset), rownames(pheno))
      dataset=dataset[,cid];pheno=pheno[cid,]
      if(sum(pheno$Female==1)>0){
         dat <- dataset[,pheno$Female == 1]
         dat=rbind(dat,matrix(rep(DNAmFitnessModels$Female_Medians_All[cpgs_toadd],ncol(dat)),ncol=ncol(dat),dimnames=list(cpgs_toadd,colnames(dat))))
         output=cbind(output,dat)
      }
      if(sum(pheno$Female==0)>0){
         dat <- dataset[,pheno$Female == 0]
         dat=rbind(dat,matrix(rep(DNAmFitnessModels$Male_Medians_All[cpgs_toadd],ncol(dat)),ncol=ncol(dat),dimnames=list(cpgs_toadd,colnames(dat))))
         output=cbind(output,dat)
      }

#      print(paste0("Total ", length(cpgs_toadd)," Missing CpGs that are assigned median values from training data"))
      dataset=output
  }
  return(dataset)
}
# Function to provide estimates for any 1 DNAm fitness models
# TidyModel is a specific model in DNAmFitnessModels list
DNAmEstimatorAnyModel <- function(dataset, TidyModel){
  intercept=matrix(rep(1.0, ncol(dataset)),ncol=ncol(dataset),dimnames=list("(Intercept)",colnames(dataset)))
  dataset=rbind(intercept,dataset)
  dataset <- dataset[as.vector(TidyModel$term),]
  dm=dimnames(dataset)
  dataset=matrix(as.numeric(dataset),nrow=nrow(dataset))
  dimnames(dataset)=dm
  estimate <- colSums(TidyModel$estimate  * dataset)
  return(estimate)
}

# Function to calculate all DNAm fitness estimates #
DNAmFitnessEstimators <- function(data, pheno){
  rownames(pheno)=pheno$SampleID
  pheno=pheno[colnames(data),]
  data=rbind(t(pheno),data)

  if(sum(pheno$Female ==1)>0){
  data_fem <- data[,pheno$Female ==1]
  fem_est1 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_noAge_Females) # gait without age
  fem_est2 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_noAge_Females) # grip
  fem_est3 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$VO2maxModel) # vo2max
  fem_est4 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_wAge_Females) # gait w age
  fem_est5 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_wAge_Females) # grip w age
  fem_est6 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$FEV1_wAge_Females) # fev1 w age
  }
  if(sum(pheno$Female ==0)>0){
  data_male <- data[,pheno$Female ==0]
  male_est1 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_noAge_Males) # gait
  male_est2 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_noAge_Males) # grip
  male_est3 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$VO2maxModel) # vo2max
  male_est4 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_wAge_Males) # gait
  male_est5 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_wAge_Males) # grip
  male_est6 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$FEV1_wAge_Males) # fev1
  }
  if(sum(pheno$Female ==1)==0){est1=male_est1;est2=male_est2;est3=male_est3;est4=male_est4;est5=male_est5;est6=male_est6
    }else if(sum(pheno$Female ==0)==0){est1=fem_est1;est2=fem_est2;est3=fem_est3;est4=fem_est4;est5=fem_est5;est6=fem_est6
    }else{
          est1 <- c(fem_est1, male_est1)
          est2 <- c(fem_est2, male_est2)
          est3 <- c(fem_est3, male_est3)
          est4 <- c(fem_est4, male_est4)
          est5 <- c(fem_est5, male_est5)
          est6 <- c(fem_est6, male_est6)
    }
  data_and_est=data.frame(DNAmGait_noAge=est1,DNAmGrip_noAge=est2,DNAmVO2max=est3,DNAmGait_wAge=est4,DNAmGrip_wAge=est5,DNAmFEV1_wAge=est6)
  data_and_est$SampleID=rownames(data_and_est)
  return(data_and_est)
}


# Function to calculate DNAmFitAge
FitAgeEstimator <- function(data){
  # prep dataset for fitage- dataset by sex
  fem <- data[data$Female == 1,]
  male <- data[data$Female == 0,]

  # can only estimate FitAge if all variables are present, remove those without them #
  fem_comcase <- fem[complete.cases(fem), ]
  male_comcase <- male[complete.cases(male), ]

  # Female FitAge
  female_fitest <- 0.1044232 * ((fem_comcase$DNAmVO2max - 46.825091) / (-0.13620215)) +
                0.1742083 * ((fem_comcase$DNAmGrip_noAge - 39.857718) / (-0.22074456)) +
                0.2278776 * ((fem_comcase$DNAmGait_noAge - 2.508547) / (-0.01245682))  +
                0.4934908 * ((fem_comcase$DNAmGrimAge - 7.978487) / (0.80928530))

  # Male FitAge
  male_fitest <- 0.1390346 * ((male_comcase$DNAmVO2max - 49.836389) / (-0.141862925)) +
                0.1787371 * ((male_comcase$DNAmGrip_noAge - 57.514016) / (-0.253179827)) +
                0.1593873 * ((male_comcase$DNAmGait_noAge - 2.349080) / (-0.009380061))  +
                0.5228411 * ((male_comcase$DNAmGrimAge - 9.549733) / (0.835120557))

  fem <- data.frame(fem_comcase, DNAmFitAge = female_fitest)
  male <- data.frame(male_comcase, DNAmFitAge = male_fitest)

  resu <- rbind(fem, male)
  resu=resu[,!names(resu) %in% c("Age","Female","DNAmGrimAge")]
  return(resu)
}

sample_data_prep <- data_prep(dataset = datMeth,pheno=pheno)
sample_data_FitnessEst <- DNAmFitnessEstimators(sample_data_prep, pheno)
sample_data_FitAge_prep <- merge(pheno, sample_data_FitnessEst, by = "SampleID")
FitAge <- FitAgeEstimator(sample_data_FitAge_prep)
return(FitAge)
}
##### end DNAmFitAge

