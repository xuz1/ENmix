
ctrlsva<-function(rgSet,percvar=0.9,npc=1,flag=1)
{
    if(!is(rgSet, "RGChannelSet")){stop("Input needs to an RGChannelSet\n")}
    if(percvar<0 | percvar>1){stop("Percentage of variation threshold needs 
    to be between 0 and 1\n")}

    ctrls<-getProbeInfo(rgSet,type="Control")
    ctrls <- ctrls[ctrls$Address %in% featureNames(rgSet),]
    ctrl_r <- getRed(rgSet)[ctrls$Address[!(ctrls$Type %in% "NEGATIVE")],]
    ctrl_g <- getGreen(rgSet)[ctrls$Address[!(ctrls$Type %in% "NEGATIVE")],]
    ctrl_nneg=rbind(ctrl_r,ctrl_g)
#PCA of control-probe intensities
    pca <- prcomp(t(ctrl_nneg))
    eigenvalue=pca$sdev^2
    perc=eigenvalue/sum(eigenvalue)
    if(flag==1){
    npc=1
    while(sum(perc[1:npc])<percvar){npc <- npc+1}
    npc
    ctrlsva=pca$x[,1:npc]
    }else{ctrlsva=pca$x[,1:npc]}
    cat(npc," surrogate variables explain ",sum(perc[1:npc])*100, "% of 
    data variation\n")
    ctrlsva
}


