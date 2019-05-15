#get probe info from rgDataSet
getCGinfo <-function(rgSet,type="IandII"){
    if(!is(rgSet, "rgDataSet")){stop("[getCGinfo] The input object must be a rgDataSet")}
    if(!(type %in% c("I","II","IandII","ctrl"))){stop("[getRawMeth] Type must be I,II,IandII,or ctrl")}
    probeinfo=rowData(rgSet)
probeIA=probeinfo[probeinfo$Infinium_Design_Type %in% c("IA"),]
probeIB=probeinfo[probeinfo$Infinium_Design_Type %in% c("IB"),]
names(probeIA)[which(names(probeIA)=="Address")]="AddressA"
probeIA$Infinium_Design_Type="I"
probeIB=probeIB[,c("Name","Address")]
names(probeIB)=c("Name","AddressB")
probeI=merge(probeIA,probeIB,by="Name")
var1=c("Name","AddressA","AddressB","Infinium_Design_Type","Color_Channel","chr","pos","Relation_to_Island")
var2=names(probeinfo)[!(names(probeinfo) %in% c(var1,"Address"))]
var=c(var1,var2)
probeI=probeI[,var]

probesnpIA=probeinfo[probeinfo$Infinium_Design_Type %in% c("snpIA"),]
probesnpIB=probeinfo[probeinfo$Infinium_Design_Type %in% c("snpIB"),]
names(probesnpIA)[which(names(probesnpIA)=="Address")]="AddressA"
probesnpIB=probesnpIB[,c("Name","Address")]
names(probesnpIB)=c("Name","AddressB")
probesnpI=merge(probesnpIA,probesnpIB,by="Name")
probesnpI=probesnpI[,var]
probesnpI$Infinium_Design_Type="snpI"

probeI=rbind(probeI,probesnpI)

probeII=probeinfo[probeinfo$Infinium_Design_Type %in% c("II","snpII"),]
names(probeII)[which(names(probeII)=="Address")]="AddressA"
probeII$AddressB=probeII$AddressA
probeII=probeII[,var]

probectrl=probeinfo[probeinfo$Infinium_Design_Type %in% c("ctrl"),]
probectrl=merge(probectrl,metadata(rgSet)$ictrl,by="Address")

if(type=="I"){
return(probeI)
}else if(type=="II"){
return(probeII)
}else if (type=="ctrl"){
return(probectrl)
}else if (type=="IandII"){
return(rbind(probeI,probeII))
}
}
