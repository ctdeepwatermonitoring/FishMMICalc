setwd("")

library(dplyr)

siteinfo<-read.csv("SiteInfo.csv",header=TRUE)
siteinfo$km2<-siteinfo$CUMDRAINAGE_SQMI*2.58999
siteinfo<-siteinfo[,c(1,8)]
names(siteinfo)<-c("SID","km2")

ecochar<-read.csv("Fish_EcoChar.csv",header=TRUE)
mastertaxa<- read.csv("fishmastertaxa_echochar.csv",header=TRUE)

fish<-read.csv("WPLR_FishData_2014.csv",header=TRUE)##CSV output from AWQ


#Fields needed in Fish data to run calcs
#"STA_SEQ","Station_Name","ï..Taxon_Name","Result_Value","Result_UOM","Pass_Count",
#"IsStocked","Reach_Length","Reach_Length_UOM","Width_Length","Width_Length_UOM","Lab_Accession"
#Fields renamed below for consistency amoung datasets
fish<-fish[,c(24,25,1:3,7,13:18)]
names(fish)<-c("SID","name","taxa_name","result","resultUOM","pass_cnt",
               "isstocked","reach_len","reachUOM","width_len","widthUOM","lab_accession")


#Get unique lab accesssion in and merge with needed site information to obtain drainage
site<- unique(fish[c("lab_accession", "SID","reach_len","width_len")])
site<- merge(site,siteinfo,by="SID")

##Make sure that pass count is 1, length measures are in meters and remove stocked fish
fish<-fish[fish$pass_cnt==1&fish$isstocked=="False"&fish$reachUOM=="m"&fish$widthUOM=="m",]
fish<- fish %>%
        group_by(lab_accession,taxa_name) %>%
        summarise(TotalTaxa=sum(result))

##Merge total ind fish datawith needed info for calcs:  US Watershed, length, width, autecology
fish<-merge(fish,site,by="lab_accession")
fish<-merge(fish,ecochar,by.x="taxa_name",all.x=TRUE)
fish$LW<-(fish$reach_len*fish$width_len)/100
fish<-fish[fish$TotalTaxa>0,]
site<-unique(fish[c(2,4:7)])

#Total Fish Calc used in percentage calcs
TotFish<- fish %>%
  group_by(lab_accession) %>%
  summarise(TotalFish=sum(TotalTaxa))


############################################################################
#############COLD WATER MMI CALCULATIONS####################################
############################################################################

##CW1 N of WBK Per 100 m2
WBK<-fish[fish$taxa_name == "Salvelinus fontinalis",]
WBK<-merge(site[1],WBK,by="lab_accession",all.x=TRUE)
WBK<-WBK[,c(1:3,14)]
WBK[is.na(WBK)] <- 0
WBK$WBK<-WBK$TotalTaxa/WBK$LW
CW1<-WBK[,c(1,5)]
CW1[is.na(CW1)]<-0
CW1<-CW1[order(CW1$lab_accession),]
CW1$CW1<-(CW1$WBK/60.6)*100

##CW2  Percent fluvial dependant individuals
FD<- fish[fish$Stream.flow=="FD",] %>%
      group_by(lab_accession) %>%
      summarise (FDTot=sum(TotalTaxa))
FD<-merge(TotFish,FD,by="lab_accession",all.x=TRUE)
FD[is.na(FD)]<-0
FD$FD<-(FD$FDTot/FD$TotalFish)
CW2<-FD[,c(1,4)]
CW2[is.na(CW2)]<-0
CW2<-CW2[order(CW2$lab_accession),]
CW2$CW2<-((0.717-CW2$FD)/0.717)*100

##CW3  N warm water spp corrected for area
WW<-fish[which(fish$Temp =="W"),]
WW<- WW %>%
      group_by(lab_accession) %>%
      summarize(n = n())
WW<- merge(site[,c(1,5)],WW,by="lab_accession",all.x=TRUE)
WW[is.na(WW)] <- 0
WW$E <- 0.0229+(0.2148*WW$km2)
WW$R <- WW$n - WW$E
CW3<- WW[,c(1,5)]
CW3[is.na(CW3)]<-0
CW3<-CW3[order(CW3$lab_accession),]
CW3$CW3<-((3.0612-CW3$R)/5.453)*100

##CW4  Percent of warm water individuals
PctWW<- fish[which(fish$Temp =="W"),] %>%
          group_by(lab_accession) %>%
          summarise(WWTot = sum(TotalTaxa))
PctWW<-merge(TotFish,PctWW,by="lab_accession",all.x=TRUE)
PctWW[is.na(PctWW)]<-0
PctWW$PctWW<-(PctWW$WWTot/PctWW$TotalFish)
CW4<-PctWW[,c(1,4)]
CW4[is.na(CW4)]<-0
CW4<-CW4[order(CW4$lab_accession),]
CW4$CW4<-((0.875-CW4$PctWW)/0.875)*100

##CW5  Percent of wild brook trout individuals
PctWBK<- fish[fish$taxa_name == "Salvelinus fontinalis",] %>%
  group_by(lab_accession) %>%
  summarise(WBKTot = sum(TotalTaxa))
PctWBK<-merge(TotFish,PctWBK,by="lab_accession",all.x=TRUE)
PctWBK[is.na(PctWBK)]<-0
PctWBK$PctWBK<-(PctWBK$WBKTot/PctWBK$TotalFish)
CW5<-PctWBK[,c(1,4)]
CW5[is.na(CW5)]<-0
CW5<-CW5[order(CW5$lab_accession),]
CW5$CW5<-(CW5$PctWBK/0.863)*100

##Join Metric Calcs together by lab accession
CW_MMI<- Reduce(function(DF1,DF2)merge(DF1,DF2,by="lab_accession",all.x=TRUE),list(CW1,CW2,CW3,CW4,CW5,TotFish))
CW_MMI$CW1<-ifelse(CW_MMI$CW1<0,0,ifelse(CW_MMI$CW1>100,100,CW_MMI$CW1))
CW_MMI$CW2<-ifelse(CW_MMI$CW2<0,0,ifelse(CW_MMI$CW2>100,100,CW_MMI$CW2))
CW_MMI$CW3<-ifelse(CW_MMI$CW3<0,0,ifelse(CW_MMI$CW3>100,100,CW_MMI$CW3))
CW_MMI$CW4<-ifelse(CW_MMI$CW4<0,0,ifelse(CW_MMI$CW4>100,100,CW_MMI$CW4))
CW_MMI$CW5<-ifelse(CW_MMI$CW5<0,0,ifelse(CW_MMI$CW5>100,100,CW_MMI$CW5))

CW_MMI$CW_MMI<-(CW_MMI$CW1+CW_MMI$CW2+CW_MMI$CW3+CW_MMI$CW4+CW_MMI$CW5)/5
CW_MMI<-merge(CW_MMI,site[,c(1:2,5)],by="lab_accession",all.x=TRUE)
CW_MMI$Flag<- ifelse(CW_MMI$TotalFish<20,1,0)

############################################################################
#############MIXED WATER MMI CALCULATIONS###################################
############################################################################

##MW1  Percent white sucker individuals
PctWS<- fish[which(fish$taxa_name=="Catostomus commersoni"|fish$taxa_name=="Catostomus commersonii"),] %>%
          group_by(lab_accession)%>%
          summarise(WSTot = sum(TotalTaxa))
PctWS<-merge(TotFish,PctWS,by="lab_accession",all.x=TRUE)
PctWS[is.na(PctWS)]<-0
PctWS$PctWS<-(PctWS$WSTot/PctWS$TotalFish)
MW1<-PctWS[,c(1,4)]
MW1[is.na(MW1)]<-0
MW1<-MW1[order(MW1$lab_accession),]
MW1$MW1<-((0.439-MW1$PctWS)/0.439)*100


##MW2  Percent Cyprinidae individuals
CYP<- mastertaxa[which(mastertaxa$FAMILY=="CYPRINIDAE"),]
CYP<-CYP[,1:2]
PctCYP<-  merge(fish,CYP,by="taxa_name") %>%
          group_by(lab_accession)%>%
          summarise(CYPTot = sum(TotalTaxa))
PctCYP<-merge(TotFish,PctCYP,by="lab_accession",all.x=TRUE)
PctCYP[is.na(PctCYP)]<-0
PctCYP$PctCYP<-(PctCYP$CYPTot/PctCYP$TotalFish)
MW2<-PctCYP[,c(1,4)]
MW2[is.na(MW2)]<-0
MW2<-MW2[order(MW2$lab_accession),]
MW2$MW2<-((MW2$PctCYP-0.002)/0.935)*100

##MW3  Percent Fluvial specialist individuals without blacknose dace
PctFS<- fish[which(fish$Stream.flow =="FS"&fish$taxa_name!="Rhinichthys atratulus"),] %>%
  group_by(lab_accession) %>%
  summarise(FSTot = sum(TotalTaxa))
PctFS<-merge(TotFish,PctFS,by="lab_accession",all.x=TRUE)
PctFS[is.na(PctFS)]<-0
PctFS$PctFS<-(PctFS$FSTot/PctFS$TotalFish)
MW3<-PctFS[,c(1,4)]
MW3[is.na(MW3)]<-0
MW3<-MW3[order(MW3$lab_accession),]
MW3$MW3<-(MW3$PctFS/0.647)*100

##MW4  Percent Nontolerant generalist feeders individuals
PctIGF<- fish[which(fish$Tolerance =="I"|fish$Tolerance=="M"& fish$Trophic.class=="GF"),]%>%
        group_by(lab_accession)%>%
        summarise(IGFTot = sum(TotalTaxa))
PctIGF<-merge(TotFish,PctIGF,by="lab_accession",all.x=TRUE)
PctIGF[is.na(PctIGF)]<-0
PctIGF$PctIGF<-(PctIGF$IGFTot/PctIGF$TotalFish)
MW4<-PctIGF[,c(1,4)]
MW4[is.na(MW4)]<-0
MW4<-MW4[order(MW4$lab_accession),]
MW4$MW4<-(MW4$PctIGF/0.516)*100

##MW5  Percent Native Warm Water Individuals
PctNWW<- fish[which(fish$Origin =="N"& fish$Temp=="W"),]%>%
  group_by(lab_accession)%>%
  summarise(NWWTot = sum(TotalTaxa))
PctNWW<-merge(TotFish,PctNWW,by="lab_accession",all.x=TRUE)
PctNWW[is.na(PctNWW)]<-0
PctNWW$PctNWW<-(PctNWW$NWWTot/PctNWW$TotalFish)
MW5<-PctNWW[,c(1,4)]
MW5[is.na(MW5)]<-0
MW5<-MW5[order(MW5$lab_accession),]
MW5$MW5<-((0.679-MW5$PctNWW)/0.679)*100

##MW6  Percent Intolerant Individuals
PctITL<- fish[which(fish$Tolerance =="I"),]%>%
  group_by(lab_accession)%>%
  summarise(ITLTot = sum(TotalTaxa))
PctITL<-merge(TotFish,PctITL,by="lab_accession",all.x=TRUE)
PctITL[is.na(PctITL)]<-0
PctITL$PctITL<-(PctITL$ITLTot/PctITL$TotalFish)
MW6<-PctITL[,c(1,4)]
MW6[is.na(MW6)]<-0
MW6<-MW6[order(MW6$lab_accession),]
MW6$MW6<-(MW6$PctITL/0.381)*100

##  N Fluvial specialist species
FS<-fish[which(fish$Stream.flow =="FS"),]
FS<- FS %>%
  group_by(lab_accession) %>%
  summarize(FSn = n())
FS<-merge(TotFish,FS,by="lab_accession",all.x=TRUE)
FS[is.na(FS)]<-0
MW7<-FS[,c(1,3)]
MW7<-MW7[order(MW7$lab_accession),]
MW7$MW7<-((MW7$FSn-1)/4)*100

##Join Metric Calcs together by lab accession
MW_MMI<- Reduce(function(DF1,DF2)merge(DF1,DF2,by="lab_accession",all.x=TRUE),
                list(MW1,MW2,MW3,MW4,MW5,MW6,MW7,TotFish))
MW_MMI$MW1<-ifelse(MW_MMI$MW1<0,0,ifelse(MW_MMI$MW1>100,100,MW_MMI$MW1))
MW_MMI$MW2<-ifelse(MW_MMI$MW2<0,0,ifelse(MW_MMI$MW2>100,100,MW_MMI$MW2))
MW_MMI$MW3<-ifelse(MW_MMI$MW3<0,0,ifelse(MW_MMI$MW3>100,100,MW_MMI$MW3))
MW_MMI$MW4<-ifelse(MW_MMI$MW4<0,0,ifelse(MW_MMI$MW4>100,100,MW_MMI$MW4))
MW_MMI$MW5<-ifelse(MW_MMI$MW5<0,0,ifelse(MW_MMI$MW5>100,100,MW_MMI$MW5))
MW_MMI$MW6<-ifelse(MW_MMI$MW6<0,0,ifelse(MW_MMI$MW6>100,100,MW_MMI$MW6))
MW_MMI$MW7<-ifelse(MW_MMI$MW7<0,0,ifelse(MW_MMI$MW7>100,100,MW_MMI$MW7))

MW_MMI$MW_MMI<-(MW_MMI$MW1+MW_MMI$MW2+MW_MMI$MW3+MW_MMI$MW4+MW_MMI$MW5+MW_MMI$MW6+MW_MMI$MW7)/7
MW_MMI<-merge(MW_MMI,site[,c(1:2,5)],by="lab_accession",all.x=TRUE)
MW_MMI$Flag<- ifelse(MW_MMI$TotalFish<20,1,0)

##Merge Cold Water and Mixed Water for Comparison##
Fish_MMI<- merge(CW_MMI,MW_MMI,by=c("lab_accession","TotalFish","SID","km2","Flag"))
write.csv(Fish_MMI,"Fish_MMItest.csv")
