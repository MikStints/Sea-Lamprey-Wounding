
#Lamprey data request 
#February 2025
#Mikayla Stinson, original codes from McKayla Jarvie

# Install the development version from GitHub
install.packages("devtools")
devtools::install_github("HoldenJe/gfsR")

library(tidyverse)
library(gfsR)
library(dplyr)
library(lubridate)
library(glfishr)
library(glisDbTools)

################### add code before completing 2025 export to check for duplicates when several tags are recorded on a single fish?
#occurs because of the way fn125 table is merged with fn125_tags table
#or just check manually after csv is saved?

#also add/update code to remove duplicates in the individual OR final gear table(s) - or do manually after csv is saved?

#adjust code to remove projects not completed in 2025/add new projects


############################# FLEN to TLEN conversions ################################################################################################
#see if lko coefficients are different from whole lake attributes

#whole lake
species_attrs <- get_species(list(detail = TRUE))

#lko only
fishtest = get_FN125(list(lake="ON",flen__not_null =TRUE,tlen__not_null=TRUE))
unique(fishtest$SPC)
unique(fishtest$PRJ_CD)
library(dplyr)
library(purrr)
library(tidyr)

get_coefficients <- function(data) {
  model <- lm(TLEN ~ FLEN, data = data)
  coef <- coef(model)
  return(data.frame(Intercept = round(coef[1],5), Slope = round(coef[2],5)))
}

results <- fishtest %>%
  group_by(SPC) %>%
  nest() %>%
  mutate(Coefficients = map(data, get_coefficients)) %>%
  unnest(Coefficients) %>%
  dplyr::select(SPC, Intercept, Slope)

results = results%>%
  rename(FLEN2TLEN_ALPHA_LKO=Intercept, FLEN2TLEN_BETA_LKO=Slope)

#combine to see diffs
compare.attr = left_join(results,species_attrs,by="SPC")
compare.attr = compare.attr %>%
  dplyr::select(SPC,FLEN2TLEN_ALPHA_LKO,FLEN2TLEN_ALPHA,FLEN2TLEN_BETA_LKO,FLEN2TLEN_BETA)
#note that some species have slope=NA because LKO only has one record where both TLEN and FLEN were recorded... (i.e., 271,380,106,173)
#write.csv(compare.attr, "./data/exports/FLEN_TLEN_coefficients.csv")

#FLEN to TLEN conversion function from Adam Cottrill if needed for certain programs - using LKO coefficients only
#had to edit the original code as it didn't return anything.

estimate_tlen2 <- function(GL1na_rows, results) {
  # flag which rows are being estimated (based on original TLEN)
  est_flag <- is.na(GL1na_rows$TLEN)
  alpha_beta <- results[, c("SPC", "FLEN2TLEN_ALPHA_LKO", "FLEN2TLEN_BETA_LKO")]
  alpha_beta <- unique(alpha_beta)
  GL1na_rows <- merge(GL1na_rows, alpha_beta, by = "SPC", all.x = TRUE)
  GL1na_rows$TLEN <- ifelse(
    is.na(GL1na_rows$TLEN),
    GL1na_rows$FLEN2TLEN_ALPHA_LKO + GL1na_rows$FLEN2TLEN_BETA_LKO * GL1na_rows$FLEN,
    GL1na_rows$TLEN
  )
  GL1na_rows$EstTL <- est_flag
  GL1na_rows <- GL1na_rows[, !names(GL1na_rows) %in% c("FLEN2TLEN_ALPHA_LKO", "FLEN2TLEN_BETA_LKO")]
  return(GL1na_rows)
}
##################### GL1 ############################ #check if any TLENs are missing -  yes there are lengths missing for 2025########################################
#pull in from glis
Lamprey <- get_FN125_Lamprey(list(prj_cd=c("LOA_IA25_GL1")))
fn121 <- get_FN121(list(prj_cd=c("LOA_IA25_GL1")))
fn122 <- get_FN122(list(prj_cd=c("LOA_IA25_GL1")))
fn125 <- get_FN125(list(prj_cd=c("LOA_IA25_GL1")))
fn125tags <- get_FN125_Tags(list(prj_cd=c("LOA_IA25_GL1"))) #empty
##no tags recorded for 2025


#this pulls out the rows with NA values in TLEN
GL1na_rows = fn125[is.na(fn125$TLEN),]

#this completely removes the NA values in TLEN 
fn125_noNA = fn125 %>%
  filter(!is.na(TLEN))

#If there are NA values then convert FLEN to TLEN
fn125.convert3 = estimate_tlen2(GL1na_rows, results)

#rbind the converted table with the larger data set
updated_fn125 = bind_rows(fn125_noNA, fn125.convert3)
##there are two rows that still have NA because there was no FLEN to convert (2025)


#biodataexport
fn125.convert<-append.spc.names(updated_fn125) #if no conversion needed, this would just be fn125

#if there are tags then run the next line of code
#fn125edit<-left_join(fn125.convert,fn125tags) #no tags present for 2025

fn125edit<-fn125.convert%>%dplyr::select(PRJ_CD,SAM,EFF,GRP,SPC,SPC_NM,FISH,TLEN,RWT,AGEST,SEX,MAT,CLIPC) 
#if there are tags for that year, replace fn125.convert with fn125edit. (line above)


#adding the lampray data
lampreyspread<-Lamprey%>%
  group_by(PRJ_CD,SAM,SPC,EFF,GRP,FISH,LAMIJC_TYPE,COMMENT_LAM)%>%
  dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),Lake="Ontario",Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),"R/D"="D","A1-A3" = A1+A2+A3) 

biodataexport<-left_join(lampreyspread,fn125edit)%>%
  dplyr::mutate(EFF=as.numeric(EFF),
                SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),
                Age="",SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                Comments=COMMENT_LAM) 



biodataexport$B1 = 0 #add these strings if there are no values for the given year
biodataexport$B4 = 0 #we add these to create columns that aren't pulled from the first lampray table

biodataexport1 = biodataexport%>%
  dplyr::select(FishID,
                LiftID,
                Lake,Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                "A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)

write.csv(biodataexport1,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/GL1/LOA_IA25_GL1_bio.csv", 
          row.names = FALSE)


#now do "gear" table export
tem<-fn122%>%
  dplyr::select(PRJ_CD,EFF,SAM,GRTEM0)
fn121<-left_join(fn121,tem)%>%
  dplyr::rename(BottomTempC=GRTEM0)
gearexport1<-fn121%>%
  dplyr::filter(SAM%in%lampreyspread$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
         SurveyType="Fish Community Assessment",
         SurveyDescription="Fish Community Assessment",
         Gear="Gillnet",
         Nights=1,
         NetMaterial="M",
         NetLengthkm=ifelse(MODE=="01",0.142,0.152),
         MinMeshmm=38.1,
         MaxMeshmm=152.4,
         MU=NA,
         Grid=NA)%>%
  dplyr::select(c(LiftID,Lake,Agency,SUBSPACE,DD_LAT0,DD_LON0,MU,Grid,Year,Month,Day,SurveyType,SurveyDescription,Gear,Nights,NetMaterial,MinMeshmm,MaxMeshmm,GRDEPMIN,GRDEPMAX,GRDEPMID,BottomTempC,SITEM0,COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,Latitude=DD_LAT0,Longitude=DD_LON0,Depth1m=GRDEPMIN,Depth2m=GRDEPMAX,AvgDepthm=GRDEPMID,SurfaceTempC=SITEM0,Comments=COMMENT1)

#write.csv

write.csv(gearexport1,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/GL1/LOA_IA25_GL1_gear.csv", 
          row.names = FALSE)

##################### TW1 - 2025 ############################ #check if any TLENs are missing -  #############################################
#pull in from glis
Lamprey_trawl <- get_FN125_Lamprey(list(prj_cd=c("LOA_IA25_TW1")))
fn121_trawl <- get_FN121(list(prj_cd=c("LOA_IA25_TW1")))
fn122_trawl <- get_FN122(list(prj_cd=c("LOA_IA25_TW1")))
fn125_trawl <- get_FN125(list(prj_cd=c("LOA_IA25_TW1")))
fn125tags_trawl <- get_FN125_Tags(list(prj_cd=c("LOA_IA25_TW1"))) #empty

#check if NA values in TLEN
TW1na_rows = fn125_trawl[is.na(fn125_trawl$TLEN),]
#one fish with NA but doesn't have a FLEN so unable to convert

#biodataexport
fn125_trawl <-append.spc.names(fn125_trawl)

#run next code if there's any tags
#fn125edit<-left_join(fn125,fn125tags)

fn125edit<-fn125_trawl%>%dplyr::select(PRJ_CD,SAM,EFF,GRP,SPC,SPC_NM,FISH,TLEN,RWT,AGEST,SEX,MAT,CLIPC)#,TAGID)

lampreyspread<-Lamprey_trawl%>%
  group_by(PRJ_CD,SAM,SPC,EFF,GRP,FISH,LAMIJC_TYPE,COMMENT_LAM)%>%
  dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
         Lake="Ontario",Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),"R/D"="D","A1-A3" = 0) #no A1/A2/A3 present this year, otherwise A1+A2+A3

biodataexport_trawl<-left_join(lampreyspread,fn125edit)%>%
  dplyr::mutate(EFF=as.numeric(EFF),
                SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),Age="",
                SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                Comments=COMMENT_LAM)

biodataexport_trawl$CWTAgency = ""
biodataexport_trawl$A1 = 0 #add since there weren't occurrences this year
biodataexport_trawl$A2 = 0 #add since there weren't occurrences this year
biodataexport_trawl$A3 = 0 #add since there weren't occurrences this year
biodataexport_trawl$A4 = 0 #add since there weren't occurrences this year
biodataexport_trawl$B1 = 0 #add since there weren't occurrences this year
biodataexport_trawl$B2 = 0 #add since there weren't occurrences this year
biodataexport_trawl$B3 = 0 #add since there weren't occurrences this year
biodataexport_trawl$B4 = 0 #add since there weren't occurrences this year

biodataexport_trawl2 = biodataexport_trawl%>%
  dplyr::select(FishID,
                LiftID,
                Lake,
                Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                CWTAgency,
                "A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)


#write.csv(biodataexport2,"./data/exports/LOA_IA25_TW1_bio.csv")
write.csv(biodataexport_trawl2,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/TW1/LOA_IA25_Tw1_bio.csv", 
          row.names = FALSE)


#now do "gear" table export
tem<-fn122_trawl%>%dplyr::select(PRJ_CD,EFF,SAM,GRTEM0)

fn121_trawl<-left_join(fn121_trawl,tem)%>%dplyr::rename(BottomTempC=GRTEM0)

gearexport_trawl<-fn121_trawl%>%
  dplyr::filter(SAM%in%lampreyspread$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                SurveyType="Fish Community Assessment",
                SurveyDescription="Fish Community Assessment",
                Gear="Trawl",
                Nights=NA,
                NetMaterial="M",
                NetLengthkm=NA,
                MinMeshmm=88.9,
                MaxMeshmm=101.6,
                MU=NA,Grid=NA)%>%
  dplyr::select(c(LiftID,
                  Lake,Agency,
                  SUBSPACE,
                  DD_LAT0,
                  DD_LON0,
                  MU,
                  Grid,
                  Year,
                  Month,
                  Day,
                  SurveyType,
                  SurveyDescription,
                  Gear,
                  Nights,
                  NetMaterial,
                  MinMeshmm,
                  MaxMeshmm,
                  GRDEPMIN,
                  GRDEPMAX,
                  GRDEPMID,
                  BottomTempC,
                  SITEM0,
                  COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,
                Latitude=DD_LAT0,
                Longitude=DD_LON0,
                Depth1m=GRDEPMIN,
                Depth2m=GRDEPMAX,
                AvgDepthm=GRDEPMID,
                SurfaceTempC=SITEM0,
                Comments=COMMENT1)
#write.csv(gearexport_trawl,"./data/exports/LOA_IA25_TW1_gear.csv")
write.csv(gearexport_trawl,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/TW1/LOA_IA25_Tw1_gear.csv", 
          row.names = FALSE)

##################### NSH/NSW/NSE/NSA - 2025 ############################ #check if any TLENs are missing - yes this year ###########################################
#pull in from glis
Lamprey_NSCIN <- get_FN125_Lamprey(list(prj_cd=c("LOA_IA25_NSH","LOA_IA25_NSW","LOA_IA25_NSE","LOA_IA25_NSA")))
fn121_NSCIN <- get_FN121(list(prj_cd=c("LOA_IA25_NSH","LOA_IA25_NSW","LOA_IA25_NSE","LOA_IA25_NSA")))
fn122_NSCIN <- get_FN122(list(prj_cd=c("LOA_IA25_NSH","LOA_IA25_NSW","LOA_IA25_NSE","LOA_IA25_NSA")))
fn125_NSCIN <- get_FN125(list(prj_cd=c("LOA_IA25_NSH","LOA_IA25_NSW","LOA_IA25_NSE","LOA_IA25_NSA")))
fn125tags_NSCIN <- get_FN125_Tags(list(prj_cd=c("LOA_IA25_NSH","LOA_IA25_NSW","LOA_IA25_NSE","LOA_IA25_NSA")))


#check if NA values in TLEN
NSCINna_rows = fn125_NSCIN[is.na(fn125_NSCIN$TLEN),]
#multiple missing TLEN

#this completely removes the NA values in TLEN 
NSCIN_fn125_noNA = fn125_NSCIN %>%
  filter(!is.na(TLEN))

#If there are NA values then convert FLEN to TLEN
NSCIN_fn125.convert = estimate_tlen2(NSCINna_rows, results)
#8 records have no FLEN or TLEN

#rbind the converted table with the larger data set
NSCIN_updated_fn125 = bind_rows(NSCIN_fn125_noNA, NSCIN_fn125.convert)




#biodataexport
NSCIN_fn125.convert<-append.spc.names(NSCIN_updated_fn125) 

NSCIN_fn125edit<-left_join(NSCIN_fn125.convert,fn125tags_NSCIN)

NSCIN_fn125edit2<-NSCIN_fn125edit%>%
  dplyr::select(PRJ_CD,SAM,EFF,GRP,SPC,SPC_NM,FISH,TLEN,RWT,AGEST,SEX,MAT,CLIPC,TAGDOC) 

NSCIN_lampreyspread<-Lamprey_NSCIN%>%
  group_by(PRJ_CD,SAM,SPC,EFF,GRP,FISH,LAMIJC_TYPE,COMMENT_LAM)%>%
  dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
         Lake="Ontario",
         Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),
         "R/D"="D","A1-A3" = A1) #no A2/A3 this year, otherwise A1+A2+A3

biodataexport_NSCIN<-left_join(NSCIN_lampreyspread,NSCIN_fn125edit2)%>%
  dplyr::mutate(EFF=as.numeric(EFF),
                SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),
                Age="",SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                CWTAgency=TAGDOC,
                Comments=COMMENT_LAM) #if no conversion, Length = TLEN

biodataexport_NSCIN$A2 = 0 #add since there weren't occurrences this year
biodataexport_NSCIN$A3 = 0 #add since there weren't occurrences this year
biodataexport_NSCIN$A4 = 0 #add since there weren't occurrences this year
biodataexport_NSCIN$B1 = 0 #add since there weren't occurrences this year
biodataexport_NSCIN$B2 = 0 #add since there weren't occurrences this year
biodataexport_NSCIN$B4 = 0 #add since there weren't occurrences this year

biodataexport_NSCIN_final = biodataexport_NSCIN%>%
  dplyr::select(FishID,
                LiftID,
                Lake,
                Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                CWTAgency,"A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)

#write.csv
write.csv(biodataexport_NSCIN_final,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/NSCIN/LOA_IA25_NSH_NSW_NSE_NSA_bio.csv", 
          row.names = FALSE)

#now do "gear" table export
tem2<-fn122_NSCIN%>%dplyr::select(PRJ_CD,EFF,SAM,GRTEM0)
fn121_NSCIN<-left_join(fn121_NSCIN,tem2)%>%dplyr::rename(BottomTempC=GRTEM0)
gearexport_NSCIN<-fn121_NSCIN%>%dplyr::filter(SAM%in%NSCIN_lampreyspread$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                SurveyType="Fish Community Assessment",
                SurveyDescription="Nearshore Community Index Netting Program",
                Gear="Trapnet",
                Nights=NA,
                NetMaterial=NA,
                NetLengthkm=NA,
                MinMeshmm=NA,
                MaxMeshmm=NA,
                MU=NA,
                Grid=NA)%>%
  dplyr::select(c(LiftID,
                  Lake,
                  Agency,
                  SUBSPACE,
                  DD_LAT0,
                  DD_LON0,
                  MU,
                  Grid,
                  Year,
                  Month,
                  Day,
                  SurveyType,
                  SurveyDescription,
                  Gear,
                  Nights,
                  NetMaterial,
                  MinMeshmm,
                  MaxMeshmm,
                  GRDEPMIN,
                  GRDEPMAX,
                  GRDEPMID,
                  BottomTempC,
                  SITEM1,
                  COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,
                Latitude=DD_LAT0,
                Longitude=DD_LON0,
                Depth1m=GRDEPMIN,
                Depth2m=GRDEPMAX,
                AvgDepthm=GRDEPMID,
                SurfaceTempC=SITEM1,
                Comments=COMMENT1)

#write.csv-gearexport_NSCIN
write.csv(gearexport_NSCIN,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/NSCIN/LOA_IA25_NSH_NSW_NSE_NSA_gear.csv", 
          row.names = FALSE)


##################### TW2/TW4 - 2025 ############################ #check if any TLENs are missing, no this year #############################################
#pull in from glis
Lamprey_TW <- get_FN125_Lamprey(list(prj_cd=c("LOA_IA25_TW2","LOA_IA25_TW4")))
fn121_TW <- get_FN121(list(prj_cd=c("LOA_IA25_TW2","LOA_IA25_TW4")))
fn122_TW <- get_FN122(list(prj_cd=c("LOA_IA25_TW2","LOA_IA25_TW4")))
fn125_TW <- get_FN125(list(prj_cd=c("LOA_IA25_TW2","LOA_IA25_TW4")))
fn125tags_TW <- get_FN125_Tags(list(prj_cd=c("LOA_IA25_TW2","LOA_IA25_TW4"))) #empty

#this pulls out the rows with NA values in TLEN
TW_na_rows = fn125_TW[is.na(fn125_TW$TLEN),]

#biodataexport
fn125_TW<-append.spc.names(fn125_TW)
#fn125_TW_edit<-left_join(fn125,fn125tags) <- if there are tags then run this line. 
fn125_TW_edit<-fn125_TW%>%
  dplyr::select(PRJ_CD,
                SAM,
                EFF,
                GRP,
                SPC,
                SPC_NM,
                FISH,
                TLEN,
                RWT,
                AGEST,
                SEX,
                MAT,
                CLIPC) #,TAGID)

lampreyspread_TW<-Lamprey_TW%>%
  group_by(PRJ_CD,SAM,SPC,EFF,GRP,FISH,LAMIJC_TYPE,COMMENT_LAM)%>%
    dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
         Lake="Ontario",
         Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),"R/D"="D","A1-A3" = 0) #no A1/A2/A3 present this year, otherwise A1+A2+A3

biodataexport_TW<-left_join(lampreyspread_TW,fn125_TW_edit)%>%
  dplyr::mutate(EFF=as.numeric(EFF),
                SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),
                Age="",SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                Comments=COMMENT_LAM)

biodataexport_TW$A1 = 0 #add since there weren't occurrences this year
biodataexport_TW$A2 = 0 #add since there weren't occurrences this year
biodataexport_TW$A3 = 0 #add since there weren't occurrences this year
biodataexport_TW$A4 = 0 #add since there weren't occurrences this year
biodataexport_TW$B1 = 0 #add since there weren't occurrences this year
biodataexport_TW$B2 = 0 #add since there weren't occurrences this year
biodataexport_TW$B3 = 0 #add since there weren't occurrences this year
biodataexport_TW$B4 = 0 #add since there weren't occurrences this year

biodataexport_TW_final = biodataexport_TW%>%
  dplyr::select(FishID,
                LiftID,
                Lake,
                Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                "A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)

#write.csv(biodataexport_TW_final,"./data/TW2_TW4/LOA_IA25_TW2_TW4_bio.csv")
write.csv(biodataexport_TW_final,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/TW2_TW4/LOA_IA25_TW2_TW4_bio.csv", 
          row.names = FALSE)

#now do "gear" table export
tem_TW<-fn122_TW%>%
  dplyr::select(PRJ_CD,EFF,SAM,GRTEM0)

fn121_TW<-left_join(fn121_TW,tem_TW)%>%
  dplyr::rename(BottomTempC=GRTEM0)

gearexport_TW<-fn121_TW%>%
  dplyr::filter(SAM%in%lampreyspread_TW$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                SurveyType="Fish Community Assessment",
                SurveyDescription="Fish Community Assessment",
                Gear="Trawl",
                Nights=NA,
                NetMaterial="M",
                NetLengthkm=NA,
                MinMeshmm=12.7,
                MaxMeshmm=101.6,
                MU=NA,Grid=NA)%>%
  dplyr::select(c(LiftID,
                  Lake,
                  Agency,
                  SUBSPACE,
                  DD_LAT0,
                  DD_LON0,
                  MU,Grid,
                  Year,
                  Month,
                  Day,
                  SurveyType,
                  SurveyDescription,
                  Gear,
                  Nights,
                  NetMaterial,
                  MinMeshmm,
                  MaxMeshmm,
                  GRDEPMIN,
                  GRDEPMAX,
                  GRDEPMID,
                  BottomTempC,
                  SITEM0,
                  COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,
                Latitude=DD_LAT0,
                Longitude=DD_LON0,
                Depth1m=GRDEPMIN,
                Depth2m=GRDEPMAX,
                AvgDepthm=GRDEPMID,
                SurfaceTempC=SITEM0,
                Comments=COMMENT1)

#write.csv(gearexport_TW,"./data/TW2_TW4/LOA_IA25_TW2_TW4_gear.csv"))
write.csv(gearexport_TW,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/TW2_TW4/LOA_IA25_TW2_TW4_gear.csv", 
          row.names = FALSE)

##################### TC2 - 2025 ############################ #check if any TLENs are missing - no this year, also no lamprey records for 0R5 creel ##########################
#pull in from glis
Lamprey_TC <- get_SC125_Lamprey(list(prj_cd=c("LOA_SC25_TC2")))
fn121_TC <- get_SC121(list(prj_cd=c("LOA_SC25_TC2")))
fn125_TC <- get_SC125(list(prj_cd=c("LOA_SC25_TC2")))
fn125tags_TC <- get_SC125_Tags(list(prj_cd=c("LOA_SC25_TC2"))) #empty

#this pulls out the rows with NA values in TLEN
TC_na_rows = fn125_TC[is.na(fn125_TC$TLEN),]


#biodataexport
fn125_TC<-append.spc.names(fn125_TC)
#fn125edit<-left_join(fn125,fn125tags)
fn125_TC_edit<-fn125_TC%>%
  dplyr::select(PRJ_CD,SAM,GRP,SPC,SPC_NM,FISH,TLEN,RWT,AGEST,SEX,MAT,CLIPC)#,TAGID)

lampreyspread_TC<-Lamprey_TC%>%
  group_by(PRJ_CD,SAM,SPC,GRP,FISH,LAMIJC_TYPE,COMMENT_LAM)%>%dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
         Lake="Ontario",
         Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),
         "R/D"="D","A1-A3" = A2) #no A1 or A3 present this year, otherwise A1+A2+A3

biodataexport_TC<-left_join(lampreyspread_TC,fn125_TC_edit)%>%
  dplyr::mutate(EFF=NA,SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),
                Age="",SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                Comments=COMMENT_LAM)

biodataexport_TC$CWTAgency = ""
biodataexport_TC$A1 = 0 #add since there weren't occurrences this year
biodataexport_TC$A3 = 0 #add since there weren't occurrences this year
biodataexport_TC$A4 = 0 #add since there weren't occurrences this year
biodataexport_TC$B1 = 0 #add since there weren't occurrences this year
biodataexport_TC$B2 = 0 #add since there weren't occurrences this year
biodataexport_TC$B3 = 0 #add since there weren't occurrences this year

biodataexport_TC_final = biodataexport_TC%>%
  dplyr::select(FishID,
                LiftID,
                Lake,
                Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                CWTAgency,
                "A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)

#write.csv(biodataexport_TC_final,"./data/TC2/LOA_SC25_TC2_bio.csv")
write.csv(biodataexport_TC_final,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/TC2/LOA_SC25_TC2_bio.csv", 
          row.names = FALSE)

#now do "gear" table export
gearexport_TC<-fn121_TC%>%
  dplyr::filter(SAM%in%lampreyspread_TC$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                SurveyType="Angler Survey",
                SurveyDescription="Angler Survey",
                Gear=NA,
                Nights=NA,
                NetMaterial=NA,
                NetLengthkm=NA,
                MinMeshmm=NA,
                MaxMeshmm=NA,
                MU=NA,
                Grid=NA,
                Depth1m=NA,
                Depth2m=NA,
                AvgDepthm=NA,
                BottomTempC=NA,
                SurfaceTempC=NA)%>%
  dplyr::select(c(LiftID,
                  Lake,
                  Agency,
                  SUBSPACE,
                  DD_LAT,
                  DD_LON,
                  MU,
                  Grid,
                  Year,
                  Month,
                  Day,
                  SurveyType,
                  SurveyDescription,
                  Gear,Nights,
                  NetMaterial,
                  MinMeshmm,
                  MaxMeshmm,
                  Depth1m,
                  Depth2m,
                  AvgDepthm,
                  BottomTempC,
                  SurfaceTempC,
                  COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,Latitude=DD_LAT,Longitude=DD_LON,Comments=COMMENT1)


#write.csv(gearexport_TC,"./data//LOA_SC25_TC2_gear.csv")
write.csv(gearexport_TC,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/TC2/LOA_SC25_TC2_gear.csv", 
          row.names = FALSE)

##################### SLR_THI/LSF - 2025 ############################ #check if any TLENs are missing - yes #####################################################
#pull in from glis
Lamprey_SLR <- get_FN125_Lamprey(list(prj_cd=c("SLR_IA25_THI","SLR_IA25_LSF")))
fn121_SLR <- get_FN121(list(prj_cd=c("SLR_IA25_THI","SLR_IA25_LSF")))
fn122_SLR <- get_FN122(list(prj_cd=c("SLR_IA25_THI","SLR_IA25_LSF")))
fn125_SLR <- get_FN125(list(prj_cd=c("SLR_IA25_THI","SLR_IA25_LSF")))
fn125tags_SLR <- get_FN125_Tags(list(prj_cd=c("SLR_IA25_THI","SLR_IA25_LSF")))

#this pulls out the rows with NA values in TLEN 
SLR_na_rows = fn125_SLR[is.na(fn125_SLR$TLEN),] #24 missing TLEN

#this completely removes the NA values in TLEN 
fn125_SLR_noNA = fn125_SLR %>%
  filter(!is.na(TLEN))

#If there are NA values then convert FLEN to TLEN
fn125.convert_SLR = estimate_tlen2(SLR_na_rows, results)

#rbind the converted table with the larger data set
updated_fn125_SLR = bind_rows(fn125_SLR_noNA, fn125.convert_SLR)
##there are six rows that still have NA because there was no FLEN to convert (2025)

#biodataexport
updated_fn125_SLR<-append.spc.names(updated_fn125_SLR) #if no conversion needed, this would just be fn125
fn125_SLR_edit<-left_join(updated_fn125_SLR,fn125tags_SLR)

fn125_SLR_edit<-fn125_SLR_edit%>%
  dplyr::select(PRJ_CD,SAM,EFF,GRP,SPC,SPC_NM,FISH,TLEN,RWT,AGEST,SEX,MAT,CLIPC,TAGID) 

lampreyspread_SLR<-Lamprey_SLR%>%
  group_by(PRJ_CD,
           SAM,
           SPC,
           EFF,
           GRP,
           FISH,
           LAMIJC_TYPE,
           COMMENT_LAM)%>%
  dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
         Lake="Ontario",
         Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),
         "R/D"="D","A1-A3" = 0) #no A1/A2/A3 present this year, otherwise A1+A2+A3

biodataexport_SLR<-left_join(lampreyspread_SLR,fn125_SLR_edit)%>%
  dplyr::mutate(EFF=as.numeric(EFF),
                SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),
                Age="",SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                CWTAgency=TAGID,
                Comments=COMMENT_LAM) 

biodataexport_SLR$A1 = 0 #add since there weren't occurrences this year
biodataexport_SLR$A2 = 0 #add since there weren't occurrences this year
biodataexport_SLR$A3 = 0 #add since there weren't occurrences this year
biodataexport_SLR$A4 = 0 #add since there weren't occurrences this year
biodataexport_SLR$B1 = 0 #add since there weren't occurrences this year
biodataexport_SLR$B2 = 0 #add since there weren't occurrences this year
biodataexport_SLR$B3 = 0 #add since there weren't occurrences this year
biodataexport_SLR$B4 = 0 #add since there weren't occurrences this year

biodataexport_SLR_final = biodataexport_SLR%>%
  dplyr::select(FishID,
                LiftID,
                Lake,
                Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                CWTAgency,
                "A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)

#write.csv(biodataexport_SLR_final,"./data/THI_LSF/SLR_IA25_THI_LSF_bio.csv")
write.csv(biodataexport_SLR_final,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/THI_LSF/LOA_IA25_THI_LSF_bio.csv", 
          row.names = FALSE)


#now do "gear" table export
tem_SLR<-fn122_SLR%>%
  dplyr::select(PRJ_CD,EFF,SAM,GRTEM0)

fn121_SLR<-left_join(fn121_SLR,tem_SLR)%>%
  dplyr::rename(BottomTempC=GRTEM0)

gearexport_SLR<-fn121_SLR%>%
  dplyr::filter(SAM%in%lampreyspread_SLR$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                SurveyType="Fish Community Assessment",
                SurveyDescription="Fish Community Assessment",
                Gear="Gillnet",
                Nights=1,
                NetMaterial="M",
                NetLengthkm=0.061,
                MinMeshmm=38.1,
                MaxMeshmm=152.4,
                MU=NA,
                Grid=NA)%>%
  dplyr::select(c(LiftID,
                  Lake,
                  Agency,
                  SUBSPACE,
                  DD_LAT0,
                  DD_LON0,
                  MU,
                  Grid,
                  Year,
                  Month,
                  Day,
                  SurveyType,
                  SurveyDescription,
                  Gear,
                  Nights,
                  NetMaterial,
                  MinMeshmm,
                  MaxMeshmm,
                  GRDEPMIN,
                  GRDEPMAX,
                  GRDEPMID,
                  BottomTempC,
                  SITEM0,
                  COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,
                Latitude=DD_LAT0,
                Longitude=DD_LON0,
                Depth1m=GRDEPMIN,
                Depth2m=GRDEPMAX,
                AvgDepthm=GRDEPMID,
                SurfaceTempC=SITEM0,
                Comments=COMMENT1)


#write.csv(gearexport_SLR,"./data/THI_LSF/SLR_IA25_THI_LSF_gear.csv")
write.csv(gearexport_SLR,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/THI_LSF/LOA_IA25_THI_LSF_gear.csv", 
          row.names = FALSE)


##################### CF25_003 - 2025 ############################ #pulled from access - doesn't get put on GLIS ##########################################################
library(RODBC)

#the connection to the access database will change depending where you store it. I currently have it in my files. 
LakeTrout_Bycatch_db = "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/CF_003/LOA_IA25_CF25_003.accdb"
fetch_data <- function(table, prj_cd, LakeTrout_Bycatch_db){
  DBConnection <- odbcConnectAccess2007(LakeTrout_Bycatch_db, uid = "", pwd = "")
  dat <- sqlFetch(DBConnection, table, as.is=TRUE, stringsAsFactors=FALSE)
  odbcClose(DBConnection)
  return(dat)
}

#load all tables
Lamprey_LakeTrout <- fetch_data("FN125_lamprey", params$prj_cd, LakeTrout_Bycatch_db)
fn121_LakeTrout <- fetch_data("FN121", params$prj_cd, LakeTrout_Bycatch_db)
fn125_LakeTrout <- fetch_data("FN125", params$prj_cd, LakeTrout_Bycatch_db)
fn125tags_LakeTrout <- fetch_data("FN125_tags", params$prj_cd, LakeTrout_Bycatch_db) #empty

#this pulls out the rows with NA values in TLEN
LakeTrout_na_rows = fn125_LakeTrout[is.na(fn125_LakeTrout$TLEN),]
#no NA values present.


#biodataexport
fn125_LakeTrout<-append.spc.names(fn125_LakeTrout)
#fn125edit<-left_join(fn125,fn125tags) #no tags present in 2025
fn125_LakeTrout_edit<-fn125_LakeTrout%>%
  dplyr::select(
    PRJ_CD,
    SAM,
    EFF,
    GRP,
    SPC,
    SPC_NM,
    FISH,
    TLEN,
    RWT,
    AGEST,
    SEX,
    MAT,
    CLIPC) #TAGID - add if tags are recorded

lampreyspread_LakeTrout<-Lamprey_LakeTrout%>%
  group_by(PRJ_CD,
           SAM,
           SPC,
           EFF,
           GRP,
           FISH,
           LAMIJC_TYPE,
           COMMENT_LAM)%>%
  dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
         Lake="Ontario",
         Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),
         "R/D"="D",
         "A1-A3" = A1+A2+A3)

biodataexport_LakeTrout<-left_join(lampreyspread_LakeTrout,fn125_LakeTrout_edit)%>%
  dplyr::mutate(EFF=as.numeric(EFF),
                SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),
                Age="",SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                Comments=COMMENT_LAM) #CWTAgency=TAGID, add if tags recorded


biodataexport_LakeTrout$B1 = 0 #add since there weren't occurrences this year


biodataexport_LakeTrout_final = biodataexport_LakeTrout%>%
  dplyr::select(FishID,
                LiftID,
                Lake,
                Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                "A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)

#write.csv(biodataexport_LakeTrout_final,"./data/LOA_CF25_003_bio.csv")
write.csv(biodataexport_LakeTrout_final,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/CF_003/LOA_CF25_003_bio.csv", 
          row.names = FALSE)

#now do "gear" table export
gearexport_LakeTrout<-fn121_LakeTrout%>%
  dplyr::filter(SAM%in%lampreyspread_LakeTrout$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(COLLECTION_DATE, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(COLLECTION_DATE, format = "%Y-%m-%d")),
                SurveyType="Commercial Bycatch Sampling",
                SurveyDescription="Commercial Bycatch Sampling",
                Gear="Commercial Fishery (Trapnet/Gillnet)",
                Nights=NA,
                NetMaterial=NA,
                NetLengthkm=NA,
                MinMeshmm=NA,
                MaxMeshmm=NA,
                MU=NA,
                Grid=NA,
                Depth1m=NA,
                Depth2m=NA,
                AvgDepthm=NA,
                BottomTempC=NA,
                SurfaceTempC=NA, 
                COMMENT1=NA)%>%
  dplyr::select(c(LiftID,
                  Lake,
                  Agency,
                  SUBSPACE,
                  DD_LAT,
                  DD_LON,
                  MU,
                  Grid,
                  Year,
                  Month,
                  Day,
                  SurveyType,
                  SurveyDescription,
                  Gear,
                  Nights,
                  NetMaterial,
                  MinMeshmm,
                  MaxMeshmm,
                  Depth1m,
                  Depth2m,
                  AvgDepthm,
                  BottomTempC,
                  SurfaceTempC,
                  COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,Latitude=DD_LAT,Longitude=DD_LON,Comments=COMMENT1)


#write.csv(gearexport_LakeTrout,"./data//LOA_CF25_003_gear.csv")
write.csv(gearexport_LakeTrout,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/CF_003/LOA_CF25_003_gear.csv", 
          row.names = FALSE)

##################### (need to run still)CF25_001 - 2025 ############################ #pulled from access - doesn't get put on GLIS  - ##########################################################

library(RODBC)

#the connection to the access database will change depending where you store it. I currently have it in my files. 
CommCatch_db = "C:/Users/StinsoM/Documents/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/CF_001/LOA_CF25_001.accdb"
fetch_data <- function(table, prj_cd, CommCatch_db){
  DBConnection <- odbcConnectAccess2007(CommCatch_db, uid = "", pwd = "")
  dat <- sqlFetch(DBConnection, table, as.is=TRUE, stringsAsFactors=FALSE)
  odbcClose(DBConnection)
  return(dat)
}

#load all tables
Lamprey_CommCatch <- fetch_data("FN125_lamprey", params$prj_cd, CommCatch_db)
fn121_CommCatch <- fetch_data("FN121", params$prj_cd, CommCatch_db)
fn125_CommCatch <- fetch_data("FN125", params$prj_cd, CommCatch_db)
fn125tags_CommCatch <- fetch_data("FN125_tags", params$prj_cd, CommCatch_db) 

#this pulls out the rows with NA values in TLEN
CommCatch_na_rows = fn125_CommCatch[is.na(fn125_CommCatch$TLEN),]



#biodataexport
fn125_CommCatch<-append.spc.names(fn125_CommCatch)
#fn125edit<-left_join(fn125,fn125tags) #no tags present in 2025
fn125_CommCatch_edit<-fn125_CommCatch%>%
  dplyr::select(
    PRJ_CD,
    SAM,
    EFF,
    GRP,
    SPC,
    SPC_NM,
    FISH,
    TLEN,
    RWT,
    AGEST,
    SEX,
    MAT,
    CLIPC) #TAGID - add if tags are recorded

lampreyspread_CommCatch<-Lamprey_CommCatch%>%
  group_by(PRJ_CD,
           SAM,
           SPC,
           EFF,
           GRP,
           FISH,
           LAMIJC_TYPE,
           COMMENT_LAM)%>%
  dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
         Lake="Ontario",
         Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),
         "R/D"="D",
         "A1-A3" = A1+A2+A3) ##check

biodataexport_CommCatch<-left_join(lampreyspread_CommCatch,fn125_CommCatch_edit)%>%
  dplyr::mutate(EFF=as.numeric(EFF),
                SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),
                Age="",SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                Comments=COMMENT_LAM) #CWTAgency=TAGID, add if tags recorded


biodataexport_CommCatch$B1 = 0 #add since there weren't occurrences this year


biodataexport_CommCatch_final = biodataexport_CommCatch%>%
  dplyr::select(FishID,
                LiftID,
                Lake,
                Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                "A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)

#write.csv(biodataexport_CommCatch_final,"./data/LOA_CF25_001_bio.csv")
write.csv(biodataexport_CommCatch_final,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/CF_001/LOA_CF25_001_bio.csv", 
          row.names = FALSE)

#now do "gear" table export
gearexport_CommCatch<-fn121_CommCatch%>%
  dplyr::filter(SAM%in%lampreyspread_CommCatch$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(COLLECTION_DATE, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(COLLECTION_DATE, format = "%Y-%m-%d")),
                SurveyType="Commercial Harvest Sampling",
                SurveyDescription="Commercial Harvest Sampling",
                Gear="Commercial Fishery (Trapnet/Gillnet)",
                Nights=NA,
                NetMaterial=NA,
                NetLengthkm=NA,
                MinMeshmm=NA,
                MaxMeshmm=NA,
                MU=NA,
                Grid=NA,
                Depth1m=NA,
                Depth2m=NA,
                AvgDepthm=NA,
                BottomTempC=NA,
                SurfaceTempC=NA, 
                COMMENT1=NA)%>%
  dplyr::select(c(LiftID,
                  Lake,
                  Agency,
                  SUBSPACE,
                  DD_LAT,
                  DD_LON,
                  MU,
                  Grid,
                  Year,
                  Month,
                  Day,
                  SurveyType,
                  SurveyDescription,
                  Gear,
                  Nights,
                  NetMaterial,
                  MinMeshmm,
                  MaxMeshmm,
                  Depth1m,
                  Depth2m,
                  AvgDepthm,
                  BottomTempC,
                  SurfaceTempC,
                  COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,Latitude=DD_LAT,Longitude=DD_LON,Comments=COMMENT1)

#write.csv(gearexport_CommCatch,"./data//LOA_CF25_001_gear.csv")
write.csv(gearexport_CommCatch,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/CF_001/LOA_CF25_001_gear.csv", 
          row.names = FALSE)

##################### (need to run still)SC25_002 - 2025 ############################ #check if any TLENs are missing - ##########################################################
#pull in from glis
Lamprey_Creel <- get_SC125_Lamprey(list(prj_cd=c("LOA_SC25_002")))
fn121_Creel <- get_SC121(list(prj_cd=c("LOA_SC25_002")))
fn125_Creel <- get_SC125(list(prj_cd=c("LOA_SC25_002")))
fn125tags_Creel <- get_SC125_Tags(list(prj_cd=c("LOA_SC25_002"))) 

#this pulls out the rows with NA values in TLEN
Creel_na_rows = fn125_Creel[is.na(fn125_Creel$TLEN),]


#biodataexport
fn125_Creel<-append.spc.names(fn125_Creel)
#fn125edit<-left_join(fn125,fn125tags)
fn125_Creel_edit<-fn125_Creel%>%
  dplyr::select(PRJ_CD,SAM,GRP,SPC,SPC_NM,FISH,TLEN,RWT,AGEST,SEX,MAT,CLIPC)#,TAGID)

lampreyspread_Creel<-Lamprey_Creel%>%
  group_by(PRJ_CD,SAM,SPC,GRP,FISH,LAMIJC_TYPE,COMMENT_LAM)%>%
  dplyr::summarize(Count=length(LAMIJC_TYPE))%>%
  spread(LAMIJC_TYPE,Count,fill=0)%>%ungroup()%>%
  mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
         Lake="Ontario",
         Agency="OMNR",
         FishID = paste(LiftID,SPC,FISH,sep="_"),
         "R/D"="D","A1-A3" = A2) #no A1 or A3 present this year, otherwise A1+A2+A3

biodataexport_Creel<-left_join(lampreyspread_Creel,fn125_Creel_edit)%>%
  dplyr::mutate(EFF=NA,SPC=as.numeric(SPC),
                SEX=ifelse(SEX==1,"M","F"),
                SEX=ifelse(SEX==9,"U",SEX),
                Age="",SpeciesAbbrev="")%>%
  dplyr::select(!c(PRJ_CD,SAM,FISH,GRP,`0`))%>%
  dplyr::rename(MaturityAgency=MAT,
                MeshSizemm=EFF,
                SpeciesNumber=SPC,
                SpeciesName=SPC_NM,
                "Length(mm)"=TLEN,
                "Weight(g)"=RWT,
                FinClipAgency=CLIPC,
                AgeStructure=AGEST,
                SexAgency=SEX,
                Comments=COMMENT_LAM)

biodataexport_Creel$CWTAgency = ""
biodataexport_Creel$A1 = 0 #add since there weren't occurrences this year
biodataexport_Creel$A3 = 0 #add since there weren't occurrences this year
biodataexport_Creel$A4 = 0 #add since there weren't occurrences this year
biodataexport_Creel$B1 = 0 #add since there weren't occurrences this year
biodataexport_Creel$B2 = 0 #add since there weren't occurrences this year
biodataexport_Creel$B3 = 0 #add since there weren't occurrences this year

biodataexport_Creel_final = biodataexport_Creel%>%
  dplyr::select(FishID,
                LiftID,
                Lake,
                Agency,
                SpeciesName,
                SpeciesNumber,
                SpeciesAbbrev,
                MeshSizemm,
                "Length(mm)",
                "Weight(g)",
                "R/D",
                Age,
                AgeStructure,
                SexAgency,
                MaturityAgency,
                FinClipAgency,
                CWTAgency,
                "A1-A3",A1,A2,A3,A4,B1,B2,B3,B4)

#write.csv(biodataexport_Creel_final,"./data/SC25_002/LOA_SC25_002_bio.csv")
write.csv(biodataexport_Creel_final,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/SC25_002/LOA_SC25_002_bio.csv", 
          row.names = FALSE)

#now do "gear" table export
gearexport_Creel<-fn121_Creel%>%
  dplyr::filter(SAM%in%lampreyspread_Creel$SAM)%>%
  dplyr::mutate(LiftID=paste("OMNR",PRJ_CD,SAM,sep="_"),
                Lake="Ontario",
                Agency="OMNR",
                Year=2025,
                Month=month(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                Day=day(as.POSIXlt(EFFDT0, format = "%Y-%m-%d")),
                SurveyType="Angler Survey",
                SurveyDescription="Angler Survey",
                Gear=NA,
                Nights=NA,
                NetMaterial=NA,
                NetLengthkm=NA,
                MinMeshmm=NA,
                MaxMeshmm=NA,
                MU=NA,
                Grid=NA,
                Depth1m=NA,
                Depth2m=NA,
                AvgDepthm=NA,
                BottomTempC=NA,
                SurfaceTempC=NA)%>%
  dplyr::select(c(LiftID,
                  Lake,
                  Agency,
                  SUBSPACE,
                  DD_LAT,
                  DD_LON,
                  MU,
                  Grid,
                  Year,
                  Month,
                  Day,
                  SurveyType,
                  SurveyDescription,
                  Gear,Nights,
                  NetMaterial,
                  MinMeshmm,
                  MaxMeshmm,
                  Depth1m,
                  Depth2m,
                  AvgDepthm,
                  BottomTempC,
                  SurfaceTempC,
                  COMMENT1))%>%
  dplyr::rename(Location=SUBSPACE,Latitude=DD_LAT,Longitude=DD_LON,Comments=COMMENT1)


#write.csv(gearexport_Creel,"./data//LOA_SC25_002_gear.csv")
write.csv(gearexport_Creel,"~/Program Folders (R etc)/Sea Lampray/SeaLampray_2025/SC25_002/LOA_SC25_002_gear.csv", 
          row.names = FALSE)




#####################################################################################################################################################



bioexportcompile = rbind(biodataexport1,
                         biodataexport_trawl2,
                         biodataexport_NSCIN_final,
                         biodataexport_TW_final,
                         biodataexport_TC_final,
                         biodataexport_SLR_final,
                         biodataexport_LakeTrout_final,
                         biodataexport_CommCatch_final,
                         biodataexport_Creel_final)
gearexportcompile = rbind(gearexport1,
                          gearexport_trawl,
                          gearexport_NSCIN,
                          gearexport_TW,
                          gearexport_TC,
                          gearexport_SLR,
                          gearexport_LakeTrout,
                          gearexport_CommCatch,
                          gearexport_Creel)

write_xlsx(bioexportcompile,"./data/exports/MNR_biodata.xlsx")
write_xlsx(gearexportcompile,"./data/exports/MNR_gear.xlsx")

#alter maturityagency to character codes from GLIS codes (1=I=immature, 2=M=mature, 9=U=unknown)
#manually combined the biodata and gear into a single file with two sheets
#manually combined tom's fishway/non-GLIS projects
#added a codereference sheet for maturity/age structure if needed
