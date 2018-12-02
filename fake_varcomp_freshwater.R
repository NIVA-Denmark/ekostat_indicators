df0 <- distinct(df.indicators,Measurement,Parameter)

df1 <- df.variances %>% 
  filter(Water_type=="Coast") %>%
  group_by(Water_type,Measurement) %>%
  summarise(V_station=mean(V_station,na.rm=T),  
            V_obspoint=mean(V_obspoint,na.rm=T), 
            V_year=mean(V_year,na.rm=T), 
            V_yearmonth=mean(V_yearmonth,na.rm=T), 
            V_monthhour=mean(V_monthhour,na.rm=T), 
            V_yearmonthhour=mean(V_yearmonthhour,na.rm=T), 
            V_stationdate=mean(V_stationdate,na.rm=T), 
            V_stationyear=mean(V_stationyear,na.rm=T), 
            V_stationmonth=mean(V_stationmonth,na.rm=T), 
            V_institution=mean(V_institution,na.rm=T), 
            V_replication=mean(V_replication,na.rm=T)) %>%
  ungroup() %>%
  left_join(df0,
            by=c("Measurement"))

df2 <- dfc %>% 
  summarise(secchi=mean(secchi,na.rm=T), 
            DIN=mean(DIN,na.rm=T), 
            DIP=mean(DIP,na.rm=T), 
            TN=mean(TN,na.rm=T), 
            TP=mean(TP,na.rm=T), 
            chla=mean(chla,na.rm=T), 
            biovol=mean(biovol,na.rm=T), 
            O2=mean(O2,na.rm=T), 
            BQI=mean(BQI,na.rm=T), 
            MSMDI=mean(MSMDI,na.rm=T)) %>%
  gather(key="Parameter",value="mean")

df3 <- df1 %>% left_join(df2,by=c("Parameter")) %>%
  mutate(V_station=ifelse(is.nan(V_station),NA,V_station/mean),  
         V_obspoint=ifelse(is.nan(V_obspoint),NA,V_obspoint/mean), 
         V_year=ifelse(is.nan(V_year),NA,V_year/mean), 
         V_yearmonth=ifelse(is.nan(V_yearmonth),NA,V_yearmonth/mean), 
         V_monthhour=ifelse(is.nan(V_monthhour),NA,V_monthhour/mean), 
         V_yearmonthhour=ifelse(is.nan(V_yearmonthhour),NA,V_yearmonthhour/mean), 
         V_stationdate=ifelse(is.nan(V_stationdate),NA,V_stationdate/mean), 
         V_stationyear=ifelse(is.nan(V_stationyear),NA,V_stationyear/mean), 
         V_stationmonth=ifelse(is.nan(V_stationmonth),NA,V_stationmonth/mean), 
         V_institution=ifelse(is.nan(V_institution),NA,V_institution/mean), 
         V_replication=ifelse(is.nan(V_replication),NA,V_replication/mean))

df4 <- df3 %>% summarise(V_station=mean(V_station,na.rm=T),  
                         V_obspoint=mean(V_obspoint,na.rm=T), 
                         V_year=mean(V_year,na.rm=T), 
                         V_yearmonth=mean(V_yearmonth,na.rm=T), 
                         V_monthhour=mean(V_monthhour,na.rm=T), 
                         V_yearmonthhour=mean(V_yearmonthhour,na.rm=T), 
                         V_stationdate=mean(V_stationdate,na.rm=T), 
                         V_stationyear=mean(V_stationyear,na.rm=T), 
                         V_stationmonth=mean(V_stationmonth,na.rm=T), 
                         V_institution=mean(V_institution,na.rm=T), 
                         V_replication=mean(V_replication,na.rm=T)) %>%
  mutate(X=1)

dflmean<-dfl %>% 
  summarise(
    Nspecies_phytoplankton=mean(Nspecies_phytoplankton,na.rm=T), 
    biovol=mean(biovol,na.rm=T), 
    Proportion_cyanobacteria=mean(Proportion_cyanobacteria,na.rm=T), 
    TrophicPlanktonIndex=mean(TrophicPlanktonIndex,na.rm=T), 
    TrophicMacrophyteIndex=mean(TrophicMacrophyteIndex,na.rm=T), 
    BenthicInvertebratesASPT=mean(BenthicInvertebratesASPT,na.rm=T), 
    BenthicInvertebratesMILA=mean(BenthicInvertebratesMILA,na.rm=T), 
    BenthicInvertebratesBQI=mean(BenthicInvertebratesBQI,na.rm=T), 
    chla=mean(chla,na.rm=T), 
    EQR8=mean(EQR8,na.rm=T),
    AindexW5=mean(AindexW5,na.rm=T), 
    EindexW3=mean(EindexW3,na.rm=T), 
    SecchiDepth=mean(SecchiDepth,na.rm=T), 
    Oxygen=mean(Oxygen,na.rm=T), 
    TP=mean(TP,na.rm=T))%>%
  gather(key="Parameter",value="mean") %>%
  mutate(Water_type="Lake") %>%
  mutate(mean=ifelse(is.nan(mean),NA,mean*ifelse(mean<0,-1,1))) %>%
  left_join(df0,by=c("Parameter"))

dfrmean<-dfr %>% 
  summarise(
    BenthicDiatomsIPS=mean(BenthicDiatomsIPS,na.rm=T), 
    BenthicDiatomsACID=mean(BenthicDiatomsACID,na.rm=T), 
    BenthicDiatomsPctPT=mean(BenthicDiatomsPctPT,na.rm=T), 
    BenthicDiatomsTDI=mean(BenthicDiatomsTDI,na.rm=T), 
    BenthicInvertebratesASPT=mean(BenthicInvertebratesASPT,na.rm=T), 
    BenthicInvertebratesDJ=mean(BenthicInvertebratesDJ,na.rm=T), 
    BenthicInvertebratesMISA=mean(BenthicInvertebratesMISA,na.rm=T), 
    VIXsm=mean(VIXsm,na.rm=T), 
    VIXh=mean(VIXh,na.rm=T), 
    VIX=mean(VIX,na.rm=T), 
    SecchiDepth=mean(SecchiDepth,na.rm=T), 
    Oxygen=mean(Oxygen,na.rm=T), 
    TP=mean(TP,na.rm=T)) %>%
  gather(key="Parameter",value="mean") %>%
  mutate(Water_type="River")%>%
  mutate(mean=ifelse(is.nan(mean),NA,mean*ifelse(mean<0,-1,1))) %>%
  left_join(df0,by=c("Parameter"))

dfmean <- bind_rows(dflmean,dfrmean) %>% 
  select(Water_type,Measurement,mean)

dfvarnew<- df.variances %>%
  filter(Water_type %in% c("Lake","River")) %>%
  select(Water_type,Measurement,Unit,Transform,Region,Type,Typename) %>%
  left_join(dfmean,by=c("Water_type","Measurement")) %>%
  mutate(X=1) %>%
  left_join(df4,by=c("X")) %>%
  mutate(V_station=mean*V_station,  
         V_obspoint=mean*V_obspoint, 
         V_year=mean*V_year, 
         V_yearmonth=mean*V_yearmonth, 
         V_monthhour=mean*V_monthhour, 
         V_yearmonthhour=mean*V_yearmonthhour, 
         V_stationdate=mean*V_stationdate, 
         V_stationyear=mean*V_stationyear, 
         V_stationmonth=mean*V_stationmonth, 
         V_institution=mean*V_institution, 
         V_replication=mean*V_replication) %>%
  select(-c(X,mean))

write.table(dfvarnew,file="varnew.txt",row.names=F,col.names=T,sep="\t",quote=F,na="")
