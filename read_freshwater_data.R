
# ------- Get type information and other properties to be used in filtering/selecting WBs -----------
df_wb_sjo <- read.table("data/waterbodies_sjo_old.txt",fileEncoding = "UTF-8",
                        sep = "\t",stringsAsFactors = F,header = T, 
                        quote = "\"")

#df_wb_sjo <- df_wb_sjo[,c("Vatten.ID","Namn.Vatten","Huvudavrinningsområde","Län","Kommun.er.","Distrikt","Vattenkategori","Vattentyp...Sjö","Limnisk.vattentypsregion","Medeldjup..m.","Alkalinitet..mekv.l.","Humus..mg.Pt.l.")]
df_wb_sjo <- df_wb_sjo[,c("EU_CD.Vatten","Namn.Vatten","Huvudavrinningsområde","Län","Kommun.er.","Distrikt","Vattenkategori","Vattentyp...Sjö","Limnisk.ekoregion.Kustvattentyp","Djupkategori","Bakgrundsalkalinitet","Färg..Humus.")]
names(df_wb_sjo)<-c("WB_ID","WB_Name","Catchment","Län","Municipality","District","Category","typology","typology_varcomp","TypeDepth","TypeAlkalinity","TypeHumus")
df_wb_sjo$CLR<-"Lake"

df_wb_vdr <- read.table("data/waterbodies_vdr.txt",fileEncoding = "UTF-8",
                        sep = "\t",stringsAsFactors = F,header = T, 
                        quote = "\"")

df_wb_vdr <- df_wb_vdr[,c("EU_CD.Vatten","Namn.Vatten","Huvudavrinningsområde","Län","Kommun.er.","Distrikt","Vattenkategori","Vattentyp...Vattendrag","Limnisk.ekoregion.Kustvattentyp","Avrinningsområde","Färg..Humus.","Bakgrundsalkalinitet")]
names(df_wb_vdr) <- c("WB_ID","WB_Name","Catchment","Län","Municipality","District","Category","typology","typology_varcomp","DrainageBasin","ColourHumus","BackgroundAlkalinity")
df_wb_vdr$CLR<-"River"

# list of freshwater WBs with typology information
df_wb <- bind_rows(df_wb_sjo,df_wb_vdr)

# -----------------------------------------------------------------------
# function to remove attributes carried from SAS data files
stripattr<-function(df,attrlist){
  for (attrname in attrlist){
    for (var in colnames(df)) {
      attr(df[[deparse(as.name(var))]], attrname) <- NULL
      }
    attr(df, attrname) <- NULL
    }
  return(df)
}
# -----------------------------------------------------------------------
# load data files
dflbio <- read_sas("data/lake_biolindex.sas7bdat")  # Lake biology data
dflwq <- read_sas("data/lake_WQdata.sas7bdat")      # Lake WQ data
dfrbio <- read_sas("data/river_biolindex.sas7bdat") # River biology data
dfrwq <- read_sas("data/river_WQdata.sas7bdat")     # River WQ data

#filter out older data
dflbio <- dflbio %>% filter(year>2006)
dflwq <- dflwq %>% filter(year>2006)
dfrbio <-dfrbio %>% filter(year>2006)
dfrwq <- dfrwq %>% filter(year>2006)

# remove SAS label and format attributes
dflbio <- stripattr(dflbio,c("label","format.sas"))
dflwq <- stripattr(dflwq,c("label","format.sas"))
dfrbio <- stripattr(dfrbio,c("label","format.sas"))
dfrwq <- stripattr(dfrwq,c("label","format.sas"))

# select and merge lake data
dflbio <- dflbio %>% select(WB_ID,station,obspoint,Nspecies_phytoplankton,biovol,
                            Proportion_cyanobacteria,TrophicPlanktonIndex,
                            TrophicMacrophyteIndex,BenthicInvertebratesASPT,
                            BenthicInvertebratesMILA,BenthicInvertebratesBQI,
                            institution,date,water_type,Altitude,LakeArea,
                            LakeDepth,EQR8,AindexW5,EindexW3,year,month,date,sali)
dflwq <- dflwq %>% select(WB_ID,station,obspoint,institution,chla,SecchiDepth,O2_surf,O2_bot,TP,year,month,date) # Oxygen

dfl <- bind_rows(dflbio,dflwq)

# select and merge river data
dfrbio <- dfrbio %>% select(WB_ID,station,obspoint,BenthicDiatomsIPS,BenthicDiatomsACID,
                            BenthicDiatomsPctPT,BenthicDiatomsTDI,
                            BenthicInvertebratesASPT,BenthicInvertebratesDJ,
                            BenthicInvertebratesMISA,VIXsm,VIXh,VIX,institution,year,month,date,sali)
dfrwq <- dfrwq %>% select(WB_ID,station,institution,obspoint,SecchiDepth,TP,year,month,date) # O2_surf, O2_bot
dfr <- bind_rows(dfrbio,dfrwq) 

# join with typology information from the WB list
dfr <- dfr %>% 
  left_join(select(df_wb,WB_ID,typology,typology_varcomp,CLR),by = "WB_ID") %>%
  rename(WB=WB_ID) %>% 
  filter(WB!="") %>% 
  filter(typology!="") %>% 
  mutate(typology_varcomp = paste0("S",as.character(typology_varcomp))) %>%
  mutate(period=ifelse(year<2007,NA,ifelse(year<2013,"2007-2012","2013-2018"))) %>% 
  filter(!is.na(period)) %>%
  mutate(obspoint=ifelse(obspoint=="",station,obspoint))
  
dfl <- dfl %>% 
  left_join(select(df_wb,WB_ID,typology,typology_varcomp,CLR),by = "WB_ID") %>%
  rename(WB=WB_ID) %>% 
  filter(WB!="") %>% 
  filter(typology!="") %>% 
  mutate(typology_varcomp = paste0("S",as.character(typology_varcomp))) %>%
  mutate(period=ifelse(year<2007,NA,ifelse(year<2013,"2007-2012","2013-2018"))) %>% 
  filter(!is.na(period)) %>%
  mutate(obspoint=ifelse(obspoint=="",station,obspoint))





