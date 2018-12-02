library(tidyverse)
library(RSQLite)
rm(list=ls())
#LanID,Lan_name,WB_ID,Name,District,Type

# ------- Get type information and other properties to be used in filtering/selecting WBs -----------

# Lakes
df_wb_sjo <- read.table("data/waterbodies_sjo_old.txt",fileEncoding = "UTF-8",
                        sep = "\t",stringsAsFactors = F,header = T, 
                        quote = "\"")

df_wb_sjo <- df_wb_sjo[,c("EU_CD.Vatten","Namn.Vatten","Huvudavrinningsområde","Län","Kommun.er.","Distrikt","Vattenkategori","Vattentyp...Sjö","Limnisk.ekoregion.Kustvattentyp","Djupkategori","Bakgrundsalkalinitet","Färg..Humus.")]
names(df_wb_sjo)<-c("WB_ID","WB_Name","Catchment","Lan","Municipality","District","Category","typology","typology_varcomp","TypeDepth","TypeAlkalinity","TypeHumus")
df_wb_sjo$CLR<-"Lake"

# Rivers
df_wb_vdr <- read.table("data/waterbodies_vdr.txt",fileEncoding = "UTF-8",
                        sep = "\t",stringsAsFactors = F,header = T, 
                        quote = "\"")

df_wb_vdr <- df_wb_vdr[,c("EU_CD.Vatten","Namn.Vatten","Huvudavrinningsområde","Län","Kommun.er.","Distrikt","Vattenkategori","Vattentyp...Vattendrag","Limnisk.ekoregion.Kustvattentyp","Avrinningsområde","Färg..Humus.","Bakgrundsalkalinitet")]
names(df_wb_vdr) <- c("WB_ID","WB_Name","Catchment","Lan","Municipality","District","Category","typology","typology_varcomp","DrainageBasin","ColourHumus","BackgroundAlkalinity")
df_wb_vdr$CLR<-"River"


# Coast
df_wb_c <- read.table("data/waterbodies.txt",fileEncoding = "UTF-8",
                    sep = "\t",stringsAsFactors = F,header = T)
df_wb_c <- df_wb_c %>% 
  rename(WB_ID=WaterbodyID,
         WB_Name=WaterbodyName,
         Lan=County,
         Municipality=Municipalities,
         typology=TypeID,
         District=DistrictID) %>%
  mutate(CLR="Coast")

# ---- list of WBs with typology information
df_wb <- bind_rows(df_wb_sjo,df_wb_vdr,df_wb_c)

# ---- Län table -------------------------------------
df_lan <- df_wb %>% 
  distinct(WB_ID,Lan) 

colsLan<-c("Lan1", "Lan2", "Lan3", "Lan4")
df_lan <- df_lan %>%
  separate(col=Lan,into=colsLan,sep = ",") %>%
  gather(key="n",value="Lan",colsLan) %>%
  select(-n) %>%
  filter(!is.na(Lan))

df_lan <- df_lan %>%
  separate(col=Lan,into = c("LanName","LanID"),sep = " - ",remove=F) %>%
  arrange(LanName)

# ---- Municipality table -------------------------------------

df_mun <- df_wb %>% 
  distinct(WB_ID,Municipality) 

colsMun<-c("Mun1", "Mun2", "Mun3","Mun4", "Mun5", "Mun6","Mun7", "Mun8", "Mun9")
df_mun <- df_mun %>%
  separate(col=Municipality, into=colsMun, sep=",") %>%
  gather(key="n",value="Mun",colsMun) %>%
  select(-n) %>%
  filter(!is.na(Mun))

df_mun <- df_mun %>%
  separate(col=Mun,into = c("MunName","MunID"),sep = " - ",remove=F) %>%
  arrange(MunName)

# ----- WB info w/o Län, Municipality --------------------------------------

df_wb <- df_wb %>% 
  arrange(WB_Name)
  #select(-c(Lan,Municipality)) %>%


# ----- write to DB --------------------------------------

dbpath<-"output/ekostat_coast.db"
#dbpath<-"../efs/ekostat/ekostat.db"

db <- dbConnect(SQLite(), dbname=dbpath)
dbWriteTable(conn=db,name="WB_info",df_wb,overwrite=T,append=F,row.names=FALSE)
dbWriteTable(conn=db,name="WB_Mun",df_mun,overwrite=T,append=F,row.names=FALSE)
dbWriteTable(conn=db,name="WB_Lan",df_lan,overwrite=T,append=F,row.names=FALSE)
dbDisconnect(db)



# 
# 
# strSQL<-paste0("SELECT * FROM WB_Lan")
# dfwb_lan <- dbGetQuery(db, strSQL)
# strSQL<-paste0("SELECT * FROM WB")
# wb <- dbGetQuery(db, strSQL)
# dbDisconnect(db)
