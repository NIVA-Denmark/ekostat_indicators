library(tidyverse)
library(haven)
library(lme4)
library(lubridate)
library(prodlim)
library(matrixStats)
library(stats)

library(extrafont)
loadfonts(quiet = TRUE)
theme_set(theme_bw(base_size = 12, base_family = "Open Sans"))

# source("R/extract_modis.R") # ATTENTION: takes few minutes

# Clean the workspace
rm(list = ls())

source("ReadMonitoringData.R")
source("ReadIndicatorParms.R")
#source("CalculateIndicator.R")

source("IndicatorSelectionSweden.R")

# TESTING COASTAL INDICATORS
# Testing the BQI indicator
variance_list <- list(V_station=2.5,V_obspoint=0.64,
                      V_year=0.10,V_yearmonth=0,
                      V_stationyear=0.63,V_stationmonth=0,V_stationdate=0.5,
                      V_institution=0.5,V_replication=0)
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
ParameterVector <- c(0)
MinObsList <- list(MinYear=3,MinObsPerYear=1)
df <- ReadMonitoringDataSMHI("data/Coast_bqi.sas7bdat") %>% filter(EU_CD=="SE580325-113500")
df <- ReadMonitoringDataSMHI("data/Coast_bqi.sas7bdat") %>% filter(EU_CD=="SE652920-222650")
CalculateIndicator("CoastBQI",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)

# MSMDI indicator
variance_list <- list(V_station=0.5211,V_obspoint=0,
                      V_year=0.1911,V_yearmonth=0,
                      V_stationdate=0.3301,
                      V_stationyear=0.3927,V_stationmonth=0,
                      V_institution=0,V_replication=0)
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
ParameterVector <- c(0)
MinObsList <- list(MinYear=3,MinObsPerYear=1)
df <- ReadMonitoringDataSMHI("data/Coast_msmdi.sas7bdat") %>% filter(EU_CD=="SE622011-146303")
CalculateIndicator("CoastMSMDI",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)

# Testing the biovolume indicator
variance_list <- list(V_station=0.23,V_obspoint=0,
                      V_year=0.10,V_yearmonth=0.05,
                      V_stationyear=0.05,V_stationmonth=0.08,V_stationdate=0.05,
                      V_institution=0.05,V_replication=0)
df <- ReadMonitoringDataSMHI("data/Coast_biovol.sas7bdat") %>% filter(EU_CD=="SE552170-130626")
MonthInclude <- c(6,7,8)
MinObsList <- list(MinYear=3,MinObsPerYear=1)
CalculateIndicator("CoastBiovol",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
ParameterVector <- c(35.0,-2.856,7,0.18,0,1)
MinObsList <- list(MinYear=3,MinObsPerYear=1)
CalculateIndicator("CoastBiovolEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)

# Testing the Chlorophyll a indicator
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(EU_CD=="SE552170-130626")
MonthInclude <- c(6,7,8)
MinObsList <- list(MinYear=3,MinObsPerYear=1)
CalculateIndicator("CoastChla",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# type without salinity correction, e.g. type=6
MonthInclude <- c(6,7,8)
ParameterVector <- c(17.0,-0.250,7,0.94,0,1)
MinObsList <- list(MinYear=3,MinObsPerYear=1)
CalculateIndicator("CoastChlaEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# type with salinity correction, e.g. type=8
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(EU_CD=="SE560900-145280")
ParameterVector <- c(33.8,-2.679,7,0.26,0.0051,1.9974)
CalculateIndicator("CoastChlaEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)

df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(EU_CD=="SE560900-145280")
# Testing the indicator for Secchi depth
MonthInclude <- c(6,7,8)
CalculateIndicator("CoastSecchi",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# type without salinity correction, e.g. type=7
ParameterVector <- c(35.0,-2.856,7,10,0,1)
CalculateIndicator("CoastSecchiEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# type with salinity correction, e.g. type=8
ParameterVector <- c(33.8,-2.679,7,0,1023.3,-1.696)
CalculateIndicator("CoastSecchiEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)

# Testing the indicator for nutrients
MonthInclude <- c(6,7,8)
CalculateIndicator("CoastTNsummer",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("CoastTPsummer",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
MonthInclude <- c(11,12,1,2)
CalculateIndicator("CoastTNwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("CoastTPwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("CoastDINwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("CoastDIPwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(23.1,-0.150,7)
CalculateIndicator("CoastTNwinterEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(1.03,-0.076,7)
CalculateIndicator("CoastTPwinterEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(46.9,-6.371,7)
CalculateIndicator("CoastDINwinterEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(0.32,-0.007,7)
CalculateIndicator("CoastDIPwinterEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
MonthInclude <- c(6,7,8)
ParameterVector <- c(35.0,-2.856,7)
CalculateIndicator("CoastTNsummerEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(0.62,-0.046,7)
CalculateIndicator("CoastTPsummerEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)

# Testing coastal oxygen indices
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
ParameterVector <- c(0)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(station=="INSTÖ RÄNNA")
CalculateIndicator("CoastBottomOxygen",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(station=="U9 / LÅNGALMAFJ")
CalculateIndicator("CoastBottomOxygen",df,ParameterVector,MinObsList,variance_list,MonthInclude,2010,2015)
MonthInclude <- c(6,7,8,9,10,11,12)
CalculateIndicator("CoastHypoxicArea",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(station=="SKÄLDERVIKEN")
CalculateIndicator("CoastHypoxicArea",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(station=="BYFJORDEN")
CalculateIndicator("CoastHypoxicArea",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)


# Testing indicator calculation for non-monitored water bodies
variance_list <- list(V_station=0.0472620391974166,V_obspoint=0,
                      V_year=0.000240047417516435,V_yearmonth=0.0192565837622898,
                      V_stationdate=0.0403816847055466,
                      V_stationyear=0,V_stationmonth=0.0205176400753212,
                      V_institution=0.0160297370140365,V_replication=0)
MonthInclude <- c(11,12,1,2)
ParameterVector <- c(0)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(EU_CD=="SE560900-145280")
un1 <-CalculateIndicator("CoastTNwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(EU_CD=="SE574050-114780")
un2 <-CalculateIndicator("CoastTNwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(EU_CD=="SE580688-114860")
un3 <-CalculateIndicator("CoastTNwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(EU_CD=="SE581740-114820")
un4 <-CalculateIndicator("CoastTNwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
df <- ReadMonitoringDataSMHI("data/Coast_watersamples.sas7bdat") %>% filter(EU_CD=="SE581853-112736")
un5 <-CalculateIndicator("CoastTNwinter",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,20128)
un <- list(un1,un2,un3,un4,un5)

CalculateIndicatorType("CoastTNwinter",un,list(V_WBperiod=9,V_WBannual=12),30,2007,2012)

# TESTING LAKE INDICATORS
# Testing of freshwater indicators in lakes
variance_list <- list(V_station=0.5211,V_obspoint=0,
                      V_year=0.1911,V_yearmonth=0,
                      V_stationdate=0.3301,
                      V_stationyear=0.3927,V_stationmonth=0,
                      V_institution=0,V_replication=0)
# Phytoplankton
MonthInclude <- c(7,8)
df <- filter(ReadMonitoringDataSMHI("data/lake_biolindex.sas7bdat"),EU_CD == 'SE615375-137087')
# For LakeBiovol, LakeBiovolEQR, LakeChla, LakeChlaEQR use Gony boundaries, if biovol Gony >5% of biovol total
Gonytest <- mean(df$biovolGony,na.rm=TRUE)/mean(df$biovol,na.rm=TRUE)
CalculateIndicator("LakeBiovol",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(0.46,16)
CalculateIndicator("LakeBiovolEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("LakeBiovolEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
CalculateIndicator("LakePTI",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
ParameterVector <- c(-1,-0.5,rep(0.5,34))
CalculateIndicator("LakePTIEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
CalculateIndicator("LakeNphytspec",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
ParameterVector <- c(rep(45,36))
CalculateIndicator("LakeNphytspecEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
df <- filter(ReadMonitoringDataSMHI("data/lake_wqdata.sas7bdat"),EU_CD == 'SE615375-137087')
CalculateIndicator("LakeChla",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(1.7,61)
CalculateIndicator("LakeChlaEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)

# Macrophytes
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
df <- filter(ReadMonitoringDataSMHI("data/lake_biolindex.sas7bdat"),EU_CD == 'SE615365-134524')
df <- filter(ReadMonitoringDataSMHI("data/lake_biolindex.sas7bdat"),EU_CD == 'SE652177-159038')
ParameterVector <- c(8.27,1,rep(0,34))
CalculateIndicator("LakeTMIEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# Benthic diatoms, should only be used if TP>6 mg/l (months 7 and 8)
df <- filter(ReadMonitoringDataSMHI("data/lake_biolindex.sas7bdat"),EU_CD == 'SE652852-155412')
CalculateIndicator("LakeIPS",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
ParameterVector <- c(19.6,rep(0,35))
CalculateIndicator("LakeIPSEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
ParameterVector <- c(9.6,rep(0,35))
CalculateIndicator("LakeACID",df,ParameterVector,MinObsList,variance_list,MonthInclude,2013,2018)
# Benthic invertebrates
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
df <- filter(ReadMonitoringDataSMHI("data/lake_biolindex.sas7bdat"),EU_CD == 'SE638665-129243')
df <- filter(ReadMonitoringDataSMHI("data/lake_biolindex.sas7bdat"),EU_CD == 'SE652177-159038')
ParameterVector <- c(rep(3.25,36))
CalculateIndicator("LakeBQIEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(rep(5.6,36))
CalculateIndicator("LakeASPTEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(rep(77.5,36))
CalculateIndicator("LakeMILAEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# Fish
MonthInclude <- c(7,8)
df <- filter(ReadMonitoringDataSMHI("data/lake_biolindex.sas7bdat"),EU_CD == 'SE615375-137087')
df <- filter(ReadMonitoringDataSMHI("data/lake_biolindex.sas7bdat"),EU_CD == 'SE664197-149337')
CalculateIndicator("LakeEQR8",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("LakeAindexW5",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("LakeEindexW3",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)

# Nutrients - Total phosphorus
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
df <- filter(ReadMonitoringDataSMHI("data/lake_WQdata.sas7bdat"),EU_CD == 'SE655455-136581')
df <- filter(ReadMonitoringDataSMHI("data/lake_WQdata.sas7bdat"),EU_CD == 'SE670275-146052')
# RefCond for TP: Function depends on availability of data (autumn circulation, annual or August)
AugustOnly <- length(unique(df$month)) == 1 && unique(df$month)[1] == 8
ParameterVector <- c(RefCond_LakeTP(mean(df$Abs_F420,na.rm=TRUE),100,0.6,AugustOnly))
CalculateIndicator("LakeTPEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# Secchi depth
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
df <- filter(ReadMonitoringDataSMHI("data/lake_WQdata.sas7bdat"),EU_CD == 'SE655455-136581')
df <- filter(ReadMonitoringDataSMHI("data/lake_WQdata.sas7bdat"),EU_CD == 'SE670275-146052')
# RefCond for Secchi depth: LakeChla should be second parameter in ParameterVector
ParameterVector <- c(RefCond_LakeSecchiDepth(mean(df$Abs_F420,na.rm=TRUE),8),8)
CalculateIndicator("LakeSecchiDepthEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# Oxygen
MonthInclude <- c(7,8)
CalculateIndicator("LakeOxygenSummer",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)

# pH change

# Testing of freshwater indicators in rivers
variance_list <- list(V_station=0.5211,V_obspoint=0,
                      V_year=0.1911,V_yearmonth=0,
                      V_stationdate=0.3301,
                      V_stationyear=0.3927,V_stationmonth=0,
                      V_institution=0,V_replication=0)
# Benthic diatoms
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
df <- filter(ReadMonitoringDataSMHI("data/river_biolindex.sas7bdat"),EU_CD == 'SE625195-141220')
df <- filter(ReadMonitoringDataSMHI("data/river_biolindex.sas7bdat"),EU_CD == 'SE654141-124734')
CalculateIndicator("RiverIPS",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(19.6)
CalculateIndicator("RiverIPSEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("RiverPctPT",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("RiverTDI",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("RiverACID",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# Benthic invertebrates
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
df <- filter(ReadMonitoringDataSMHI("data/river_biolindex.sas7bdat"),EU_CD == 'SE632137-147160')
df <- filter(ReadMonitoringDataSMHI("data/river_biolindex.sas7bdat"),EU_CD == 'SE638475-137575')
ParameterVector <- c(6.53)
CalculateIndicator("RiverASPTEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(14)
CalculateIndicator("RiverDJEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
ParameterVector <- c(47.5)
CalculateIndicator("RiverMISAEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# Fish
MonthInclude <- c(8,9,10)
df <- filter(ReadMonitoringDataSMHI("data/river_biolindex.sas7bdat"),EU_CD == 'SE624900-134390')
df <- filter(ReadMonitoringDataSMHI("data/river_biolindex.sas7bdat"),EU_CD == 'SE694939-144561')
CalculateIndicator("RiverVIX",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("RiverVIXh",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
CalculateIndicator("RiverVIXsm",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)
# Nutrients - Total phosphorus
MonthInclude <- c(1,2,3,4,5,6,7,8,9,10,11,12)
df <- filter(ReadMonitoringDataSMHI("data/river_WQdata.sas7bdat"),EU_CD == 'SE648930-125964')
df <- filter(ReadMonitoringDataSMHI("data/river_WQdata.sas7bdat"),EU_CD == 'SE632601-145366')
ParameterVector <- c(RefCond_RiverTP(mean(df$Abs_F420,na.rm=TRUE),100,1),rep(0,35))
CalculateIndicator("RiverTPEQR",df,ParameterVector,MinObsList,variance_list,MonthInclude,2007,2012)





