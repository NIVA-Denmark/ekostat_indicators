# AGGREGATION ROUTINES BASED ON MEASUREMENTS
# Aggregation principle used for e.g. coastal chlorophyll (over time points) and BQI (over obspoint)
# Aggregate over stations to yearly means then over years
Aggregate_year_station <- function(df) {
  yearmeans <- df %>%    group_by(year,station) %>%
                         summarise(xvar = mean(xvar,na.rm = TRUE)) %>%
                         group_by(year) %>%
                         summarise(xvar = mean(xvar,na.rm = TRUE))
                         
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Aggregation principle used for e.g. nutrients
# Aggregate over years and then period
Aggregate_year <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvar,na.rm = TRUE))
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}
Min_year <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = min(xvar,na.rm = TRUE))
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}
Max_year <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = max(xvar,na.rm = TRUE))
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}


# Aggregation principle used for e.g. Secchi depth
# Aggregate over entire period
Aggregate_period <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvar,na.rm = TRUE))
  
  periodmean <- mean(df$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# AGGREGATION PRINCIPLES BASED ON EQR VALUES
# Aggregation principle used for e.g. biovolume in lakes
# Aggregate over entire period and then calculate EQR value
Aggregate_period_P_EQR <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(RefCond,na.rm = TRUE)/mean(xvar,na.rm = TRUE))
  
  periodmean <- mean(df$RefCond)/mean(df$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Aggregation principle used for e.g. biovolume in lakes
# Aggregate over entire period and then calculate EQR value
Aggregate_period_N_EQR <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvar,na.rm = TRUE)/mean(RefCond,na.rm = TRUE))
  
  periodmean <- mean(df$xvar)/mean(df$RefCond)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Aggregation principle used for lakes in remiss
# Aggregate over entire period and then calculate EQR value based on transformation (x-max)/(ref-max)
Aggregate_period_RefMax_EQR <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = (mean(xvar,na.rm = TRUE)-mean(MaxCond,na.rm = TRUE))/(mean(RefCond,na.rm = TRUE)-mean(MaxCond,na.rm = TRUE)))
  
  periodmean <- (mean(df$xvar)-mean(df$MaxCond))/(mean(df$RefCond)-mean(df$MaxCond))
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Aggregation principle used for lakes in remiss
# Aggregate over entire period and then calculate EQR value based on transformation (x-max)/(ref-max), truncation of values <0 and >1
Aggregate_period_RefMax_EQRtrunc <- function(df) {
  df <- mutate(df,xvarEQR = ifelse(RefCond<MaxCond,ifelse(xvar<RefCond,1,ifelse(xvar>MaxCond,0,(xvar-MaxCond)/(RefCond-MaxCond))),ifelse(xvar>RefCond,1,ifelse(xvar<MaxCond,0,(xvar-MaxCond)/(RefCond-MaxCond)))))
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = (mean(xvarEQR,na.rm = TRUE)))
  
  periodmean <- mean(df$xvarEQR)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Aggregation principle used for DJ-index in rivers
# Aggregate over entire period and then calculate EQR value after subtracting the value 5
Aggregate_period_N_EQR5 <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = (mean(xvar,na.rm = TRUE)-5)/(mean(RefCond,na.rm = TRUE)-5))
  
  periodmean <- mean(df$xvar)/mean(df$RefCond)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Aggregation principle used for proportions with EQR calculation
# Aggregate over entire period and then calculate EQR value
Aggregate_period_Prop_EQR <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = (1-mean(xvar,na.rm = TRUE))/(1-mean(RefCond,na.rm = TRUE)))

  periodmean <- (1-mean(df$xvar))/(1-mean(df$RefCond))
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Aggregation principle used for proportions with EQR calculation
# Aggregate over entire period and then calculate EQR value
Aggregate_period_TPI_EQR <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(x = mean(xvar,na.rm = TRUE),r50 = mean(RefCond),r75 = mean(HG_boundary),xvar = (r75-r50)/(x+r75-2*r50))
  yearmeans <- yearmeans %>% mutate(x = NULL,r50 = NULL, r75 = NULL)
  
  r50 = mean(df$RefCond)
  r75 = mean(df$HG_boundary)
  periodmean <- (r75-r50)/(mean(df$xvar)+r75-2*r50)
  
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}


# AGGREGATION PRINCIPLES BASED ON EQR OBSERVATIONS
# Aggregation principle used for e.g. coastal chlorophyll as EQR
# Compute EQR values and then aggregate over stations, years and period
AggregateEQRtrunc_year_station <- function(df) {
  
  df <- mutate(df,xvarEQR = ifelse(xvar<RefCond,1,RefCond/xvar))

  yearmeans <- df %>%    group_by(year,station) %>%
    summarise(xvarEQR = mean(xvarEQR,na.rm = TRUE)) %>%
    group_by(year) %>%
    summarise(xvar = mean(xvarEQR,na.rm = TRUE))   # should be returned in xvar
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)  
}

# Aggregation principle used for e.g. nutrients as EQR
# Aggregate over years and then period
AggregateEQR_year <- function(df) {
  
  df <- mutate(df,xvarEQR = RefCond/xvar)
  
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvarEQR,na.rm = TRUE))   # should be returned in xvar
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)  
}
MaxEQR_year <- function(df) {
  
  df <- mutate(df,xvarEQR = RefCond/xvar)
  
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = max(xvarEQR,na.rm = TRUE))   # should be returned in xvar
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)  
}

# Aggregate over entire period
# Indicator response positive to degradation, i.e. chlorophyll
AggregateEQR_P_period <- function(df) {
  
  df <- mutate(df,xvarEQR = RefCond/xvar)
  
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvarEQR,na.rm = TRUE))   # should be returned in xvar
  
  periodmean <- mean(df$xvarEQR)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)  
}

# Aggregate over entire period
# Indicator response negative to degradation, i.e. Secchi depth
AggregateEQR_N_period <- function(df) {
  
  df <- mutate(df,xvarEQR = xvar/RefCond)
  
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvarEQR,na.rm = TRUE))   # should be returned in xvar
  
  periodmean <- mean(df$xvarEQR)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)  
}

# Aggregate over entire period with truncation of values above 1
# Indicator response negative to degradation, e.g. LakeFishAindexW5 and LakeFishEindexW3
AggregateEQRtrunc_N_period <- function(df) {
  
  df <- mutate(df,xvarEQR = ifelse(xvar>RefCond,1,xvar/RefCond))
  
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = mean(xvarEQR,na.rm = TRUE))   # should be returned in xvar
  
  periodmean <- mean(df$xvarEQR)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)  
}

# SPECIAL CASES FOR COASTAL INDICATORS
# Calculation of BQI indicator according to Handbook
BQIbootstrap <- function(df) {
  error_code <- 0
  yearstatmeans <- df %>%    group_by(year,station) %>%
    summarise(xvar = mean(xvar,na.rm = TRUE))
  Nyearstat <- yearstatmeans %>% group_by(year) %>% summarise(n_station = length(xvar))
  # Check if fewer than 5 stations and years - error_code=-93
  if (sum(Nyearstat$n_station)<5) error_code <- -93
  BQIsimyear <- mat.or.vec(length(Nyearstat$n_station),9999)
  BQIsimperiod <- mat.or.vec(9999,1)
  for (isim in 1:9999) {
    for(i in 1:length(Nyearstat$n_station)) {
      BQIsim <- trunc(runif(Nyearstat$n_station[i],1,Nyearstat$n_station[i]+1))
      BQIsim <- yearstatmeans$xvar[yearstatmeans$year == Nyearstat$year[i]][BQIsim]
      BQIsimyear[i,isim] <- mean(BQIsim)
    }
    BQIsimperiod[isim] <- mean(BQIsimyear[,isim])
  }
  periodmean <- quantile(BQIsimperiod,probs=0.2)
  yearmeans <- data.frame(year=Nyearstat$year,xvar = apply(BQIsimyear,1,quantile,probs=0.2,na.rm=TRUE))
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=error_code)
  return(res)  
}

# Oxygen in bottom water - calculate lower quartile
OxygenLowerQuartile <- function(df) {
  yearmeans <- df %>% group_by(year) %>%
    summarise(xvar = quantile(xvar,probs=c(0.25),na.rm = TRUE))
  
  periodmean <- mean(yearmeans$xvar)
  res <- list(periodmean=periodmean,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Calculation of Oxygen indicator according to Handbook
# First test calculates the average of O2 observations below the 25%-percentile - threshold is 3.5 ml/l
# Second test calculates the average of O2 observations (Jan-May) below the 25%-percentile 
OxygenTest1 <- function(df) {
  years <- df %>% group_by(year) %>% summarise()
  df1 <- filter(df,xvar<=quantile(xvar,na.rm=TRUE)[2])
  O2_test1 <- mean(df1$xvar)
  O2_test1_yearmeans <- df1 %>% group_by(year) %>% summarise(xvar = mean(xvar))
  O2_test1_yearmeans <- left_join(years,O2_test1_yearmeans,c("year"))
  df2 <- df1 %>% filter(month %in% c(1,2,3,4,5)) 
  O2_test2 <- mean(df2$xvar)
  O2_test2_yearmeans <- df2 %>% group_by(year) %>% summarise(xvar = mean(xvar))
  O2_test2_yearmeans <- left_join(years,O2_test2_yearmeans,c("year"))
  yearmeans <- data.frame(year=O2_test1_yearmeans$year,xvar = O2_test1_yearmeans$xvar)
  res <- list(periodmean=O2_test1,yearmeans=yearmeans,error_code=0)
  return(res)
}

# Calculation of Oxygen indicator according to Handbook
# First test calculates the average of O2 observations below the 25%-percentile - threshold is 3.5 ml/l
# Second test calculates the average of O2 observations (Jan-May) below the 25%-percentile 
# df contains oxygen profiles 
OxygenTest2 <- function(df) {
  df <- filter(df,!is.na(xvar))
  # create list of years for producing vectors of similar length in the returned list
  years <- df %>% group_by(year) %>% summarise()
  # Ensure that the profiles are sorted by depth and contain oxygen concentrations
  df <- df %>% group_by(station,date,time,station_depth,depth) %>% filter(is.na(xvar) == FALSE)
  # Use measured profile from surface until O2 concentrations has decreased >1 ml/l
  O2surface <- df %>% group_by(station,date,time,station_depth) %>% summarise(O2surface = O2[which.min(depth)])
  df <- full_join(df,O2surface,c("station","date","time","station_depth"))
  df$xvar <- ifelse(df$O2<df$O2surface-1,df$xvar,df$O2)
  # find depth for 3.5 ml/l by extrapolation, when this threshold is not in the profile
  O2bottom_ext <- df %>% group_by(station,date,time,station_depth) %>% summarise(n_obs =n(),
                                                                                 O2bottom1 = xvar[which.max(depth)],O2bottom2 = ifelse(n_obs>1,xvar[which.max(depth)-1],NA),
                                                                                 depth1 = depth[which.max(depth)],depth2 =ifelse(n_obs>1,depth[which.max(depth)-1],NA),
                                                                                 O2clinedepth_max=ifelse(n_obs>1,ifelse(O2bottom1>3.5 && O2bottom2-O2bottom1>0,depth1+(3.5-O2bottom1)/(O2bottom2-O2bottom1)*(depth2-depth1),1000),NA))  # find profile statistics for tests and indicator calculation
  O2bottom <- df %>% group_by(station,date,time,station_depth) %>% summarise(n_obs =n(), O2range = range(xvar)[2]-range(xvar)[1],
                                                                             max_depth=max(depth),O2bottom = xvar[which.max(depth)],O2clinedepth = ifelse(n_obs>1 && O2range>0,approx(xvar,depth,c(3.5))$y,NA))
  # If O2clinedepth is not in the profile then use the extrapolated value
  O2bottom <- full_join(O2bottom,O2bottom_ext,c("station","date","time","station_depth"))
  O2bottom <- O2bottom %>% mutate(O2clinedepth = ifelse(is.na(O2clinedepth),O2clinedepth_max,O2clinedepth),depth1 = NULL, depth2 = NULL, O2bottom1 = NULL, O2bottom2 = NULL, O2clinedepth_max = NULL)
  # Find the percent area affected by O2 concentrations <3.5 ml/l
  O2bottom <- O2bottom %>% mutate(area_hyp = 100-approx(WB_bathymetry$depth,WB_bathymetry$area_pct,O2clinedepth,yleft=0,yright=100)$y)
  # Calculate test1 as average of O2 observations at bottom (<1.5 m from bottom depth) below the 25-percentile for Jan-Dec
  lower_quantile <- quantile(O2bottom$O2bottom,na.rm=TRUE)[2]
  df1 <- O2bottom %>% filter(station_depth-max_depth<1.5) %>% filter(O2bottom<=lower_quantile) %>% mutate(year=lubridate::year(date))
  #  df1 <- df %>% filter(station_depth-depth<1.5) %>% filter(xvar<=quantile(xvar,na.rm=TRUE)[2])
  O2_test1 <- mean(df1$O2bottom)
  O2_test1_yearmeans <- df1 %>% group_by(year) %>% summarise(O2bottom = mean(O2bottom))
  # Complete with all years having O2 data
  O2_test1_yearmeans <- left_join(years,O2_test1_yearmeans,c("year"))
  # Calculate EQR from Table 7.1 in Handbook
  EQR_test1 <- approx(c(-5.0,0.0,1.0,2.1,3.5,7.0),c(0,0.2,0.4,0.6,0.8,1),O2_test1,yleft=0,yright=1)$y
  EQR_test1_yearmeans <- approx(c(-10.0,0.0,1.0,2.1,3.5,7.0),c(0,0.2,0.4,0.6,0.8,1),O2_test1_yearmeans$O2bottom,yleft=0,yright=1)$y
  # Calculate test2 as average of O2 concentrations below the 25-percentile for Jan-May
  df2 <- O2bottom %>% mutate(month = lubridate::month(date),year = lubridate::year(date)) %>% filter(month %in% c(1,2,3,4,5))
  # Return from function if no observations available (Jan-May) for calculation of O2_test2
  if (nrow(df2) == 0) {
    yearmeans <- df %>% group_by(year) %>% summarise(xvar = NA)  # Return NA values in res
    res <- list(periodmean=NA,yearmeans=yearmeans,error_code=-91)
    return(list(error_code=-91))
  }
  O2_test2 <- mean(df2$O2bottom)
  O2_test2_yearmeans <- df2 %>% group_by(year) %>% summarise(O2bottom = mean(O2bottom))
  # Complete with all years having O2 data
  O2_test2_yearmeans <- left_join(years,O2_test2_yearmeans,c("year"))
  # Calculate indicator for percent area affected by <3.5 ml/l
  df2 <- O2bottom %>% mutate(month = lubridate::month(date),year = lubridate::year(date)) %>% filter(month %in% c(6,7,8,9,10,11,12))
  # Return from function if no observations available (Jun-Dec) for calculation of hypoxic area
  if (nrow(df2) == 0) {
    yearmeans <- df %>% group_by(year) %>% summarise(xvar = NA)  # Return NA values in res
    res <- list(periodmean=NA,yearmeans=yearmeans,error_code=-92)
    return(list(error_code=-92))
  }
  hyparea <- mean(df2$area_hyp)
  hyparea_yearmeans <- df2 %>% group_by(year) %>% summarise(area_hyp = mean(area_hyp))
  # Complete with all years having O2 data
  hyparea_yearmeans <- left_join(years,hyparea_yearmeans,c("year"))
  # Calculate EQR from Table 7.1 in Handbook
  EQR_test2 <- approx(BoundariesHypoxicArea,c(0,0.2,0.4,0.6,0.8,1),hyparea,yleft=0,yright=1)$y
  EQR_test2_yearmeans <- approx(BoundariesHypoxicArea,c(0,0.2,0.4,0.6,0.8,1),hyparea_yearmeans$area_hyp,yleft=0,yright=1)$y
  O2_test1<-ifelse(is.nan(O2_test1),0,O2_test1) # Set O2_test1 to zero if no data to complete the if-clause below
  if (O2_test1>3.5 || O2_test2>3.5) {
    res <- list(periodmean=EQR_test1,yearmeans=data.frame(year=O2_test1_yearmeans$year,xvar = EQR_test1_yearmeans),error_code=0)
  } else {
    res <- list(periodmean=EQR_test2,yearmeans=data.frame(year=O2_test1_yearmeans$year,xvar = EQR_test2_yearmeans),error_code=0)
  }
  return(res)
}






