# Include supporting routines for indicator calculations
source("CalculateIndicatorSupport.R")
# Include routines for calculating indicators
source("IndicatorFunctions.R")


#' Generic routine for calculating indicator statistics
#' 
#' @param Indicator The name identifier for the indicator to be calculated. 
#' @param df A dataframe with monitoring data from the Swedish Monitoring program. 
#'   The dataframe should contain the
#'   following variables:
#'   
#'   \describe{ 
#'   \item{station}{An identifier for the monitoring station.} 
#'   \item{date}{Date of the observation.} 
#'   \item{institution}{Provider of the observation.} 
#'   \item{chla}{Chlorophyll a concentration in sample.} 
#'   
#' @param ParameterVector A vector of parameters to be parsed to the indicator calcutions (Ref Cond, salinity correction, etc.)
#' @param MinObsList A parameter list with official data requirements that will flag the indicator calculation 
#' @param var_list List of variance components
#' @param MonthInclude A list of month numbers for filtering data
#' @param startyear The first year in the indicator calculation
#' @param endyear The last year in the indicator calculation
#' @param n_iter Number of iterations for Monte Carlo simulation
#'   
#' @return
#' @export
#' 
#' @examples 

CalculateIndicator <-
  function(Indicator,df,ParameterVector,MinObsList,var_list,MonthInclude,startyear,endyear,n_iter=1000) {
    # Set flag to zero and change it for error handling below
    flag <- 0
    # Select the observation variable for the indicator
    xvar <- switch(Indicator,
                   CoastBQI          = df$BQI,
                   CoastMSMDI        = df$MSMDI,
                   CoastChla         = df$chla,
                   CoastChlaEQR      = df$chla,
                   CoastBiovol       = df$biovol,
                   CoastBiovolEQR    = df$biovol,
                   CoastSecchi       = df$secchi,
                   CoastSecchiEQR    = df$secchi,
                   CoastDINwinter    = df$DIN,
                   CoastDINwinterEQR = df$DIN,
                   CoastDIPwinter    = df$DIP,
                   CoastDIPwinterEQR = df$DIP,
                   CoastTNsummer     = df$TN,
                   CoastTNsummerEQR  = df$TN,
                   CoastTNwinter     = df$TN,
                   CoastTNwinterEQR  = df$TN,
                   CoastTPsummer     = df$TP,
                   CoastTPsummerEQR  = df$TP,
                   CoastTPwinter     = df$TP,
                   CoastTPwinterEQR  = df$TP,
                   CoastOxygen       = df$O2,
                   CoastBottomOxygen = df$O2_bot,
                   CoastHypoxicArea  = df$HypoxicAreaPct,
                   LakeBiovol        = df$biovol,
                   LakeBiovolEQR     = df$biovol,
                   LakeChla          = df$chla,
                   LakeChlaEQR       = df$chla,
                   LakePropCyano     = df$Proportion_cyanobacteria,
                   LakePropCyanoEQR  = df$Proportion_cyanobacteria,
                   LakePTI           = df$PhytoplanktonTrophicIndex,
                   LakePTIEQR        = df$PhytoplanktonTrophicIndex,
                   LakeNphytspec     = df$Nspecies_phytoplankton,
                   LakeNphytspecEQR  = df$Nspecies_phytoplankton,
                   LakeTMIEQR        = df$TrophicMacrophyteIndex,
                   LakeIPS           = df$BenthicDiatomsIPS,
                   LakeIPSEQR        = df$BenthicDiatomsIPS,
                   LakeACID          = df$BenthicDiatomsACID,
                   LakeACIDEQR       = df$BenthicDiatomsACID,
                   LakeASPTEQR       = df$BenthicInvertebratesASPT,
                   LakeBQIEQR        = df$BenthicInvertebratesBQI,
                   LakeMILAEQR       = df$BenthicInvertebratesMILA,
                   LakeEQR8          = df$EQR8,
                   LakeAindexW5      = df$AindexW5,
                   LakeEindexW3      = df$EindexW3,
                   LakeTPEQR         = df$TP,
                   LakeSecchiDepthEQR= df$SecchiDepth,
                   LakeOxygenSummer  = df$O2_bot,
                   RiverIPS          = df$BenthicDiatomsIPS,
                   RiverIPSEQR       = df$BenthicDiatomsIPS,
                   RiverPctPT        = df$BenthicDiatomsPctPT,
                   RiverTDI          = df$BenthicDiatomsTDI,
                   RiverACID         = df$BenthicDiatomsACID,
                   RiverASPTEQR      = df$BenthicInvertebratesASPT,
                   RiverDJEQR        = df$BenthicInvertebratesDJ,
                   RiverMISAEQR      = df$BenthicInvertebratesMISA,
                   RiverVIX          = df$VIX,
                   RiverVIXh         = df$VIXh,
                   RiverVIXsm        = df$VIXsm,
                   RiverTPEQR        = df$TP,
                   RiverOxygenSummer = df$O2_bot
    )
    df <- mutate(df,xvar=xvar)
    # Associating indicators with transformation from observations
    f_fun <- switch(Indicator,
                    CoastBQI          = Aggregate_year_station,
                    CoastMSMDI        = Aggregate_period,
                    CoastChla         = Aggregate_year_station,
                    CoastChlaEQR      = AggregateEQRtrunc_year_station,
                    CoastBiovol       = Aggregate_year_station,
                    CoastBiovolEQR    = AggregateEQRtrunc_year_station,
                    CoastSecchi       = Aggregate_period,
                    CoastSecchiEQR    = AggregateEQR_N_period,
                    CoastDINwinter    = Aggregate_year, # Replace with Max_year?
                    CoastDINwinterEQR = AggregateEQR_year, # Replace with MaxEQR_year?
                    CoastDIPwinter    = Aggregate_year, # Replace with Max_year?
                    CoastDIPwinterEQR = AggregateEQR_year, # Replace with MaxEQR_year?
                    CoastTNsummer     = Aggregate_year, 
                    CoastTNsummerEQR  = AggregateEQR_year,
                    CoastTNwinter     = Aggregate_year, # Replace with Max_year?
                    CoastTNwinterEQR  = AggregateEQR_year, # Replace with MaxEQR_year?
                    CoastTPsummer     = Aggregate_year,
                    CoastTPsummerEQR  = AggregateEQR_year,
                    CoastTPwinter     = Aggregate_year, # Replace with Max_year?
                    CoastTPwinterEQR  = AggregateEQR_year, # Replace with MaxEQR_year?
                    CoastOxygen       = OxygenTest2,
                    CoastBottomOxygen = OxygenLowerQuartile,
                    CoastHypoxicArea  = Aggregate_year,
                    LakeBiovol        = Aggregate_period,
                    LakeBiovolEQR     = Aggregate_period_RefMax_EQRtrunc,
                    LakeChla          = Aggregate_period,
                    LakeChlaEQR       = Aggregate_period_RefMax_EQRtrunc,
                    LakePropCyano     = Aggregate_period,
                    LakePropCyanoEQR  = Aggregate_period_Prop_EQR,
                    LakePTI           = Aggregate_period,
                    LakePTIEQR        = Aggregate_period_RefMax_EQRtrunc,
                    LakeNphytspec     = Aggregate_period,
                    LakeNphytspecEQR  = Aggregate_period_N_EQR,
                    LakeTMIEQR        = Aggregate_period_RefMax_EQRtrunc,
                    LakeIPS           = Aggregate_period,
                    LakeIPSEQR        = Aggregate_period_N_EQR,
                    LakeACID          = Aggregate_period,
                    LakeACIDEQR       = Aggregate_period_N_EQR,
                    LakeASPTEQR       = Aggregate_period_N_EQR,
                    LakeBQIEQR        = Aggregate_period_N_EQR,
                    LakeMILAEQR       = Aggregate_period_N_EQR,
                    LakeEQR8          = Aggregate_period,
                    LakeAindexW5      = AggregateEQRtrunc_N_period,
                    LakeEindexW3      = AggregateEQRtrunc_N_period,
                    LakeTPEQR         = Aggregate_period_P_EQR,
                    LakeSecchiDepthEQR= Aggregate_period_N_EQR,
                    LakeOxygenSummer  = Min_year,
                    RiverIPS          = Aggregate_period,
                    RiverIPSEQR       = Aggregate_period_N_EQR,
                    RiverPctPT        = Aggregate_period,
                    RiverTDI          = Aggregate_period,
                    RiverACID         = Aggregate_period,
                    RiverASPTEQR      = Aggregate_period_N_EQR,
                    RiverDJEQR        = Aggregate_period_N_EQR5,
                    RiverMISAEQR      = Aggregate_period_N_EQR,
                    RiverVIX          = Aggregate_period,
                    RiverVIXh         = Aggregate_period,
                    RiverVIXsm        = Aggregate_period,
                    RiverTPEQR        = Aggregate_period_P_EQR,
                    RiverOxygenSummer = Min_year
    )
    # Assigning transformations for measurements to obtain normal distributed variates
    g_fun <- switch(Indicator,
                    CoastBQI          = identity,
                    CoastMSMDI        = logit_w_replace,
                    CoastChla         = log,
                    CoastChlaEQR      = log,
                    CoastBiovol       = log,
                    CoastBiovolEQR    = log,
                    CoastSecchi       = log,
                    CoastSecchiEQR    = log,
                    CoastDINwinter    = log,
                    CoastDINwinterEQR = log,
                    CoastDIPwinter    = log,
                    CoastDIPwinterEQR = log,
                    CoastTNsummer     = log,
                    CoastTNsummerEQR  = log,
                    CoastTNwinter     = log,
                    CoastTNwinterEQR  = log,
                    CoastTPsummer     = log,
                    CoastTPsummerEQR  = log,
                    CoastTPwinter     = log,
                    CoastTPwinterEQR  = log,
                    CoastOxygen       = identity,
                    CoastBottomOxygen = identity,
                    CoastHypoxicArea  = logit_w_replace,
                    LakeBiovol        = log,
                    LakeBiovolEQR     = log,
                    LakeChla          = log,
                    LakeChlaEQR       = log,
                    LakePropCyano     = logit_w_replace,
                    LakePropCyanoEQR  = logit_w_replace,
                    LakePTI           = identity,
                    LakePTIEQR        = identity,
                    LakeNphytspec     = log,
                    LakeNphytspecEQR  = log,
                    LakeTMIEQR        = identity,
                    LakeIPS           = identity,
                    LakeIPSEQR        = identity,
                    LakeACID          = identity,
                    LakeACIDEQR       = identity,
                    LakeASPTEQR       = identity,
                    LakeBQIEQR        = identity,
                    LakeMILAEQR       = identity,
                    LakeEQR8          = logit_w_replace,
                    LakeAindexW5      = logit_w_replace,
                    LakeEindexW3      = logit_w_replace,
                    LakeTPEQR         = log,
                    LakeSecchiDepthEQR= identity,
                    LakeOxygenSummer  = identity,
                    RiverIPS          = identity,
                    RiverIPSEQR       = identity,
                    RiverPctPT        = logit_w_replace,
                    RiverTDI          = identity,
                    RiverACID         = identity,
                    RiverASPTEQR      = identity,
                    RiverDJEQR        = identity,
                    RiverMISAEQR      = identity,
                    RiverVIX          = logit_w_replace,
                    RiverVIXh         = logit_w_replace,
                    RiverVIXsm        = logit_w_replace,
                    RiverTPEQR        = log,
                    RiverOxygenSummer = identity
    )    
    # Assigning inverse transformations of g_fun
    g_fun_inv <- switch(Indicator,
                        CoastBQI          = identity,
                        CoastMSMDI        = plogis,
                        CoastChla         = exp,
                        CoastChlaEQR      = exp,
                        CoastBiovol       = exp,
                        CoastBiovolEQR    = exp,
                        CoastSecchi       = exp,
                        CoastSecchiEQR    = exp,
                        CoastDINwinter    = exp,
                        CoastDINwinterEQR = exp,
                        CoastDIPwinter    = exp,
                        CoastDIPwinterEQR = exp,
                        CoastTNsummer     = exp,
                        CoastTNsummerEQR  = exp,
                        CoastTNwinter     = exp,
                        CoastTNwinterEQR  = exp,
                        CoastTPsummer     = exp,
                        CoastTPsummerEQR  = exp,
                        CoastTPwinter     = exp,
                        CoastTPwinterEQR  = exp,
                        CoastOxygen       = identity,
                        CoastBottomOxygen = identity,
                        CoastHypoxicArea  = plogis,
                        LakeBiovol        = exp,
                        LakeBiovolEQR     = exp,
                        LakeChla          = exp,
                        LakeChlaEQR       = exp,
                        LakePropCyano     = plogis,
                        LakePropCyanoEQR  = plogis,
                        LakePTI           = identity,
                        LakePTIEQR        = identity,
                        LakeNphytspec     = exp,
                        LakeNphytspecEQR  = exp,
                        LakeTMIEQR        = identity,
                        LakeIPS           = identity,
                        LakeIPSEQR        = identity,
                        LakeACID          = identity,
                        LakeACIDEQR       = identity,
                        LakeASPTEQR       = identity,
                        LakeBQIEQR        = identity,
                        LakeMILAEQR       = identity,
                        LakeEQR8          = plogis,
                        LakeAindexW5      = plogis,
                        LakeEindexW3      = plogis,
                        LakeTPEQR         = exp,
                        LakeSecchiDepthEQR= identity,
                        LakeOxygenSummer  = identity,
                        RiverIPS          = identity,
                        RiverIPSEQR       = identity,
                        RiverPctPT        = plogis,
                        RiverTDI          = identity,
                        RiverACID         = identity,
                        RiverASPTEQR      = identity,
                        RiverDJEQR        = identity,
                        RiverMISAEQR      = identity,
                        RiverVIX          = plogis,
                        RiverVIXh         = plogis,
                        RiverVIXsm        = plogis,
                        RiverTPEQR        = exp,
                        RiverOxygenSummer = identity
    ) 
    # Add month and year to df
    df <- mutate(df,month=month(date))
    df <- mutate(df,year=year(date))
    # Switch year for winter months (Nov+Dec) to include together with (Jan+Feb)
    if (Indicator %in% c("CoastDINwinterEQR","CoastDIPwinterEQR","CoastTNwinter","CoastTNwinterEQR","CoastTPwinter","CoastTPwinterEQR")) {
      df <- mutate(df,year=ifelse(month %in% c(11,12),year+1,year))
    }
    # Filter dataframe to include observations used in indicator only
    df <- Filter_df(df,MonthInclude,startyear,endyear) 
    # setting RefCond depending on salinity for indicators with salinity correction
    RefCond <- mat.or.vec(nrow(df), 1)
    if (Indicator %in% c("CoastChlaEQR","CoastBiovolEQR","CoastSecchiEQR","CoastDINwinterEQR","CoastDIPwinterEQR","CoastTNsummerEQR","CoastTPsummerEQR","CoastTNwinterEQR","CoastTPwinterEQR")) {
      df <- filter(df,!is.na(sali))
      df$sali <- ifelse(df$sali<2,2,df$sali) # Set RefCond for sali<2 to the value at sali=2
      df$sali <- ifelse(df$sali>ParameterVector[3],ParameterVector[3],df$sali) # Set RefCond for sali higher than outer to the value at sali=outer boundary
      RefCond <- ParameterVector[1]+ParameterVector[2]*df$sali
      if (Indicator %in% c("CoastChlaEQR","CoastBiovolEQR","CoastSecchiEQR"))
         RefCond <- ParameterVector[4]+ParameterVector[5]*RefCond^ParameterVector[6]
      df <- mutate(df,RefCond = RefCond)
    }
    # setting RefCond and MaxCond depending for lake and river indicators, i.e. RefCond in ParameterVector[1] and MaxCond in ParameterVector[2]
    if (Indicator %in% c("LakeBiovolEQR","LakeChlaEQR","LakePTIEQR","LakeNphytspecEQR","LakeTMIEQR","LakeIPSEQR","LakeASPTEQR","LakeBQIEQR","LakeMILAEQR","LakeEQR8","LakeAindexW5","LakeEindexW3","LakeTPEQR","LakeSecchiDepthEQR",
                         "RiverIPSEQR","RiverASPTEQR","RiverDJEQR","RiverMISAEQR","RiverTPEQR")) {
      RefCond <- mat.or.vec(nrow(df), 1)
      MaxCond <- mat.or.vec(nrow(df), 1)
      for (i in 1:nrow(df)) {
        RefCond[i] <- ParameterVector[1]
        MaxCond[i] <- ParameterVector[2]
      }
      df <- mutate(df,RefCond = RefCond,MaxCond = MaxCond) 
    }
    # Adding H-G boundary to data for TPI calculations of EQR
    if (Indicator %in% c("LakeTPIEQR")) {
      HG_boundary <- mat.or.vec(nrow(df), 1)
      for (i in 1:nrow(df)) {HG_boundary[i] <- ParameterVector[2]}
      df <- mutate(df,HG_boundary = HG_boundary) 
    }
    # Calculate number of years, stations, months, institutions and combinations thereof in df 
    ndf <- DF_Ncalculation(df)
    # Return from function if no observations for calculation
    if (ndf$n_obs == 0) return(list(result_code=-90))
    # Check data requirements and flag if not fulfilled
    # Check min number year with at least MinObsPerYear per year
    flag <- ifelse(sum((df %>% group_by(year) %>% summarise(n_year=length(xvar)))$n_year>MinObsList$MinObsPerYear-1)<MinObsList$MinYear,-2,flag)
    # Check min number of years
    flag <- ifelse(ndf$n_year<MinObsList$MinYear,-1,flag)
    # Estimate mean of the transformed observation for simulation
    alpha <- df %>% group_by(year) %>% summarise(mean = mean(g_fun(xvar),na.rm=TRUE))
    # Calculate indicator
    mu_indicator <- f_fun(df)
    # Return from function if no observations in Jan-May (result_code=-91) or Jun-Dec (result_code=-92) or less than 5 obs for BQI (result_code=-93)
    if (mu_indicator$error_code != 0) return(list(result_code=mu_indicator$error_code))
    # Simulate system with random variables for estimating the variance of the indicator
    simres <- vector("numeric",n_iter)
    simresyear <- matrix(nrow=ndf$n_year,ncol=n_iter)
    simrescode <- vector("numeric",n_iter)
    
    # simulation loop - simres contains the residuals from n_iter simulations
    for (isim in 1:n_iter) {
      # simulate variations in the random factors using the data structure
      if (Indicator == "CoastOxygen") {
        simulobs <- SetVector_IndicatorSimO2(alpha$mean,ndf,var_list,df,length(MonthInclude))
      } else {
        simulobs <- SetVector_IndicatorSim(alpha$mean,ndf,var_list,df,length(MonthInclude))
      }
      # backtransform simulations from transformed domain to original domain
      simulobs <- g_fun_inv(simulobs)
      # add simulated observation to df
      simul_df <- df %>% mutate(xvar = NULL, xvar=simulobs)
      # Calculate indicator value for each year and period
      simul_indicator <- f_fun(simul_df)
      # Check for errors in simulations
      simrescode[isim] <- simul_indicator$error_code
      if (simul_indicator$error_code == 0) {
        simresyear[,isim]=simul_indicator$yearmeans$xvar
        simres[isim] <- simul_indicator$periodmean
      } else{
        simresyear[,isim] = NA
        simres[isim] = NA
      }
    } # end simulation loop
    
    # Adjust simulations to have zero mean and then add indicator means - bias correction
    # For some indicators the bias correction needs to be adjusted, because the transformation approximation does not hold (Aggregate_period_RefMax_EQRtrunc)
    if(Indicator %in% c("LakeBiovolEQR","LakeChlaEQR","LakePTIEQR","LakeTMIEQR")) {
      simres <- simres-mean(simres)+mu_indicator$periodmean
      simres <- ifelse(simres<0,0,simres)
      simres <- ifelse(simres>1,1,simres)
      simresyear <- simresyear-apply(simresyear,1,mean)+mu_indicator$yearmeans$xvar
    }
    else {    
      simres <- g_fun_inv(g_fun(simres)-g_fun(mean(simres))+g_fun(mu_indicator$periodmean))
      simresyear <- g_fun_inv(g_fun(simresyear)-g_fun(apply(simresyear,1,mean))+g_fun(mu_indicator$yearmeans$xvar))
    }
    
    # Calculate statistics
    period <- data.frame(mean=mean(simres),stderr=sd(simres),
                         lower_1  = quantile(simres,probs=0.005,na.rm=TRUE),
                         lower_5  = quantile(simres,probs=0.025,na.rm=TRUE),
                         lower_10 = quantile(simres,probs=0.05,na.rm=TRUE),
                         upper_10 = quantile(simres,probs=0.95,na.rm=TRUE),
                         upper_5  = quantile(simres,probs=0.975,na.rm=TRUE),
                         upper_1  = quantile(simres,probs=0.995,na.rm=TRUE),row.names = NULL)
    annual <- data.frame(year = mu_indicator$yearmeans$year,mean = apply(simresyear,1,mean),stderr = apply(simresyear,1,sd),
                         lower_1  = apply(simresyear,1,quantile,probs=0.005,na.rm=TRUE),
                         lower_5  = apply(simresyear,1,quantile,probs=0.025,na.rm=TRUE),
                         lower_10 = apply(simresyear,1,quantile,probs=0.05,na.rm=TRUE),
                         upper_10 = apply(simresyear,1,quantile,probs=0.95,na.rm=TRUE),
                         upper_5  = apply(simresyear,1,quantile,probs=0.975,na.rm=TRUE),
                         upper_1  = apply(simresyear,1,quantile,probs=0.995,na.rm=TRUE))

    res <- list(period=period,annual=annual,indicator_sim=simres,result_code=flag)
    return(res)
  }

#' Generic routine for calculating indicator statistics for non-monitored water body
#' 
#' @param unc_list A list of uncertainty objects from indicator calculations 
#'   All uncertainty objects contain the following variables:
#'   \describe{ 
#'   \item{period}{Mean for the entire assessment period} 
#'   \item{annual}{Annual means for the assessment period} 
#'   \item{indicator_sim}{The simulated distribution of the indicator} 
#'   \item{n_list}{List of dimensions for variance components}} 
#' @param var_list List of variance components containing
#'   \describe{ 
#'   \item{V_WBperiod}{Variance among water bodies for period indicator values} 
#'   \item{V_WBannual}{Variance among water bodies for annual indicator values}} 
#' @param var_WB The variance for the variation among water bodies within a type 
#' @param ntype_WB The total number of waterbodies within the type
#' @param startyear The first year in the indicator calculation
#' @param endyear The last year in the indicator calculation
#' @param n_iter Number of iterations for Monte Carlo simulation
#'   
#' @return An uncertainty object
#' @export
#' 
#' @examples 
CalculateIndicatorType <-
  function(Indicator,unc_list,var_list,ntype_WB,startyear,endyear,n_iter=1000) {
    # Set flag to zero and change it for error handling below
    flag <- 0
    # Calculations on list of uncertainty objects
    n_WB <- length(unc_list)
    n_yearlist <- vector("numeric",n_WB)
    for (i in 1:n_WB) {n_yearlist[i] <- nrow(unc_list[[i]]$annual)}
    n_year <- max(n_yearlist)
    n_WByear <- sum(n_yearlist)
    # Return from function if no information from other WBs is available for calculation
    if (n_WB == 0)  return(list(result_code=-80))
    # Organise data from uncertainty objects into vectors and matrices
    periodWB_mean <- vector("numeric",n_WB)
    periodWB_stderr <- vector("numeric",n_WB)
    annualWB_year <- vector("numeric",n_WByear)
    annualWB_yearmean <- vector("numeric",n_WByear)
    annualWB_yearstderr <- vector("numeric",n_WByear)
    #    annual_unc_list <- matrix(n_year,n_WB)
    icount <- 0
    for (i in 1:n_WB) {
      periodWB_mean[i] <- unc_list[[i]]$period$mean
      periodWB_stderr[i] <- unc_list[[i]]$period$stderr
      for (j in 1:n_yearlist[i]) {
        icount <- icount+1
        annualWB_year[icount] <- unc_list[[i]]$annual$year[j]
        annualWB_yearmean[icount] <- unc_list[[i]]$annual$mean[j]
        annualWB_yearstderr[icount] <- unc_list[[i]]$annual$stderr[j]
      }
    }
    # Find distributions for period and annual means, add variance contribution from variation among WBs
    period_mean <- mean(periodWB_mean)
    period_stderr <- sqrt(stderr_aggr(periodWB_stderr)^2+var_list$V_WBperiod*(1-n_WB/ntype_WB))
    annual_WB <- data.frame(year=annualWB_year,mean=annualWB_yearmean,stderr=annualWB_yearstderr)
    annual <- annual_WB %>% group_by(year) %>% summarise(mean = mean(mean),stderr = sqrt(stderr_aggr(periodWB_stderr)^2+var_list$V_WBannual*(1-n_WB/ntype_WB)))
    # Simulate distribution - the aggregate across WBs is assumed normal distributed
    simres <- rnorm(n_iter,mean=period_mean,sd=period_stderr)
    # calculate information for uncertainty object
    period <- data.frame(mean=period_mean,stderr=period_stderr,
                         lower_1  = qnorm(0.005,period_mean,period_stderr),
                         lower_5  = qnorm(0.025,period_mean,period_stderr),
                         lower_10 = qnorm(0.05,period_mean,period_stderr),
                         upper_10 = qnorm(0.95,period_mean,period_stderr),
                         upper_5  = qnorm(0.975,period_mean,period_stderr),
                         upper_1  = qnorm(0.995,period_mean,period_stderr),row.names = NULL)
    annual <- data.frame(year = annual$year ,mean = annual$mean,stderr = annual$stderr,
                         lower_1  = qnorm(0.005,annual$mean,annual$stderr),
                         lower_5  = qnorm(0.025,annual$mean,annual$stderr),
                         lower_10 = qnorm(0.05,annual$mean,annual$stderr),
                         upper_10 = qnorm(0.95,annual$mean,annual$stderr),
                         upper_5  = qnorm(0.975,annual$mean,annual$stderr),
                         upper_1  = qnorm(0.995,annual$mean,annual$stderr))
    
    res <- list(period=period,annual=annual,indicator_sim=simres,result_code=flag)
    return(res)
  }    

