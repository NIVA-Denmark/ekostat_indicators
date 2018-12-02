#' LIBRARY OF ROUTINES FOR CALCULATING SWEDISH WFD INDICATORS
#'
#'
#' Swedish indicator for Chla
#' 
#' 
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
#' @param MonthInclude A list of months to be included in the indicator
#' @param var_list List of variance components
#' @param n_iter Number of iterations for Monte Carlo simulation
#'   
#' @return
#' @export
#' 
#' @examples
CalculateIndicator_Chla <-
  function(df,MonthInclude,var_list,n_iter=1000) {
# Set flag to zero and change it for error handling below
    flag <- 0
# Filter dataframe to include observations used in indicator only
    df <- Filter_df(df,MonthInclude,"chla")
# Calculate number of years, stations, months, institutions and combinations thereof in df 
    ndf <- DF_Ncalculation(df)
    if (ndf$n_obs==0) {
      flag <- -2
      stop("There are no valid observations for indicator computation!")}
    df <- mutate(df,yvar = chla)
# Estimate mean of the log-transformed chla obs
    alpha <- mean(log(df$chla))
# Simulate system with random variables for estimating the variance of the indicator
    simres <- vector("numeric",n_iter)
    simresyear <- matrix(nrow=n_iter,ncol=ndf$n_year)
# simulation loop - simres contains the residuals from n_iter simulations
    for (isim in 1:n_iter) {
    # simulate variations in the random factors using the data structure
          simulobs <- SetVector_IndicatorSim(alpha,ndf,var_list,df)
    # backtransform simulations from log domain to original domain
          simulobs <- exp(simulobs)
    # add simulated observation to df
          simul_df <- data.frame(df$year,df$station,simulobs)
    # Calculate indicator value for each year
          simresyearx <- simul_df %>%
            group_by(df.year,df.station) %>%
            summarise(simulobs = mean(simulobs)) %>%
            group_by(df.year) %>%
            summarise(simulobs = mean(simulobs))
          simresyear[isim,]=simresyearx$simulobs
    # Calculate indicator for period
          simres[isim] <- mean(simresyearx$simulobs)
         } # end simulation loop
# Estimate statistics from original data and simulations
    res <- Indicator_statistics(df,simres,simresyear,simresyearx$df.year,ndf,flag)
    return(res)
  }

#' Swedish indicator for Chla, based on EQR calculations
#' 
#' 
#' @param df A dataframe with monitoring data from the Swedish Monitoring program. 
#'   The dataframe should contain the
#'   following variables:
#'   
#'   \describe{ 
#'   \item{station}{An identifier for the monitoring station.} 
#'   \item{date}{Date of the sampling.} 
#'   \item{institution}{Provider of the observation.} 
#'   \item{chla}{Chlorophyll a concentration in sample.} 
#'   \item{sali}{Salinity for the measured chla observation.} 
#'   
#' @param RefCond_sali Vector of RC values (n=36) for different salinity levels
#' @param MonthInclude A list of months to be included in the indicator
#' @param var_list List of variance components
#' @param n_iter Number of iterations for Monte Carlo simulation
#'   
#' @return
#' @export
#' 
#' @examples
CalculateIndicator_ChlaEQR <-
  function(df,RefCond_sali,MonthInclude,var_list,n_iter=1000) {
# Set flag to zero and change it for error handling below
    flag <- 0
# Filter dataframe to include observations used in indicator only
    df <- Filter_df(df,MonthInclude,"chla")
# Calculate number of years, stations, months, institutions and combinations thereof in df 
    ndf <- DF_Ncalculation(df)
# setting RefCond depending on salinity and calculate chlaEQR
    RefCond_chla <- mat.or.vec(ndf$n_obs, 1)
    sali_class <- findInterval(df$sali, c(seq(0, 35)))
    for (i in 1:ndf$n_obs) {RefCond_chla[i] <- RefCond_sali[sali_class[i]]}
    df <- mutate(df,RefCond_chla = RefCond_chla)
    df <- mutate(df,chlaEQR = ifelse(chla<RefCond_chla,1,RefCond_chla/chla))
    df <- mutate(df,yvar = chlaEQR)
# Estimate mean of the log-transformed chla obs
    alpha <- mean(log(df$chla))
# Calculate number of years, stations, year*stations, institutions in df 
    ndf <- DF_Ncalculation(df)
# Simulate system with random variables for estimating the variance of the indicator
    simres <- vector("numeric",n_iter)
    simresyear <- matrix(nrow=n_iter,ncol=ndf$n_year)
# simulation loop - simres contains the residuals from n_iter simulations
    for (isim in 1:n_iter) {
    # simulate variations in the random factors using the data structure
      simulobs <- SetVector_IndicatorSim(alpha,ndf,var_list,df)
    # backtransform simulations from log domain to original domain
      simulobs <- exp(simulobs)
    # transform simulations to EQR scale with truncation
      simulobs <- ifelse(simulobs<RefCond_chla,1,RefCond_chla/simulobs)
    # add simulated observation to df
      simul_df <- data.frame(df$year,df$station,simulobs)
    # Calculate indicator value for each year
      simresyearx <- simul_df %>%
        group_by(df.year,df.station) %>%
        summarise(simulobs = mean(simulobs)) %>%
        group_by(df.year) %>%
        summarise(simulobs = mean(simulobs))
      simresyear[isim,]=simresyearx$simulobs
    # Calculate indicator for period
      simres[isim] <- mean(simresyearx$simulobs)
     } # end simulation loop

# Estimate statistics from original data and simulations
    res <- Indicator_statistics(df,simres,simresyear,simresyearx$df.year,ndf,flag)
    return(res)
  }

#' Swedish indicator for nutrient concentrations, based on EQR calculations
#' 
#' 
#' @param NutrientIndicator A string giving the name of the nutrient indicator to be calculated
#' @param df A dataframe with monitoring data from the Swedish Monitoring program. 
#'   The dataframe should contain the
#'   following variables:
#'   
#'   \describe{ 
#'   \item{station}{An identifier for the monitoring station.} 
#'   \item{date}{Date of the sampling.} 
#'   \item{institution}{Provider of the observation.} 
#'   \item{DIN}{DIN concentration in sample.} 
#'   \item{TN}{TN concentration in sample.} 
#'   \item{DIP}{DIP concentration in sample.} 
#'   \item{TP}{TP concentration in sample.} 
#'   \item{sali}{Salinity for the measured nutrient observation.} 
#'   
#' @param RefCond_sali Vector of RC values (n=36) for different salinity levels
#' @param MonthInclude A list of months to be included in the indicator
#' @param var_list List of variance components
#' @param n_iter Number of iterations for Monte Carlo simulation
#'   
#' @return
#' @export
#' 
#' @examples
CalculateIndicator_nutrientEQR <-
  function(NutrientIndicator,df,RefCond_sali,MonthInclude,var_list,n_iter=1000) {
# Set flag to zero and change it for error handling below
    flag <- 0
# Select the indicator response variable
    nutrient <- switch(NutrientIndicator,
                       DINsummer = df$DIN,
                       DIPsummer = df$DIP,
                       TNsummer = df$TN,
                       TNwinter = df$TN,
                       TPsummer = df$TP,
                       TPwinter = df$TP)
    df <- mutate(df,nutrient=nutrient)  
# Filter dataframe to include observations used in indicator only
    df <- Filter_df(df,MonthInclude,"nutrient")
# Switch year for winter months (Nov+Dec) to include in winter indicators
    df <- mutate(df,year=ifelse(month %in% c(11,12),year+1,year))
# Calculate number of years and stations in df 
    ndf <- DF_Ncalculation(df)
# setting RefCond depending on salinity and calculate chlaEQR
    RefCond_nutrient <- mat.or.vec(ndf$n_obs, 1)
    sali_class <- findInterval(df$sali, c(seq(0, 35)))
    for (i in 1:ndf$n_obs) {RefCond_nutrient[i] <- RefCond_sali[sali_class[i]]}
    df <- mutate(df,RefCond_nutrient = RefCond_nutrient)
    df <- mutate(df,nutrientEQR = RefCond_nutrient/nutrient)
    df <- mutate(df,yvar = nutrientEQR)
    # Estimate mean of the log-transformed nutrient obs
    alpha <- mean(log(df$nutrient))
# Simulate system with random variables for estimating the variance of the indicator
    simres <- vector("numeric",n_iter)
    simresyear <- matrix(nrow=n_iter,ncol=ndf$n_year)
# simulation loop - simres contains the residuals from n_iter simulations
    for (isim in 1:n_iter) {
      # simulate variations in the random factors using the data structure
      simulobs <- SetVector_IndicatorSim(alpha,ndf,var_list,df)
      # backtransform simulations from log domain to original domain
      simulobs <- exp(simulobs)
      # transform simulations to EQR scale without truncation
      simulobs <- RefCond_nutrient/simulobs
      # add simulated observation to df
      simul_df <- data.frame(df$year,df$station,simulobs)
      # Calculate indicator value for each year
      simresyearx <- simul_df %>%
        group_by(df.year,df.station) %>%
        summarise(simulobs = mean(simulobs)) %>%
        group_by(df.year) %>%
        summarise(simulobs = mean(simulobs))
      simresyear[isim,]=simresyearx$simulobs
      # Calculate indicator for period
      simres[isim] <- mean(simresyearx$simulobs)
    } # end simulation loop

# Estimate statistics from original data and simulations
    res <- Indicator_statistics(df,simres,simresyear,simresyearx$df.year,ndf,flag)
    return(res)
  }

#' Swedish indicator for Secchi depth
#' 
#' 
#' @param df A dataframe with monitoring data from the Swedish Monitoring program. 
#'   The dataframe should contain the
#'   following variables:
#'   
#'   \describe{ 
#'   \item{station}{An identifier for the monitoring station.} 
#'   \item{date}{Date of the observation.} 
#'   \item{institution}{Provider of the observation.} 
#'   \item{secchi}{Secchi depth at the sampling location.} 
#'   
#' @param MonthInclude A list of months to be included in the indicator
#' @param var_list List of variance components
#' @param n_iter Number of iterations for Monte Carlo simulation
#'   
#' @return
#' @export
#' 
#' @examples
CalculateIndicator_Secchi <-
  function(df,MonthInclude,var_list,n_iter=1000) {
# Set flag to zero and change it for error handling below
    flag <- 0
# Filter dataframe to include observations used in indicator only
    df <- Filter_df(df,MonthInclude,"secchi")
# Calculate number of years, stations, months, institutions and combinations thereof in df 
    ndf <- DF_Ncalculation(df)
    if (ndf$n_obs==0) {
      flag <- -2
      stop("There are no valid observations for indicator computation!")}
    df <- mutate(df,yvar = secchi)
# Estimate mean of the log-transformed chla obs
    alpha <- mean(log(df$secchi))
# Simulate system with random variables for estimating the variance of the indicator
    simres <- vector("numeric",n_iter)
    simresyear <- matrix(nrow=n_iter,ncol=ndf$n_year)
# simulation loop - simres contains the residuals from n_iter simulations
    for (isim in 1:n_iter) {
      # simulate variations in the random factors using the data structure
      simulobs <- SetVector_IndicatorSim(alpha,ndf,var_list,df)
      # backtransform simulations from log domain to original domain
      simulobs <- exp(simulobs)
      # add simulated observation to df
      simul_df <- data.frame(df$year,df$station,simulobs)
      # Calculate indicator value for each year
      simresyearx <- simul_df %>%
        group_by(df.year,df.station) %>%
        summarise(simulobs = mean(simulobs)) %>%
        group_by(df.year) %>%
        summarise(simulobs = mean(simulobs))
      simresyear[isim,]=simresyearx$simulobs
      # Calculate indicator for period
      simres[isim] <- mean(simresyearx$simulobs)
    } # end simulation loop
# Estimate statistics from original data and simulations
    res <- Indicator_statistics(df,simres,simresyear,simresyearx$df.year,ndf,flag)
    return(res)
  }


#' Swedish indicator for Secchi depth, based on EQR calculations
#' 
#' 
#' @param df A dataframe with monitoring data from the Swedish Monitoring program. 
#'   The dataframe should contain the
#'   following variables:
#'   
#'   \describe{ 
#'   \item{station}{An identifier for the monitoring station.} 
#'   \item{date}{Date of the sampling.} 
#'   \item{institution}{Provider of the observation.} 
#'   \item{secchi}{Secchi depth at the sampling location.} 
#'   \item{sali}{Salinity for the measured chla observation.} 
#'   
#' @param RefCond_sali Vector of RC values (n=36) for different salinity levels
#' @param MonthInclude A list of months to be included in the indicator
#' @param var_list List of variance components
#' @param n_iter Number of iterations for Monte Carlo simulation
#'   
#' @return
#' @export
#' 
#' @examples
CalculateIndicator_secchiEQR <-
  function(df,RefCond_sali,MonthInclude,var_list,n_iter=1000) {
# Set flag to zero and change it for error handling below
    flag <- 0
# Filter dataframe to include observations used in indicator only
    df <- Filter_df(df,MonthInclude,"secchi")
# Calculate number of years, stations, months, institutions and combinations thereof in df 
    ndf <- DF_Ncalculation(df)
# setting RefCond depending on salinity and calculate secchiEQR
    RefCond_secchi <- mat.or.vec(ndf$n_obs, 1)
    sali_class <- findInterval(df$sali, c(seq(0, 35)))
    for (i in 1:ndf$n_obs) {RefCond_secchi[i] <- RefCond_sali[sali_class[i]]}
    df <- mutate(df,RefCond_secchi = RefCond_secchi)
    df <- mutate(df,secchiEQR = secchi/RefCond_secchi)
    df <- mutate(df,yvar = secchiEQR)
# Estimate mean of the log-transformed chla obs
    alpha <- mean(log(df$secchi))
# Calculate number of years, stations, year*stations, institutions in df 
    ndf <- DF_Ncalculation(df)
    # Simulate system with random variables for estimating the variance of the indicator
    simres <- vector("numeric",n_iter)
    simresyear <- matrix(nrow=n_iter,ncol=ndf$n_year)
    # simulation loop - simres contains the residuals from n_iter simulations
    for (isim in 1:n_iter) {
      # simulate variations in the random factors using the data structure
      simulobs <- SetVector_IndicatorSim(alpha,ndf,var_list,df)
      # backtransform simulations from log domain to original domain
      simulobs <- exp(simulobs)
      # transform simulations to EQR scale with truncation
      simulobs <- simulobs/RefCond_secchi
      # add simulated observation to df
      simul_df <- data.frame(df$year,df$station,simulobs)
      # Calculate indicator value for each year
      simresyearx <- simul_df %>%
        group_by(df.year,df.station) %>%
        summarise(simulobs = mean(simulobs)) %>%
        group_by(df.year) %>%
        summarise(simulobs = mean(simulobs))
      simresyear[isim,]=simresyearx$simulobs
      # Calculate indicator for period
      simres[isim] <- mean(simresyearx$simulobs)
    } # end simulation loop
    
    # Estimate statistics from original data and simulations
    res <- Indicator_statistics(df,simres,simresyear,simresyearx$df.year,ndf,flag)
    return(res)
  }

