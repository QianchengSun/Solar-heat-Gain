# This is the revised verison of R-script based on:
# "Smart Wi-Fi physics-informed thermostat enabled estimation of residential passive solar heat gain for any residence"
# Introduction:
# By using this package (function) can calculate the hourly solar radiation for the residential houses.
# The output of the dataset will have:
# 1. Solar time
# 2. Solar angles
# 3. GHI 
# 4. Beam solar radiation on each orientation
# 5. Diffuse solar radiation on each orientation
# 5. Ground reflective on each orientation
# More information please check the reference link:
# https://doi.org/10.1016/j.enbuild.2022.111934

# Access the data from NASA Power website: 
# https://power.larc.nasa.gov/data-access-viewer/
# Data explanation:
# ALLSKY_SFC_DWN : ALL Sky Surface Shortwave Downward Irradiance (Wh/m^2) # will be used to calculate the solar radiation under all sky condition
# CLRSKY_SFC_DWN : Clear Sky Surface Shortwave Downward Irradiance (Wh/m^2) # could be used to calculate the solar radiation under clear shy condition
# T2M : Tmeperature at 2 Meters (C) # will be used during the solar radiation calculation
# ALLSKY_KT : All Insolation Clearness Index (dimensionless) # could be used as the cloud cover data. Note: In old version of NASA power website, there is no information about the cloud cover information
# RH2M : Relative Humidity at 2 Meters (%) 
# If you are looking for more features that related daily/hourly solar information, please visit:
# https://power.larc.nasa.gov/data-access-viewer/

# load packages
library(MASS)
library(boot)
library(CircStats)
library(grid) # required package for REdaS
library(REdaS) # use to calculate rad2deg
# define path which accessed for NASA power website
file_path <- "/Users/qianchengsun/Desktop/MEE460 teaching recording/Melink/POWER_Point_Hourly_20210630_20220131_039d3569N_084d2976W_LST.csv"
# read csv solar file
df <- read.csv(file_path)
# define the input information
GMT <- -5 # Greenwhich time for Mason, Ohio
longitude <- -84.2976 # longitude for Mason, Ohio
latitude <- 39.3569 # latitude for Mason, Ohio
H_bar <- df$ALLSKY_SFC_SW_DWN # All Sky Surface Shortwave Downward Irradiance (Wh/m^2) 
## Surface Azimuth Angle degree
SurfaceAzimuthAngle_south <- 0
SurfaceAzimuthAngle_north <- 180
SurfaceAzimuthAngle_west <- 90 
SurfaceAzimuthAngle_east <- -90
### north 180, east -90, south 0, west 90
TiltAngle_deg_west <- 90 # the tilt angle has to be calculated in four different orientation
TiltAngle_deg_north <- 90
TiltAngle_deg_east <- 90
TiltAngle_deg_south <- 90
# define the air temperature at 2 meters
T_air <- df$T2M
# set up the output path 
output_path  <- "/Users/qianchengsun/Desktop/MEE460 teaching recording/Melink/solar_data.csv"

# calculate the day of year 
# After using this function there will generate a warning about the dimension of n, 
# in function define, I simple define n <- NA, I haven't define the shape of n, in order to fix this issue
# future version of solar calculation update will fix this issue,
# however, the result will not be effect by this warning
# Note: On day of year calculation, in February calculate, I simply define it as 28. 
# because when February has 29 days, the total number of day of year will change into 365.
# In future version of solar calculation, this issue will be fixed.
df$day_of_year <- Day_of_year(month = df$MO,
                              day = df$DY)
# calculate the equation of time
df$equation_of_time <- EquationOfTime(n = df$day_of_year)
# calculate the solar time for each time
df$solar_time <- SolarTime(hour = df$HR,
                           GMT = GMT,
                           long_local = longitude,
                           n = df$day_of_year)
# calculate the declination angle based on the day of year
df$declination_angle <- DeclinationAngle(n = df$day_of_year)
# calculate the sunset angle 
df$sunset_angle <- SunsetAngle(n = df$day_of_year,
                               latitude = latitude)
# calculate the hour angle
df$hour_angle <- HourAngle(Solar_Time = df$solar_time)

# calculate the GHI 
df$GHI <- GHI(n = df$day_of_year,
              latitude = latitude,
              H_bar = H_bar,
              Solar_Time = df$solar_time)

# calculate the solar altitude angle
df$solar_altitude_angle <- SolarAltitudeAngle(latitude = latitude,
                                              Declination_Angle = df$declination_angle,
                                              Solar_Time = df$solar_time)
# calculate the solar azimuth angle
df$solar_azimuth_angle <- SolarAzimuthAngle(latitude = latitude,
                                            Declination_Angle = df$declination_angle,
                                            Solar_Time = df$solar_time,
                                            Solar_Altitude_Angle = df$solar_altitude_angle)
# calculate the solar incidence angle on south side
df$solar_incidence_angle_south <- SolarIncidenceAngle_south(SurfaceAzimuthAngle_south = SurfaceAzimuthAngle_south,
                                                            TiltAngle_deg_south = TiltAngle_deg_south,
                                                            Solar_Altitude_Angle = df$solar_altitude_angle,
                                                            Solar_Azimuth_Angle = df$solar_azimuth_angle)

# calculate the solar incidence angle on north side
df$solar_incidence_angle_north <- SolarIncidenceAngle_north(SurfaceAzimuthAngle_north = SurfaceAzimuthAngle_north,
                                                            TiltAngle_deg_north = TiltAngle_deg_north,
                                                            Solar_Altitude_Angle = df$solar_altitude_angle,
                                                            Solar_Azimuth_Angle = df$solar_azimuth_angle)

# calculate the solar incidence angle on east side
df$solar_incidence_angle_east <- SolarIncidenceAngle_east(SurfaceAzimuthAngle_east = SurfaceAzimuthAngle_east,
                                                            TiltAngle_deg_east = TiltAngle_deg_east,
                                                            Solar_Altitude_Angle = df$solar_altitude_angle,
                                                            Solar_Azimuth_Angle = df$solar_azimuth_angle)

# calculate the solar incidence angle on west side
df$solar_incidence_angle_west <- SolarIncidenceAngle_west(SurfaceAzimuthAngle_west = SurfaceAzimuthAngle_west,
                                                          TiltAngle_deg_west = TiltAngle_deg_west,
                                                          Solar_Altitude_Angle = df$solar_altitude_angle,
                                                          Solar_Azimuth_Angle = df$solar_azimuth_angle)

# Calculate the daily extraterrestrial solar radiation in kW/m^2/day
df$daily_extraterrestrial_solar_radiation <- H_bar_0(latitude = latitude,
                                                     n = df$day_of_year)

# Calculate Monthly average daily diffuse solar radiation on horizontal (south)
df$monthly_averge_daily_diffuse_solar_radiation_south <- H_bar_d_south(latitude = latitude,
                                                                 H_bar = H_bar,
                                                                 TiltAngle_deg_south = TiltAngle_deg_south,
                                                                 n = df$day_of_year)

# Calculate Monthly average daily diffuse solar radiation on horizontal (north)
df$monthly_averge_daily_diffuse_solar_radiation_north <- H_bar_d_north(latitude = latitude,
                                                                 H_bar = H_bar,
                                                                 TiltAngle_deg_north = TiltAngle_deg_north,
                                                                 n = df$day_of_year)
# Calculate Monthly average daily diffuse solar radiation on horizontal (east)
df$monthly_averge_daily_diffuse_solar_radiation_east <- H_bar_d_east(latitude = latitude,
                                                                       H_bar = H_bar,
                                                                       TiltAngle_deg_east = TiltAngle_deg_east,
                                                                       n = df$day_of_year)
# Calculate Monthly average daily diffuse solar radiation on horizontal (west)
df$monthly_averge_daily_diffuse_solar_radiation_west <- H_bar_d_west(latitude = latitude,
                                                                       H_bar = H_bar,
                                                                       TiltAngle_deg_west = TiltAngle_deg_west,
                                                                       n = df$day_of_year)
# calculate beam radiation on the south side of the house
df$beam_radiation_south <- I_bar_b_south(n = df$day_of_year,
                                         latitude = latitude,
                                         H_bar = H_bar,
                                         Solar_Time = df$solar_time,
                                         Solar_Altitude_Angle = df$solar_altitude_angle,
                                         Solar_Incidence_Angle_south = df$solar_incidence_angle_south)
# calculate beam radiation on the north side of the house
df$beam_radiation_north <- I_bar_b_north(n = df$day_of_year,
                                         latitude = latitude,
                                         H_bar = H_bar,
                                         Solar_Time = df$solar_time,
                                         Solar_Altitude_Angle = df$solar_altitude_angle,
                                         Solar_Incidence_Angle_north = df$solar_incidence_angle_north)

# calculate beam radiaion son the south side of the house
df$beam_radiation_east <- I_bar_b_east(n = df$day_of_year,
                                        latitude = latitude,
                                         H_bar = H_bar,
                                         Solar_Time = df$solar_time,
                                         Solar_Altitude_Angle = df$solar_altitude_angle,
                                         Solar_Incidence_Angle_east = df$solar_incidence_angle_east)

# calculate beam radiation son the west side of the house
df$beam_radiation_west <- I_bar_b_west(n = df$day_of_year,
                                         latitude = latitude,
                                         H_bar = H_bar,
                                         Solar_Time = df$solar_time,
                                         Solar_Altitude_Angle = df$solar_altitude_angle,
                                         Solar_Incidence_Angle_west = df$solar_incidence_angle_west)

# calculate the diffuse solar radiation on the south side of the house
df$diffuse_radiation_south <- I_bar_d_south(n = df$day_of_year,
                                            latitude = latitude,
                                            H_bar = H_bar,
                                            Solar_Time = df$solar_time,
                                            TiltAngle_deg_south = TiltAngle_deg_south)

# calculate the diffuse solar radiation on the north side of the house
df$diffuse_radiation_north <- I_bar_d_north(n = df$day_of_year,
                                            latitude = latitude,
                                            H_bar = H_bar,
                                            Solar_Time = df$solar_time,
                                            TiltAngle_deg_north = TiltAngle_deg_north)

# calculate the diffuse solar radiation on the east side of the house
df$diffuse_radiation_east <- I_bar_d_east(n = df$day_of_year,
                                            latitude = latitude,
                                            H_bar = H_bar,
                                            Solar_Time = df$solar_time,
                                            TiltAngle_deg_east = TiltAngle_deg_east)
# calculate the diffuse solar radiation on the west side of the house
df$diffuse_radiation_west <- I_bar_d_west(n = df$day_of_year,
                                            latitude = latitude,
                                            H_bar = H_bar,
                                            Solar_Time = df$solar_time,
                                            TiltAngle_deg_west = TiltAngle_deg_west)

# calculate the ground reflected solar radiation on the south side of the house
df$ground_reflected_south <- I_bar_g_south(n = df$day_of_year,
                                           latitude = latitude,
                                           H_bar = H_bar,
                                           TiltAngle_deg_south = TiltAngle_deg_south,
                                           Solar_Time = df$solar_time,
                                           T_air = T_air)

# calculate the ground reflected solar radiation on the north side of the house
df$ground_reflected_north <- I_bar_g_north(n = df$day_of_year,
                                           latitude = latitude,
                                           H_bar = H_bar,
                                           TiltAngle_deg_north = TiltAngle_deg_north,
                                           Solar_Time = df$solar_time,
                                           T_air = T_air)

# calculate the ground reflected solar radiation on the east side of the house
df$ground_reflected_east <- I_bar_g_east(n = df$day_of_year,
                                           latitude = latitude,
                                           H_bar = H_bar,
                                           TiltAngle_deg_east = TiltAngle_deg_east,
                                           Solar_Time = df$solar_time,
                                           T_air = T_air)

# calculate the ground reflected solar radiation on the west side of the house
df$ground_reflected_west <- I_bar_g_west(n = df$day_of_year,
                                           latitude = latitude,
                                           H_bar = H_bar,
                                           TiltAngle_deg_west = TiltAngle_deg_west,
                                           Solar_Time = df$solar_time,
                                           T_air = T_air)
# write df as csv file
write.csv(df, output_path)
###########################################################################################################################################
# solar radiation calculation
# calculate day of year
Day_of_year <- function(month,day){
  # Calculation of the day of year by using month & day 
  # Arguments:
  # Input: 
  # month, numeric 
  #        the input month for the month that you are going to calculate
  # day, numeric
  #        the input day for the day that you are going to calculate
  # Output:
  # day of year, n, numeric
  #        the day of year for specific date, for example, if using Jan, 1, as the input, then day of year will be 1.
  
  # Here in order to avoid the comparison error in the if condition, 
  # add a for loop to fix the issue
  n <- NA # create n as NA
  for (i in 1:length(month)){
    if (month[i] == 1 & TRUE ){
      n[i] <- day[i]
    }else if (month[i] == 2 & TRUE){
      n[i] <- 31 + day[i]
    }else if (month[i] == 3 & TRUE){
      n[i] <- 59 + day[i]
    }else if (month[i] == 4 & TRUE){
      n[i] <- 90 + day[i]
    }else if (month[i] == 5 & TRUE){
      n[i] <- 120 + day[i]
    }else if (month[i] == 6 & TRUE){
      n[i] <- 151 + day[i]
    }else if (month[i] == 7 & TRUE){
      n[i] <- 181 + day[i]
    }else if (month[i] == 8 & TRUE){
      n[i] <- 212 + day[i]
    }else if (month[i] == 9 & TRUE){
      n[i] <- 243 + day[i]
    }else if (month[i] == 10 & TRUE){
      n[i] <- 273 + day[i]
    }else if (month[i] == 11 & TRUE){
      n[i] <- 304 + day[i]
    }else {
      n[i] <- 334 + day[i]
    }
  }
  return(n)
} 

## EquationOfTime(n)
EquationOfTime <- function(n){
  # This function calculates the Equation of Time in minutes from the day number (n)
  # Formula is for B in degrees so must be converted to radius
  # Arguments:
  # Input: 
  # n, numeric
  #   n is the day of the year for specific date.
  
  # b = Application.WorksheetFunction.Radians((n - 1) * 360 / 365)
  b <- rad((n - 1) * 360 / 365)
  Equation_Of_Time <- 229.2 * (7.5E-05 + 0.001868 * cos(b) - 0.032077 * sin(b) - 0.014615 * cos(2 * b) - 0.04089 * sin(2 * b))
  return(Equation_Of_Time)
}

## Solar Time function
# Function to convert the regular time into solar time in order to do the calculation. 
SolarTime <- function(hour,GMT,long_local,n){
  #This function calculates the local solar time (in decimal fractions) given the local standard time.
  #Long_local is the local longitude (meridian) measured in degrees negative (west) and positive (east) of Greenwhich (i.e., -180 > Long_local > 180).
  #GMT is a negative (positive) number equal to the west (east) difference in time between Greenwich Mean Time and the local standard (not daylight savings) time.
  #n is the day of the year (1 < n < 365).
  # Arguments:
  # Input: 
  # hour, numeric
  #       input hour for the time that you are going to covert
  # GMT, numeric
  #       input Greenwhich time zone, for example, the Greenwhich in Dayton, Ohio, GMT = -5
  # Long_local, numeric
  #       local longitude of the location
  # n, numeric
  #       day of the year
  # Output:
  # solar time, numeric
  #       convert the regular time into solar time, for example, the regular time at 12:00, could be 13:00 in the solar time.
  Long_std <- GMT * 15
  #n <- Day_of_year(month,day)
  E <- EquationOfTime(n)
  Long_local <- long_local 
  Solar_Time <- hour + (4 * (Long_local - Long_std) + E) / 60
  return(Solar_Time)
}


## DeclinationAngle
DeclinationAngle <- function(n){
  # Function that used to calculate the declination angle
  # Reference website:
  # https://www.pveducation.org/pvcdrom/properties-of-sunlight/declination-angle
  
  # Arguments:
  # Input:
  # day of year, numeric
  #               day of year for specific date
  # Output:
  # declination angle, numeric
  Declination_Angle <- 23.45 * sin(2 * pi * (284 + n) / 365)
  return(Declination_Angle)
}

## SunsetAngle
SunsetAngle <- function(n, latitude){
  # Function that used to calculate the sunset angle 
  # Arguments:
  # Input:
  # n, numeric
  #    day of year for the date
  # latitude, numeric
  #     latitude of the location
  # Output:
  # sunset_angle, numeric
  #             the result of the sunset angle by using latitude and the day of year
  Delta_deg <- DeclinationAngle(n)
  delta_rad <- rad(Delta_deg)
  Lat_rad <- rad(latitude)
  Sunset_Angle <- deg(acos(-tan(Lat_rad) * tan(delta_rad)))
  return(Sunset_Angle)
}


## HourAngle
HourAngle <- function(Solar_Time){
  # Function used to calculate the hour angle
  # Arguments:
  # Input:
  # solar_time, numeric
  #           solar_time for the hour of the day based on the pervious calculation
  Hour_Angle <- (as.numeric(Solar_Time) - 12) * 15
  return(Hour_Angle)
}


## GHI function
GHI <- function(n,latitude,H_bar,Solar_Time){
  # This function calculates the total hourly solar radiation on horizontal over a specified hour of a given day of a given month
  # This is done with formula from Collares-Pereira and Rabl for global irradiate.
  # Arguments:
  # Input:
  # n = day of the year
  # Latitude = latitude of location in degrees. Positive (negative) is for northern (southern) hemisphere.
  # H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # r_h = ratio of hourly total to daily total global radiation
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <- rad(omega_deg)                          # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  a = 0.409 + 0.5016 * sin(omega_rad - pi / 3)         # constant in formula
  b = 0.6609 - 0.4767 * sin(omega_rad - pi / 3)        # constant in formula
  r_h <- ifelse(HourAngle_rad < -omega_rad | HourAngle_rad > omega_rad,
                0,
                pi / 24 * (a + b * cos(HourAngle_rad)) * (cos(HourAngle_rad) - cos(omega_rad)) / (sin(omega_rad) - omega_rad * cos(omega_rad)))
  I_bar <-  r_h * H_bar
  return(I_bar)
}


## Calculate Solar Altitude Angle
SolarAltitudeAngle <- function(latitude,Declination_Angle,Solar_Time){
  # This function calculates the solar altitude angle (in degrees) as a function of the above angles
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/solar-altitude-angle
  
  # Arguments:
  # Inputs:
  # latitude, numeric
  #           latitude of the location
  # declination_angle, numeric
  #           declination_angle of the sun
  # solar_time, numeric
  #           solar_time of the day
  lat_rad <- rad(latitude)     #Latitude in radius
  delta_rad <- rad(Declination_Angle)  #declination angle in radius
  HourAngle_rad <- rad(HourAngle(Solar_Time)) # hour angle in radius
  # obtain the solar altitude angle
  Solar_Altitude_Angle <- rad2deg(asin(cos(lat_rad)*cos(delta_rad)*cos(HourAngle_rad)+sin(lat_rad)*sin(delta_rad)))
  return(Solar_Altitude_Angle)
}



## Calculate Solar Azimuth Angle
SolarAzimuthAngle <- function(latitude,Declination_Angle,Solar_Time,Solar_Altitude_Angle){
  # Function SolarAzimuthAngle(Latitude, Delta_deg, Solar_Time, SolarAltitudeAngle_deg)
  #'This function calculates the solar azimuth angle (in degrees) as a function of the above angles
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/solar-azimuth-angle
  
  SolarAltitudeAngle_rad <- rad(Solar_Altitude_Angle)
  lat_rad <- rad(latitude)     #Latitude in radians
  delta_rad <- rad(Declination_Angle)  #declination angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time)) # hour angle in radians
  # calculate hte solar azimuth angle
  Solar_Azimuth_Angle <-ifelse(Solar_Time > 12,
                              rad2deg(acos((sin(SolarAltitudeAngle_rad)*sin(lat_rad)-sin(delta_rad))/(cos(SolarAltitudeAngle_rad)*cos(lat_rad)))),
                              rad2deg(acos((sin(SolarAltitudeAngle_rad)*sin(lat_rad)-sin(delta_rad))/(cos(SolarAltitudeAngle_rad)*cos(lat_rad))))*-1)
  return(Solar_Azimuth_Angle)
}

# calculate hte solar incidence angle in different orientation

## Calculate SolarIncidenceAngle_South
SolarIncidenceAngle_south <- function(SurfaceAzimuthAngle_south,TiltAngle_deg_south,Solar_Altitude_Angle,Solar_Azimuth_Angle){
  #'This function calculates the solar incidence angle (in degrees) as a function of the above angles
  # SurfaceAzimuthAngle_deg : the facing direction (angle in degrees) as a function of the surface (0=south, east azimuths are NEGATIVE, west azimuths are POSITVE)
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/solar-incidence-angle
  
  # Arguments:
  # Inputs: 
  # surface_azimuth_angle_south, numeric, degree
  #                             surface azimuth angle from south side 
  # tilt_angle_deg_south, numeric, degree
  #                       Tilt angle from the south side for residential building, in research the tilt angle is defined as 90, because wall is always perpendicular with the ground
  # solar_altitude_angle, numeric, degree
  #                       Please view solar altitude angle calculation function
  # solar_azimuth_angle, numeric, degree
  #                       Please view solar azimuth angle calculation function
  
  # convert deg into rad
  SolarAltitudeAngle_rad <- rad(Solar_Altitude_Angle)
  # calculate the surface solar azimuth angle and convert deg into rad
  SurfaceSolarAzimuth_south_rad <- rad(Solar_Azimuth_Angle - SurfaceAzimuthAngle_south)
  # convert tilt angle deg into rad
  TiltAngle_south_rad <- rad(TiltAngle_deg_south)
  # calculate the solar incidence angle, and convert rad into deg
  Solar_Incidence_Angle_south <- rad2deg(acos(cos(SolarAltitudeAngle_rad) * cos(SurfaceSolarAzimuth_south_rad) * sin(TiltAngle_south_rad) + sin(SolarAltitudeAngle_rad) * cos(TiltAngle_south_rad)))
  Solar_Incidence_Angle_south <-ifelse(Solar_Incidence_Angle_south < 90, Solar_Incidence_Angle_south, 90)
  return(Solar_Incidence_Angle_south)
}




## Calculate SolarIncidenceAngle_North
SolarIncidenceAngle_north <- function(SurfaceAzimuthAngle_north,TiltAngle_deg_north,Solar_Altitude_Angle,Solar_Azimuth_Angle){
  #'This function calculates the solar incidence angle (in degrees) as a function of the above angles
  # SurfaceAzimuthAngle_deg : the facing direction (angle in degrees) as a function of the surface (0=south, east azimuths are NEGATIVE, west azimuths are POSITVE)
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/solar-incidence-angle
  
  # Arguments:
  # Inputs: 
  # surface_azimuth_angle_north, numeric, degree
  #                             surface azimuth angle from north side 
  # tilt_angle_deg_north, numeric, degree
  #                       Tilt angle from the north side for residential building, in research the tilt angle is defined as 90, because wall is always perpendicular with the ground
  # solar_altitude_angle, numeric, degree
  #                       Please view solar altitude angle calculation function
  # solar_azimuth_angle, numeric, degree
  #                       Please view solar azimuth angle calculation function
  
  # convert deg into rad
  SolarAltitudeAngle_rad <- rad(Solar_Altitude_Angle)
  # calculate the surface solar azimuth angle and convert deg into rad
  SurfaceSolarAzimuth_north_rad <- rad(Solar_Azimuth_Angle - SurfaceAzimuthAngle_north)
  # convert tilt angle deg into rad
  TiltAngle_north_rad <- rad(TiltAngle_deg_north)
  # calculate the solar incidence angle, and convert rad into deg
  Solar_Incidence_Angle_north <- rad2deg(acos(cos(SolarAltitudeAngle_rad) * cos(SurfaceSolarAzimuth_north_rad) * sin(TiltAngle_north_rad) + sin(SolarAltitudeAngle_rad) * cos(TiltAngle_north_rad)))
  Solar_Incidence_Angle_north <-ifelse(Solar_Incidence_Angle_north < 90, Solar_Incidence_Angle_north, 90)
  return(Solar_Incidence_Angle_north)
}

## Calculate SolarIncidenceAngle_east
SolarIncidenceAngle_east <- function(SurfaceAzimuthAngle_east,TiltAngle_deg_east,Solar_Altitude_Angle,Solar_Azimuth_Angle){
  #'This function calculates the solar incidence angle (in degrees) as a function of the above angles
  # SurfaceAzimuthAngle_deg : the facing direction (angle in degrees) as a function of the surface (0=south, east azimuths are NEGATIVE, west azimuths are POSITVE)
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/solar-incidence-angle
  
  # Arguments:
  # Inputs: 
  # surface_azimuth_angle_east, numeric, degree
  #                             surface azimuth angle from east side 
  # tilt_angle_deg_east, numeric, degree
  #                       Tilt angle from the east side for residential building, in research the tilt angle is defined as 90, because wall is always perpendicular with the ground
  # solar_altitude_angle, numeric, degree
  #                       Please view solar altitude angle calculation function
  # solar_azimuth_angle, numeric, degree
  #                       Please view solar azimuth angle calculation function
  
  # convert deg into rad
  SolarAltitudeAngle_rad <- rad(Solar_Altitude_Angle)
  # calculate the surface solar azimuth angle and convert deg into rad
  SurfaceSolarAzimuth_east_rad <- rad(Solar_Azimuth_Angle - SurfaceAzimuthAngle_east)
  # convert tilt angle deg into rad
  TiltAngle_east_rad <- rad(TiltAngle_deg_east)
  # calculate the solar incidence angle, and convert rad into deg
  Solar_Incidence_Angle_east <- rad2deg(acos(cos(SolarAltitudeAngle_rad) * cos(SurfaceSolarAzimuth_east_rad) * sin(TiltAngle_east_rad) + sin(SolarAltitudeAngle_rad) * cos(TiltAngle_east_rad)))
  Solar_Incidence_Angle_east <-ifelse(Solar_Incidence_Angle_east < 90, Solar_Incidence_Angle_east, 90)
  return(Solar_Incidence_Angle_east)
}


## Calculate SolarIncidenceAngle_west
SolarIncidenceAngle_west <- function(SurfaceAzimuthAngle_west,TiltAngle_deg_west,Solar_Altitude_Angle,Solar_Azimuth_Angle){
  #'This function calculates the solar incidence angle (in degrees) as a function of the above angles
  # SurfaceAzimuthAngle_deg : the facing direction (angle in degrees) as a function of the surface (0=south, west azimuths are NEGATIVE, west azimuths are POSITVE)
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/solar-incidence-angle
  
  # Arguments:
  # Inputs: 
  # surface_azimuth_angle_west, numeric, degree
  #                             surface azimuth angle from west side 
  # tilt_angle_deg_west, numeric, degree
  #                       Tilt angle from the west side for residential building, in research the tilt angle is defined as 90, because wall is always perpendicular with the ground
  # solar_altitude_angle, numeric, degree
  #                       Please view solar altitude angle calculation function
  # solar_azimuth_angle, numeric, degree
  #                       Please view solar azimuth angle calculation function
  
  # convert deg into rad
  SolarAltitudeAngle_rad <- rad(Solar_Altitude_Angle)
  # calculate the surface solar azimuth angle and convert deg into rad
  SurfaceSolarAzimuth_west_rad <- rad(Solar_Azimuth_Angle - SurfaceAzimuthAngle_west)
  # convert tilt angle deg into rad
  TiltAngle_west_rad <- rad(TiltAngle_deg_west)
  # calculate the solar incidence angle, and convert rad into deg
  Solar_Incidence_Angle_west <- rad2deg(acos(cos(SolarAltitudeAngle_rad) * cos(SurfaceSolarAzimuth_west_rad) * sin(TiltAngle_west_rad) + sin(SolarAltitudeAngle_rad) * cos(TiltAngle_west_rad)))
  Solar_Incidence_Angle_west <-ifelse(Solar_Incidence_Angle_west < 90, Solar_Incidence_Angle_west, 90)
  return(Solar_Incidence_Angle_west)
}



# Calculate the daily extraterrestrial solar radiation in kW/m^2/day
H_bar_0 <- function(latitude,n){
  # Function to calculate daily extraterrestrial Ho incident on a horizontal surface
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/extraterrestrial-radiation
  # Arguments:
  # Input:
  # latitude, numeric
  #         latitude of the calculated location
  # n, numeric
  #         day of the year
  Gsc <- 1367
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  #HourAngle_rad <- rad(HourAngle(SolarTime1))           # hour angle in radians
  latitude_rad <- rad(latitude)
  DeclinationAngle_deg <- DeclinationAngle(n)
  DeclinationAngle_rad <- rad(DeclinationAngle_deg)
  H_bar_0_1 <- 86400 * 2.77778E-07 * Gsc / pi * (1 + 0.033 * cos(2 * pi * n/365)) * (cos(latitude_rad) * cos(DeclinationAngle_rad) * sin(omega_rad) + omega_rad * sin(latitude_rad) * sin(DeclinationAngle_rad))
  return(H_bar_0_1)
}



# Calculate Monthly average daily diffuse solar radiation on horizontal (south)
H_bar_d_south <- function(latitude,H_bar,TiltAngle_deg_south,n){
  # Calculate Monthly average daily diffuse solar radiation on horizontal (south)
  if (latitude<0){
    TiltAngle_deg_south <- -TiltAngle_deg_south
  }
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                    # the sunset angle in radius
  DeclinationAngle_deg <- DeclinationAngle(n)
  DeclinationAngle_rad <- rad(DeclinationAngle_deg)
  TiltAngle_south_rad <- rad(TiltAngle_deg_south)
  H_bar_0 <- H_bar_0(latitude,n)
  k_t <- H_bar / H_bar_0                # the clearness index
  H_bar_diffuse_south  <-ifelse(omega_deg< 81.4,
                           max(0,(1.391-3.56*k_t + 4.189*k_t^2-2.137*k_t^3)*H_bar*(1+cos(TiltAngle_south_rad))/2),
                           max(0,(1.311-3.022*k_t + 3.427*k_t^2-1.821*k_t^3)*H_bar*(1+cos(TiltAngle_south_rad))/2))
  
  return(H_bar_diffuse_south)
}

# Calculate Monthly average daily diffuse solar radiation on horizontal (north)
H_bar_d_north <- function(latitude,H_bar,TiltAngle_deg_north,n){
  # Calculate Monthly average daily diffuse solar radiation on horizontal (north)
  if (latitude<0){
    TiltAngle_deg_north <- -TiltAngle_deg_north
  }
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                    # the sunset angle in radius
  DeclinationAngle_deg <- DeclinationAngle(n)
  DeclinationAngle_rad <- rad(DeclinationAngle_deg)
  TiltAngle_north_rad <- rad(TiltAngle_deg_north)
  H_bar_0 <- H_bar_0(latitude,n)
  k_t <- H_bar / H_bar_0                # the clearness index
  H_bar_diffuse_north  <-ifelse(omega_deg< 81.4,
                                max(0,(1.391-3.56*k_t + 4.189*k_t^2-2.137*k_t^3)*H_bar*(1+cos(TiltAngle_north_rad))/2),
                                max(0,(1.311-3.022*k_t + 3.427*k_t^2-1.821*k_t^3)*H_bar*(1+cos(TiltAngle_north_rad))/2))
  
  return(H_bar_diffuse_north)
}



# Calculate Monthly average daily diffuse solar radiation on horizontal (east)
H_bar_d_east <- function(latitude,H_bar,TiltAngle_deg_east,n){
  # Calculate Monthly average daily diffuse solar radiation on horizontal (east)
  if (latitude<0){
    TiltAngle_deg_east <- -TiltAngle_deg_east
  }
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                    # the sunset angle in radius
  DeclinationAngle_deg <- DeclinationAngle(n)
  DeclinationAngle_rad <- rad(DeclinationAngle_deg)
  TiltAngle_east_rad <- rad(TiltAngle_deg_east)
  H_bar_0 <- H_bar_0(latitude,n)
  k_t <- H_bar / H_bar_0                # the clearness index
  H_bar_diffuse_east  <-ifelse(omega_deg< 81.4,
                                max(0,(1.391-3.56*k_t + 4.189*k_t^2-2.137*k_t^3)*H_bar*(1+cos(TiltAngle_east_rad))/2),
                                max(0,(1.311-3.022*k_t + 3.427*k_t^2-1.821*k_t^3)*H_bar*(1+cos(TiltAngle_east_rad))/2))
  
  return(H_bar_diffuse_east)
}

# Calculate Monthly average daily diffuse solar radiation on horizontal (west)
H_bar_d_west <- function(latitude,H_bar,TiltAngle_deg_west,n){
  # Calculate Monthly average daily diffuse solar radiation on horizontal (west)
  if (latitude<0){
    TiltAngle_deg_west <- -TiltAngle_deg_west
  }
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                    # the sunset angle in radius
  DeclinationAngle_deg <- DeclinationAngle(n)
  DeclinationAngle_rad <- rad(DeclinationAngle_deg)
  TiltAngle_west_rad <- rad(TiltAngle_deg_west)
  H_bar_0 <- H_bar_0(latitude,n)
  k_t <- H_bar / H_bar_0                # the clearness index
  H_bar_diffuse_west  <-ifelse(omega_deg< 81.4,
                               max(0,(1.391-3.56*k_t + 4.189*k_t^2-2.137*k_t^3)*H_bar*(1+cos(TiltAngle_west_rad))/2),
                               max(0,(1.311-3.022*k_t + 3.427*k_t^2-1.821*k_t^3)*H_bar*(1+cos(TiltAngle_west_rad))/2))
  
  return(H_bar_diffuse_west)
}



# Calculate beam solar radiation (south)
I_bar_b_south <- function(n,latitude,H_bar,Solar_Time,Solar_Altitude_Angle,Solar_Incidence_Angle_south){
  # This function calculates the diffuse component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # This is done with formula from Collares-Pereira and Rabl for global irradiance.
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/beam-radiation
  # 'H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # 'r_d = ratio of hourly total to daily total diffuse radiation
  # 'r_b = ratio of beam radiation to that on horizontal
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <- rad(omega_deg)                           # the sunset angle in radius
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radius
  SolarIncidenceAngle_south_rad <- rad(Solar_Incidence_Angle_south) #solar incidence angle in radius
  SolarZeinthAngle_rad <- rad(90 - Solar_Altitude_Angle) # solar zenith angle in radius
  
  r_b <- cos(SolarIncidenceAngle_south_rad)/cos(SolarZeinthAngle_rad)
  r_d <- ifelse(HourAngle_rad< -omega_rad|HourAngle_rad>omega_rad,  # constrain to sun hours only
                0,
                pi/24*(cos(HourAngle_rad)-cos(omega_rad))/(sin(omega_rad)-omega_rad*cos(omega_rad)) )

  H_bar_d_0_south <- H_bar_d_south(latitude,H_bar,0,n) # monthly average daily diffuse solar radiation on horizontal
  I_bar_d_horizontal_south <- r_d*H_bar_d_0_south # diffuse radiation on horizontal
  I_bar_horizontal <- GHI(n,latitude,H_bar,Solar_Time) # total radiation on horizontal
  I_bar_b_horizontal_south <- I_bar_horizontal-I_bar_d_horizontal_south # beam radiation on horizontal
  I_bar_beam_south <- I_bar_b_horizontal_south*r_b
  return(I_bar_beam_south)
}


# Calculate beam solar radiation (north)
I_bar_b_north <- function(n,latitude,H_bar,Solar_Time,Solar_Altitude_Angle,Solar_Incidence_Angle_north){
  # This function calculates the diffuse component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # This is done with formula from Collares-Pereira and Rabl for global irradiance.
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/beam-radiation
  # 'H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # 'r_d = ratio of hourly total to daily total diffuse radiation
  # 'r_b = ratio of beam radiation to that on horizontal
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <- rad(omega_deg)                           # the sunset angle in radius
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radius
  SolarIncidenceAngle_north_rad <- rad(Solar_Incidence_Angle_north) #solar incidence angle in radius
  SolarZeinthAngle_rad <- rad(90 - Solar_Altitude_Angle) # solar zenith angle in radius
  
  r_b <- cos(SolarIncidenceAngle_north_rad)/cos(SolarZeinthAngle_rad)
  r_d <- ifelse(HourAngle_rad< -omega_rad|HourAngle_rad>omega_rad,  # constrain to sun hours only
                0,
                pi/24*(cos(HourAngle_rad)-cos(omega_rad))/(sin(omega_rad)-omega_rad*cos(omega_rad)) )
  
  H_bar_d_0_north <- H_bar_d_north(latitude,H_bar,0,n) # monthly average daily diffuse solar radiation on horizontal
  I_bar_d_horizontal_north <- r_d*H_bar_d_0_north # diffuse radiation on horizontal
  I_bar_horizontal <- GHI(n,latitude,H_bar,Solar_Time) # total radiation on horizontal
  I_bar_b_horizontal_north <- I_bar_horizontal-I_bar_d_horizontal_north # beam radiation on horizontal
  I_bar_beam_north <- I_bar_b_horizontal_north*r_b
  return(I_bar_beam_north)
}

# Calculate beam solar radiation (east)
I_bar_b_east <- function(n,latitude,H_bar,Solar_Time,Solar_Altitude_Angle,Solar_Incidence_Angle_east){
  # This function calculates the diffuse component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # This is done with formula from Collares-Pereira and Rabl for global irradiance.
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/beam-radiation
  # 'H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # 'r_d = ratio of hourly total to daily total diffuse radiation
  # 'r_b = ratio of beam radiation to that on horizontal
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <- rad(omega_deg)                           # the sunset angle in radius
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radius
  SolarIncidenceAngle_east_rad <- rad(Solar_Incidence_Angle_east) #solar incidence angle in radius
  SolarZeinthAngle_rad <- rad(90 - Solar_Altitude_Angle) # solar zenith angle in radius
  
  r_b <- cos(SolarIncidenceAngle_east_rad)/cos(SolarZeinthAngle_rad)
  r_d <- ifelse(HourAngle_rad< -omega_rad|HourAngle_rad>omega_rad,  # constrain to sun hours only
                0,
                pi/24*(cos(HourAngle_rad)-cos(omega_rad))/(sin(omega_rad)-omega_rad*cos(omega_rad)) )
  
  H_bar_d_0_east <- H_bar_d_east(latitude,H_bar,0,n) # monthly average daily diffuse solar radiation on horizontal
  I_bar_d_horizontal_east <- r_d*H_bar_d_0_east # diffuse radiation on horizontal
  I_bar_horizontal <- GHI(n,latitude,H_bar,Solar_Time) # total radiation on horizontal
  I_bar_b_horizontal_east <- I_bar_horizontal-I_bar_d_horizontal_east # beam radiation on horizontal
  I_bar_beam_east <- I_bar_b_horizontal_east*r_b
  return(I_bar_beam_east)
}


# Calculate beam solar radiation (west)
I_bar_b_west <- function(n,latitude,H_bar,Solar_Time,Solar_Altitude_Angle,Solar_Incidence_Angle_west){
  # This function calculates the diffuse component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # This is done with formula from Collares-Pereira and Rabl for global irradiance.
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/beam-radiation
  # 'H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # 'r_d = ratio of hourly total to daily total diffuse radiation
  # 'r_b = ratio of beam radiation to that on horizontal
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <- rad(omega_deg)                           # the sunset angle in radius
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radius
  SolarIncidenceAngle_west_rad <- rad(Solar_Incidence_Angle_west) #solar incidence angle in radius
  SolarZeinthAngle_rad <- rad(90 - Solar_Altitude_Angle) # solar zenith angle in radius
  
  r_b <- cos(SolarIncidenceAngle_west_rad)/cos(SolarZeinthAngle_rad)
  r_d <- ifelse(HourAngle_rad< -omega_rad|HourAngle_rad>omega_rad,  # constrain to sun hours only
                0,
                pi/24*(cos(HourAngle_rad)-cos(omega_rad))/(sin(omega_rad)-omega_rad*cos(omega_rad)) )
  
  H_bar_d_0_west <- H_bar_d_west(latitude,H_bar,0,n) # monthly average daily diffuse solar radiation on horizontal
  I_bar_d_horizontal_west <- r_d*H_bar_d_0_west # diffuse radiation on horizontal
  I_bar_horizontal <- GHI(n,latitude,H_bar,Solar_Time) # total radiation on horizontal
  I_bar_b_horizontal_west <- I_bar_horizontal-I_bar_d_horizontal_west # beam radiation on horizontal
  I_bar_beam_west <- I_bar_b_horizontal_west*r_b
  return(I_bar_beam_west)
}



## calculate diffuse solar radiation (south)
I_bar_d_south <- function(n,latitude,H_bar,Solar_Time,TiltAngle_deg_south){
  # 'This function calculates the diffuse component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # 'This is done with formulae from Collares-Pereira and Rabl for global irradiance.
  # Reference website:
  # https://www.sciencedirect.com/topics/earth-and-planetary-sciences/diffuse-radiation
  
  # 'Latitude = latitude of location in degrees.
  # 'H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # 'TiltAngle_deg = tilt angle (in degrees) of the surface relative to the horizontal.
  # 'r_d = ratio of hourly total to daily total diffuse radiation
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  TiltAngle_south_rad <- rad(TiltAngle_deg_south)                  # tilt angle in radians

  r_d <- ifelse(HourAngle_rad < -omega_rad|HourAngle_rad > omega_rad,
                0,
                pi/24*(cos(HourAngle_rad)-cos(omega_rad))/(sin(omega_rad)-omega_rad*cos(omega_rad)))
 
  H_bar_d_0 <- H_bar_d_south(latitude,H_bar,0,n) # monthly average daily diffuse solar radiation on horizontal
  I_bar_d_horizontal <- r_d* H_bar_d_0 # diffuse radiation on horizontal
  I_bar_diffuse_south <- I_bar_d_horizontal*(1+cos(TiltAngle_south_rad))/2 # second term is the view factor
  return(I_bar_diffuse_south)
}



## calculate diffuse solar radiation (north)
I_bar_d_north <- function(n,latitude,H_bar,Solar_Time,TiltAngle_deg_north){
  # 'This function calculates the diffuse component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # 'This is done with formulae from Collares-Pereira and Rabl for global irradiance.
  # Reference website:
  # https://www.sciencedirect.com/topics/earth-and-planetary-sciences/diffuse-radiation
  
  # 'Latitude = latitude of location in degrees. 
  # 'H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # 'TiltAngle_deg = tilt angle (in degrees) of the surface relative to the horizontal.
  # 'r_d = ratio of hourly total to daily total diffuse radiation
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  TiltAngle_north_rad <- rad(TiltAngle_deg_north)                  # tilt angle in radians
  
  r_d <- ifelse(HourAngle_rad < -omega_rad|HourAngle_rad > omega_rad,
                0,
                pi/24*(cos(HourAngle_rad)-cos(omega_rad))/(sin(omega_rad)-omega_rad*cos(omega_rad)))
  
  H_bar_d_0 <- H_bar_d_north(latitude,H_bar,0,n) # monthly average daily diffuse solar radiation on horizontal
  I_bar_d_horizontal <- r_d* H_bar_d_0 # diffuse radiation on horizontal
  I_bar_diffuse_north <- I_bar_d_horizontal*(1+cos(TiltAngle_north_rad))/2 # second term is the view factor
  return(I_bar_diffuse_north)
}


## calculate diffuse solar radiation (east)
I_bar_d_east <- function(n,latitude,H_bar,Solar_Time,TiltAngle_deg_east){
  # 'This function calculates the diffuse component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # 'This is done with formulae from Collares-Pereira and Rabl for global irradiance.
  # Reference website:
  # https://www.sciencedirect.com/topics/earth-and-planetary-sciences/diffuse-radiation
  
  # 'Latitude = latitude of location in degrees. .
  # 'H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # 'TiltAngle_deg = tilt angle (in degrees) of the surface relative to the horizontal.
  # 'r_d = ratio of hourly total to daily total diffuse radiation
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  TiltAngle_east_rad <- rad(TiltAngle_deg_east)                  # tilt angle in radians
  
  r_d <- ifelse(HourAngle_rad < -omega_rad|HourAngle_rad > omega_rad,
                0,
                pi/24*(cos(HourAngle_rad)-cos(omega_rad))/(sin(omega_rad)-omega_rad*cos(omega_rad)))
  
  H_bar_d_0 <- H_bar_d_east(latitude,H_bar,0,n) # monthly average daily diffuse solar radiation on horizontal
  I_bar_d_horizontal <- r_d* H_bar_d_0 # diffuse radiation on horizontal
  I_bar_diffuse_east <- I_bar_d_horizontal*(1+cos(TiltAngle_east_rad))/2 # second term is the view factor
  return(I_bar_diffuse_east)
}

## calculate diffuse solar radiation (west)
I_bar_d_west <- function(n,latitude,H_bar,Solar_Time,TiltAngle_deg_west){
  # 'This function calculates the diffuse component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # 'This is done with formulae from Collares-Pereira and Rabl for global irradiance.
  # Reference website:
  # https://www.sciencedirect.com/topics/earth-and-planetary-sciences/diffuse-radiation
  
  # 'Latitude = latitude of location in degrees.
  # 'H_bar = monthly average daily solar radiation on horizontal (kWh/m^2/day)
  # 'TiltAngle_deg = tilt angle (in degrees) of the surface relative to the horizontal.
  # 'r_d = ratio of hourly total to daily total diffuse radiation
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  TiltAngle_west_rad <- rad(TiltAngle_deg_west)                  # tilt angle in radians
  
  r_d <- ifelse(HourAngle_rad < -omega_rad|HourAngle_rad > omega_rad,
                0,
                pi/24*(cos(HourAngle_rad)-cos(omega_rad))/(sin(omega_rad)-omega_rad*cos(omega_rad)))
  
  H_bar_d_0 <- H_bar_d_west(latitude,H_bar,0,n) # monthly average daily diffuse solar radiation on horizontal
  I_bar_d_horizontal <- r_d* H_bar_d_0 # diffuse radiation on horizontal
  I_bar_diffuse_west <- I_bar_d_horizontal*(1+cos(TiltAngle_west_rad))/2 # second term is the view factor
  return(I_bar_diffuse_west)
}



# Calculate ground reflected solar radiation (south)
I_bar_g_south <- function(n,latitude,H_bar,TiltAngle_deg_south,Solar_Time,T_air){
  # 'This function calculates the ground-reflected component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # 'This is done with formulae from Collares-Pereira and Rabl for global irradiance.
  
  # 'Latitude = latitude of location in degrees.
  # 'H_bar = hourly solar radiation on horizontal (kWh/m^2/day)
  # 'TiltAngle_deg = tilt angle (in degrees) of the surface relative to the horizontal.
  # 'T_air = hourly air temperature (C)
  
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/ground-reflected-radiation
  
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  TiltAngle_south_rad <- rad(TiltAngle_deg_south)                  # tilt angle in radians
  I_bar_horizontal <- GHI(n,latitude,H_bar,Solar_Time) # total radiation on horizontal
  # the following are approximations for ground reflectivity (rho_g) as a function of air temperature
  rho_g <- ifelse(T_air>0,0.2,ifelse(T_air < -5,0.7,-0.1*T_air+0.2))
  I_bar_ground_south <- I_bar_horizontal*rho_g*(1-cos(TiltAngle_south_rad))/2
  return(I_bar_ground_south)
}

# Calculate ground reflected solar radiation (north)
I_bar_g_north <- function(n, latitude,H_bar,TiltAngle_deg_north,Solar_Time,T_air){
  # 'This function calculates the ground-reflected component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # 'This is done with formulae from Collares-Pereira and Rabl for global irradiance.
  
  # 'Latitude = latitude of location in degrees.
  # 'H_bar = hourly solar radiation on horizontal (kWh/m^2/day)
  # 'TiltAngle_deg = tilt angle (in degrees) of the surface relative to the horizontal.
  # 'T_air = hourly air temperature (C)
  
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/ground-reflected-radiation
  
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  TiltAngle_north_rad <- rad(TiltAngle_deg_north)                  # tilt angle in radians
  I_bar_horizontal <- GHI(n,latitude,H_bar,Solar_Time) # total radiation on horizontal
  # the following are approximations for ground reflectivity (rho_g) as a function of air temperature
  rho_g <- ifelse(T_air>0,0.2,ifelse(T_air < -5,0.7,-0.1*T_air+0.2))
  I_bar_ground_north <- I_bar_horizontal*rho_g*(1-cos(TiltAngle_north_rad))/2
  return(I_bar_ground_north)
}


# Calculate ground reflected solar radiation (east)
I_bar_g_east <- function(n, latitude,H_bar,TiltAngle_deg_east,Solar_Time,T_air){
  # 'This function calculates the ground-reflected component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # 'This is done with formulae from Collares-Pereira and Rabl for global irradiance.
  
  # 'Latitude = latitude of location in degrees.
  # 'H_bar = hourly solar radiation on horizontal (kWh/m^2/day)
  # 'TiltAngle_deg = tilt angle (in degrees) of the surface relative to the horizontal.
  # 'T_air = hourly air temperature (C)
  
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/ground-reflected-radiation
  
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  TiltAngle_east_rad <- rad(TiltAngle_deg_east)                  # tilt angle in radians
  I_bar_horizontal <- GHI(n,latitude,H_bar,Solar_Time) # total radiation on horizontal
  # the following are approximations for ground reflectivity (rho_g) as a function of air temperature
  rho_g <- ifelse(T_air>0,0.2,ifelse(T_air < -5,0.7,-0.1*T_air+0.2))
  I_bar_ground_east <- I_bar_horizontal*rho_g*(1-cos(TiltAngle_east_rad))/2
  return(I_bar_ground_east)
}


# Calculate ground reflected solar radiation (west)
I_bar_g_west <- function(n, latitude,H_bar,TiltAngle_deg_west,Solar_Time,T_air){
  # 'This function calculates the ground-reflected component of hourly solar radiation on a surface over a specified hour of a given day of a given month
  # 'This is done with formulae from Collares-Pereira and Rabl for global irradiance.
  
  # 'Latitude = latitude of location in degrees.
  # 'H_bar = hourly solar radiation on horizontal (kWh/m^2/day)
  # 'TiltAngle_deg = tilt angle (in degrees) of the surface relative to the horizontal.
  # 'T_air = hourly air temperature (C)
  
  # Reference website:
  # https://www.sciencedirect.com/topics/engineering/ground-reflected-radiation
  
  omega_deg <- SunsetAngle(n, latitude)     # the sunset angle in degrees
  omega_rad <-rad(omega_deg)                           # the sunset angle in radians
  HourAngle_rad <- rad(HourAngle(Solar_Time))           # hour angle in radians
  TiltAngle_west_rad <- rad(TiltAngle_deg_west)                  # tilt angle in radians
  I_bar_horizontal <- GHI(n,latitude,H_bar,Solar_Time) # total radiation on horizontal
  # the following are approximations for ground reflectivity (rho_g) as a function of air temperature
  rho_g <- ifelse(T_air>0,0.2,ifelse(T_air < -5,0.7,-0.1*T_air+0.2))
  I_bar_ground_west <- I_bar_horizontal*rho_g*(1-cos(TiltAngle_west_rad))/2
  return(I_bar_ground_west)
}
