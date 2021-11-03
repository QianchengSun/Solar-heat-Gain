#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 19:52:38 2021

@author: qianchengsun
"""

import numpy as np
import numpy as np

#%%
# Day of the year
def day_of_year(month, day):
    if month == 1:
        n = day
    elif month == 2 :
        n = 31 + day
    elif month == 3 :
        n =  59 + day
    elif month == 4 :
        n = 90 + day 
    elif month == 5 :
        n = 120 + day 
    elif month == 6 :
        n = 151 + day 
    elif month == 7 :
        n = 181 + day 
    elif month == 8 :
        n = 212 + day
    elif month == 9 :
        n = 243 + day
    elif month == 10 :
        n = 273 + day
    elif month == 11 :
        n = 304 + day
    elif month == 12 :
        n = 304 + day
    else:
        n = 334 + day
    return n 

# Equation of time 
def equation_of_time(n):
    b =  np.radians((n - 1) * 360 / 365) # have to import math then use math.radians
    Equation_of_Time = 229.2 * (7.5 * 10 ** -5 + 0.001868 * np.cos(b) - 0.014615 * np.cos(2 * b) - 0.04089 * np.sin(2 * b))
    return Equation_of_Time

# Solar Time

def Solar_time(hour, GMT, long_local, n):
    Long_std = GMT * 15
    E = equation_of_time(n)
    Long_local = long_local
    SolarTime = hour + (4 * (Long_local - Long_std) + E) / 60
    return SolarTime

# Declination Angle
def Declination_Angle(n):
    DeclinationAngle = 23.45 * np.sin(2 * np.pi * (284 + n) / 365)
    return DeclinationAngle

# Sunset Angle
def Sunset_Angle(month, day, latitude):
    n = day_of_year(month, day)
    Delta_deg = Declination_Angle(n)
    delta_rad = np.radians(Delta_deg)
    Lat_rad = np.radians(latitude)
    SunsetAngle = np.degrees(np.acos(-np.tan(Lat_rad) * np.tan(delta_rad)))
    return SunsetAngle

# Hour Angle 
def Hour_angle(SolarTime):
    HourAngle = (SolarTime - 12) * 15
    return HourAngle

# Radiation on Horizontal (GHI)
def Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime):
    omega_deg = Sunset_Angle(month, day, latitude)  # the sunset angle in degrees
    omega_rad = np.radians(omega_deg)             # the sunset angle in radians
    Hour_angle_rad = np.radians(Hour_angle(SolarTime))   # hour angle in radians
    a = 0.409 + 0.5016 * np.sin(omega_rad - np.pi/3)      # constant in formula
    b = 0.6609 - 0.4767 * np.sin(omega_rad - np.pi/3)     # constant in formula
    if Hour_angle_rad < omega_rad or Hour_angle_rad > omega_rad:
        r_h = 0
    else: 
        r_h = np.pi / 24 * (a + b * np.cos(Hour_angle_rad)) * (np.cos(Hour_angle_rad) - np.cos(omega_rad))/(np.sin(omega_rad) - omega_rad * np.cos(omega_rad))
    I_bar = r_h * H_bar
    return I_bar

# Solar Alitude Angle
def Solar_alitude_angle(latitude, DeclinationAngle, SolarTime):
    lat_rad = np.radians(latitude)
    delta_rad = np.radians(DeclinationAngle)
    Hour_angle_rad = np.radians(SolarTime)
    SolarAltitudeAngle = np.rad2deg(np.asin(np.cos(lat_rad) * np.cos(Hour_angle_rad) + np.sin(lat_rad) * np.sin(delta_rad)))
    return SolarAltitudeAngle

# Solar Azimuth Angle
def Solar_azimuth_angle(latitude, DeclinationAngle, SolarTime, SolarAltitudeAngle):
    Solar_alitude_angle_rad = np.radians(SolarAltitudeAngle)
    lat_rad = np.radians(latitude)
    delta_rad = np.radians(DeclinationAngle)
    if SolarTime > 12:
        SolarAzimuthAngle = np.rad2deg(np.acos((np.sin(Solar_alitude_angle_rad) * np.sin(lat_rad) - np.sin(delta_rad)) / (np.cos(Solar_alitude_angle_rad) * np.cos(lat_rad))))
    else:
        SolarAzimuthAngle = np.rad2deg(np.acos((np.sin(Solar_alitude_angle_rad) * np.sin(lat_rad) - np.sin(delta_rad)) / (np.cos(Solar_alitude_angle_rad) * np.cos(lat_rad))))* -1
    return SolarAzimuthAngle 
#%%

# Surface Azimuth Angle define
SurfaceAzimuthAngle_south = 0
SurfaceAzimuthAngle_north = 180
SurfaceAzimuthAngle_west = 90
SurfaceAzimuthAngle_east = -90

# Solar Incidence Angle South 
def Solar_incidence_angle_south(SurfaceAzimuthAngle_south, TiltAngle_deg_south, SolarAltitudeAngle_deg, SolarAzimuthAngle, n):
    SolarAltitudeAngle_rad = np.radians(SolarAltitudeAngle_deg)
    SurfaceSolarAzimuth_south_rad = np.radians(SolarAzimuthAngle - SurfaceAzimuthAngle_south)
    TiltAngle_south_rad = np.radians(TiltAngle_deg_south)
    SolarIncidenceAngle_south = np.rad2deg(np.acos(SolarAltitudeAngle_rad * np.cos(SurfaceSolarAzimuth_south_rad) * np.sin(TiltAngle_south_rad) + np.sin(SolarAltitudeAngle_rad) * np.cos(TiltAngle_south_rad)))
    if SolarIncidenceAngle_south < 90:
        SolarIncidenceAngle_south = SolarIncidenceAngle_south
    else:
        SolarIncidenceAngle_south = 90
    return SolarIncidenceAngle_south

# Solar Incidence Angle North
def Solar_incidence_angle_north(SurfaceAzimuthAngle_north, TiltAngle_deg_north, SolarAltitudeAngle_deg, SolarAzimuthAngle, n):
    SolarAltitudeAngle_rad = np.radians(SolarAltitudeAngle_deg)
    SurfaceSolarAzimuth_north_rad = np.radians(SolarAzimuthAngle - SurfaceAzimuthAngle_north)
    TiltAngle_north_rad = np.radians(TiltAngle_deg_north)
    SolarIncidenceAngle_north = np.rad2deg(np.acos(SolarAltitudeAngle_rad * np.cos(SurfaceSolarAzimuth_north_rad) * np.sin(TiltAngle_north_rad) + np.sin(SolarAltitudeAngle_rad) * np.cos(TiltAngle_north_rad)))
    if SolarIncidenceAngle_north < 90:
        SolarIncidenceAngle_north = SolarIncidenceAngle_north
    else:
        SolarIncidenceAngle_north = 90
    return SolarIncidenceAngle_north

# Solar Incidence Angle East
def Solar_incidence_angle_east(SurfaceAzimuthAngle_east, TiltAngle_deg_east, SolarAltitudeAngle_deg, SolarAzimuthAngle, n):
    SolarAltitudeAngle_rad = np.radians(SolarAltitudeAngle_deg)
    SurfaceSolarAzimuth_east_rad = np.radians(SolarAzimuthAngle - SurfaceAzimuthAngle_east)
    TiltAngle_east_rad = np.radians(TiltAngle_deg_east)
    SolarIncidenceAngle_east = np.rad2deg(np.acos(SolarAltitudeAngle_rad * np.cos(SurfaceSolarAzimuth_east_rad) * np.sin(TiltAngle_east_rad) + np.sin(SolarAltitudeAngle_rad) * np.cos(TiltAngle_east_rad)))
    if SolarIncidenceAngle_east < 90:
        SolarIncidenceAngle_east = SolarIncidenceAngle_east
    else:
        SolarIncidenceAngle_east = 90
    return SolarIncidenceAngle_east

# Solar Incidence Angle West
def Solar_incidence_angle_west(SurfaceAzimuthAngle_west, TiltAngle_deg_west, SolarAltitudeAngle_deg, SolarAzimuthAngle, n):
    SolarAltitudeAngle_rad = np.radians(SolarAltitudeAngle_deg)
    SurfaceSolarAzimuth_west_rad = np.radians(SolarAzimuthAngle - SurfaceAzimuthAngle_west)
    TiltAngle_west_rad = np.radians(TiltAngle_deg_west)
    SolarIncidenceAngle_west = np.rad2deg(np.acos(SolarAltitudeAngle_rad * np.cos(SurfaceSolarAzimuth_west_rad) * np.sin(TiltAngle_west_rad) + np.sin(SolarAltitudeAngle_rad) * np.cos(TiltAngle_west_rad)))
    if SolarIncidenceAngle_west < 90:
        SolarIncidenceAngle_west = SolarIncidenceAngle_west
    else:
        SolarIncidenceAngle_west = 90
    return SolarIncidenceAngle_west

# Calculate Daily extraterrestrial solar radiation in kw/m^2/day
def H_bar_0(month, day, latitude, n):
    Gsc = 1367
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = np.radians(omega_deg)
    latitude_rad = np.radians(latitude)
    Declination_Angle_deg = Declination_Angle(n)
    Declination_Angle_rad = np.radians(Declination_Angle_deg)
    Hbar_0 = 86400 * 2.77778 * (10 ** -7) * Gsc / np.pi * (1 + 0.033 * np.cos(2 * np.pi * n / 365)) * (np.cos(latitude_rad) * np.cos(Declination_Angle_rad) * np.sin(omega_rad) + omega_rad * np.sin(latitude_rad) * np.sin(Declination_Angle_rad))
    return Hbar_0

 # Calculate Monthly average daily diffuse solar radiation on horizontal (south)
def H_bar_d_south(month, day, latitude, H_bar, TiltAngle_deg_south, n):
    if latitude < 0 :
        TiltAngle_deg_south = -TiltAngle_deg_south
    omega_deg = Sunset_Angle(month, day, latitude)
   # omega_rad = math.radians(omega_deg)
   # Declination_Angle_deg = Declination_Angle(n)
   # Declination_Angle_rad = math.radians(Declination_Angle_deg)
    TiltAngle_south_rad = np.radians(TiltAngle_deg_south)
    Hbar_0 = H_bar_0(month, day, latitude, n)
    k_t = H_bar / Hbar_0
    if omega_deg < 81.4 :
        H_bar_d_south = max(0, (1.391 - 3.56 * k_t + 4.189 * k_t ** 2 - 2.137 * k_t ** 3) * H_bar * (1 + np.cos(TiltAngle_south_rad))/2 )
    else:
        H_bar_d_south = max(0, (1.311 - 3.022 * k_t + 3.427 * k_t ** 2 - 1.821 * k_t ** 3) * H_bar * ( 1 + np.cos(TiltAngle_south_rad))/2 )
    return H_bar_d_south

# Calculate Monthly average daily diffuse solar radiation on horizontal (north)
def H_bar_d_north(month, day, latitude, H_bar, TiltAngle_deg_north, n):
    if latitude < 0 :
        TiltAngle_deg_north = -TiltAngle_deg_north
    omega_deg = Sunset_Angle(month, day, latitude)
   # omega_rad = math.radians(omega_deg)
   # Declination_Angle_deg = Declination_Angle(n)
   # Declination_Angle_rad = math.radians(Declination_Angle_deg)
    TiltAngle_north_rad = np.radians(TiltAngle_deg_north)
    Hbar_0 = H_bar_0(month, day, latitude, n)
    k_t = H_bar / Hbar_0
    if omega_deg < 81.4 :
        H_bar_d_north = max(0, (1.391 - 3.56 * k_t + 4.189 * k_t ** 2 - 2.137 * k_t ** 3) * H_bar * (1 + np.cos(TiltAngle_north_rad))/2 )
    else:
        H_bar_d_north = max(0, (1.311 - 3.022 * k_t + 3.427 * k_t ** 2 - 1.821 * k_t ** 3) * H_bar * ( 1 + np.cos(TiltAngle_north_rad))/2 )
    return H_bar_d_north

# Calculate Monthly average daily diffuse solar radiation on horizontal (east)
def H_bar_d_east(month, day, latitude, H_bar, TiltAngle_deg_east, n):
    if latitude < 0 :
        TiltAngle_deg_east = -TiltAngle_deg_east
    omega_deg = Sunset_Angle(month, day, latitude)
   # omega_rad = math.radians(omega_deg)
   # Declination_Angle_deg = Declination_Angle(n)
   # Declination_Angle_rad = math.radians(Declination_Angle_deg)
    TiltAngle_east_rad = np.radians(TiltAngle_deg_east)
    Hbar_0 = H_bar_0(month, day, latitude, n)
    k_t = H_bar / Hbar_0
    if omega_deg < 81.4 :
        H_bar_d_east = max(0, (1.391 - 3.56 * k_t + 4.189 * k_t ** 2 - 2.137 * k_t ** 3) * H_bar * (1 + np.cos(TiltAngle_east_rad))/2 )
    else:
        H_bar_d_east = max(0, (1.311 - 3.022 * k_t + 3.427 * k_t ** 2 - 1.821 * k_t ** 3) * H_bar * ( 1 + np.cos(TiltAngle_east_rad))/2 )
    return H_bar_d_east

# Calculate Monthly average daily diffuse solar radiation on horizontal (west)
def H_bar_d_west(month, day, latitude, H_bar, TiltAngle_deg_west, n):
    if latitude < 0 :
        TiltAngle_deg_west = -TiltAngle_deg_west
    omega_deg = Sunset_Angle(month, day, latitude)
   # omega_rad = math.radians(omega_deg)
   # Declination_Angle_deg = Declination_Angle(n)
   # Declination_Angle_rad = math.radians(Declination_Angle_deg)
    TiltAngle_west_rad = np.radians(TiltAngle_deg_west)
    Hbar_0 = H_bar_0(month, day, latitude, n)
    k_t = H_bar / Hbar_0
    if omega_deg < 81.4 :
        H_bar_d_west = max(0, (1.391 - 3.56 * k_t + 4.189 * k_t ** 2 - 2.137 * k_t ** 3) * H_bar * (1 + np.cos(TiltAngle_west_rad))/2 )
    else:
        H_bar_d_west = max(0, (1.311 - 3.022 * k_t + 3.427 * k_t ** 2 - 1.821 * k_t ** 3) * H_bar * ( 1 + np.cos(TiltAngle_west_rad))/2 )
    return H_bar_d_west

#%%
#===================================
# Beam Solar Radiation Calculation
#===================================

# Calculate beam solar radiation (south)
def Beam_solar_radiation_south(month, day, latitude, H_bar, SolarTime, SolarAltitudeAngle, SolarIncidenceAngle_south,TiltAngle_deg_south,n):
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = np.radians(omega_deg)
    HourAngle_rad = np.radians(Hour_angle(SolarTime))
    Solar_Incidence_Angle_south_rad = np.radians(SolarIncidenceAngle_south)
    SolarZeinthAngle_rad = np.radians(90 - SolarAltitudeAngle)
    
    r_b = np.cos(Solar_Incidence_Angle_south_rad)/np.cos(SolarZeinthAngle_rad)
    if HourAngle_rad < - omega_rad or HourAngle_rad > omega_rad:
        r_d = 0
    else:
        r_d = np.pi/24 * (np.cos(HourAngle_rad) - np.cos(omega_rad)) / (np.sin(omega_rad) - omega_rad * np.cos(omega_rad))
    
    H_bar_d_0_south = H_bar_d_south(month, day, latitude, H_bar, TiltAngle_deg_south, n)
    I_bar_d_horizontal_south = r_d * H_bar_d_0_south
    I_bar_horizontal = Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime)
    I_bar_b_horizontal_south = I_bar_horizontal - I_bar_d_horizontal_south
    I_bar_b_south = I_bar_b_horizontal_south * r_b
    return I_bar_b_south


# Calculate beam solar radiation (north)
def Beam_solar_radiation_north(month, day, latitude, H_bar, SolarTime, SolarAltitudeAngle, SolarIncidenceAngle_north,TiltAngle_deg_north,n):
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = np.radians(omega_deg)
    HourAngle_rad = np.radians(Hour_angle(SolarTime))
    Solar_Incidence_Angle_north_rad = np.radians(SolarIncidenceAngle_north)
    SolarZeinthAngle_rad = np.radians(90 - SolarAltitudeAngle)
    
    r_b = np.cos(Solar_Incidence_Angle_north_rad)/np.cos(SolarZeinthAngle_rad)
    if HourAngle_rad < - omega_rad or HourAngle_rad > omega_rad:
        r_d = 0
    else:
        r_d = np.pi/24 * (np.cos(HourAngle_rad) - np.cos(omega_rad)) / (np.sin(omega_rad) - omega_rad * np.cos(omega_rad))
    
    H_bar_d_0_north = H_bar_d_north(month, day, latitude, H_bar, TiltAngle_deg_north, n)
    I_bar_d_horizontal_north = r_d * H_bar_d_0_north
    I_bar_horizontal = Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime)
    I_bar_b_horizontal_north = I_bar_horizontal - I_bar_d_horizontal_north
    I_bar_b_north = I_bar_b_horizontal_north * r_b
    return I_bar_b_north


# Calculate beam solar radiation (east)
def Beam_solar_radiation_east(month, day, latitude, H_bar, SolarTime, SolarAltitudeAngle, SolarIncidenceAngle_east,TiltAngle_deg_east,n):
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = np.radians(omega_deg)
    HourAngle_rad = np.radians(Hour_angle(SolarTime))
    Solar_Incidence_Angle_east_rad = np.radians(SolarIncidenceAngle_east)
    SolarZeinthAngle_rad = np.radians(90 - SolarAltitudeAngle)
    
    r_b = np.cos(Solar_Incidence_Angle_east_rad)/np.cos(SolarZeinthAngle_rad)
    if HourAngle_rad < - omega_rad or HourAngle_rad > omega_rad:
        r_d = 0
    else:
        r_d = np.pi/24 * (np.cos(HourAngle_rad) - np.cos(omega_rad)) / (np.sin(omega_rad) - omega_rad * np.cos(omega_rad))
    
    H_bar_d_0_east = H_bar_d_east(month, day, latitude, H_bar, TiltAngle_deg_east, n)
    I_bar_d_horizontal_east = r_d * H_bar_d_0_east
    I_bar_horizontal = Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime)
    I_bar_b_horizontal_east = I_bar_horizontal - I_bar_d_horizontal_east
    I_bar_b_east = I_bar_b_horizontal_east * r_b
    return I_bar_b_east

# Calculate beam solar radiation (west)
def Beam_solar_radiation_west(month, day, latitude, H_bar, SolarTime, SolarAltitudeAngle, SolarIncidenceAngle_west,TiltAngle_deg_west,n):
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = math.radians(omega_deg)
    HourAngle_rad = math.radians(Hour_angle(SolarTime))
    Solar_Incidence_Angle_west_rad = math.radians(SolarIncidenceAngle_west)
    SolarZeinthAngle_rad = math.radians(90 - SolarAltitudeAngle)
    
    r_b = math.cos(Solar_Incidence_Angle_west_rad)/math.cos(SolarZeinthAngle_rad)
    if HourAngle_rad < - omega_rad or HourAngle_rad > omega_rad:
        r_d = 0
    else:
        r_d = pi/24 * (math.cos(HourAngle_rad) - math.cos(omega_rad)) / (math.sin(omega_rad) - omega_rad * math.cos(omega_rad))
    
    H_bar_d_0_west = H_bar_d_west(month, day, latitude, H_bar, TiltAngle_deg_west, n)
    I_bar_d_horizontal_west = r_d * H_bar_d_0_west
    I_bar_horizontal = Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime)
    I_bar_b_horizontal_west = I_bar_horizontal - I_bar_d_horizontal_west
    I_bar_b_west = I_bar_b_horizontal_west * r_b
    return I_bar_b_west

#======================================
# Diffuse Solar Radiation Calculation 
#======================================

# Caculate diffuse solar radiation (south)
def Diffuse_solar_radiation_south(month, day, latitude, H_bar, SolarTime, TiltAngle_deg_south,n):
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = np.radians(omega_deg)
    HourAngle_rad = np.radians(Hour_angle(SolarTime))
    TiltAngle_south_rad = np.radians(TiltAngle_deg_south)
    if HourAngle_rad < - omega_rad or HourAngle_rad > omega_rad:
        r_d = 0
    else:
        np.pi / 24 * (np.cos(HourAngle_rad) - np.cos(omega_rad)) / (np.sin(omega_rad) - omega_rad * (np.cos(omega_rad)))
    H_bar_d_0 = H_bar_d_south(month, day, latitude, H_bar, TiltAngle_deg_south, n)
    I_bar_d_horizontal = r_d * H_bar_d_0
    I_bar_d_south = I_bar_d_horizontal * (1 + np.cos(TiltAngle_south_rad)) / 2
    return I_bar_d_south

# Caculate diffuse solar radiation (north)
def Diffuse_solar_radiation_north(month, day, latitude, H_bar, SolarTime, TiltAngle_deg_north,n):
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = np.radians(omega_deg)
    HourAngle_rad = np.radians(Hour_angle(SolarTime))
    TiltAngle_north_rad = np.radians(TiltAngle_deg_north)
    if HourAngle_rad < - omega_rad or HourAngle_rad > omega_rad:
        r_d = 0
    else:
        np.pi / 24 * (np.cos(HourAngle_rad) - np.cos(omega_rad)) / (np.sin(omega_rad) - omega_rad * (np.cos(omega_rad)))
    H_bar_d_0 = H_bar_d_north(month, day, latitude, H_bar, TiltAngle_deg_north, n)
    I_bar_d_horizontal = r_d * H_bar_d_0
    I_bar_d_north = I_bar_d_horizontal * (1 + np.cos(TiltAngle_north_rad)) / 2
    return I_bar_d_north

# Caculate diffuse solar radiation (east)
def Diffuse_solar_radiation_east(month, day, latitude, H_bar, SolarTime, TiltAngle_deg_east,n):
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = np.radians(omega_deg)
    HourAngle_rad = np.radians(Hour_angle(SolarTime))
    TiltAngle_east_rad = np.radians(TiltAngle_deg_east)
    if HourAngle_rad < - omega_rad or HourAngle_rad > omega_rad:
        r_d = 0
    else:
        np.pi / 24 * (np.cos(HourAngle_rad) - np.cos(omega_rad)) / (np.sin(omega_rad) - omega_rad * (np.cos(omega_rad)))
    H_bar_d_0 = H_bar_d_east(month, day, latitude, H_bar, TiltAngle_deg_east, n)
    I_bar_d_horizontal = r_d * H_bar_d_0
    I_bar_d_east = I_bar_d_horizontal * (1 + np.cos(TiltAngle_east_rad)) / 2
    return I_bar_d_east


# Caculate diffuse solar radiation (west)
def Diffuse_solar_radiation_west(month, day, latitude, H_bar, SolarTime, TiltAngle_deg_west,n):
    omega_deg = Sunset_Angle(month, day, latitude)
    omega_rad = np.radians(omega_deg)
    HourAngle_rad = np.radians(Hour_angle(SolarTime))
    TiltAngle_west_rad = np.radians(TiltAngle_deg_west)
    if HourAngle_rad < - omega_rad or HourAngle_rad > omega_rad:
        r_d = 0
    else:
        np.pi / 24 * (np.cos(HourAngle_rad) - np.cos(omega_rad)) / (np.sin(omega_rad) - omega_rad * (np.cos(omega_rad)))
    H_bar_d_0 = H_bar_d_west(month, day, latitude, H_bar, TiltAngle_deg_west, n)
    I_bar_d_horizontal = r_d * H_bar_d_0
    I_bar_d_west = I_bar_d_horizontal * (1 + np.cos(TiltAngle_west_rad)) / 2
    return I_bar_d_west

#===============================
# Ground Reflective Calculation 
#===============================

# Calculate ground reflective solar radiation (south)
def Ground_reflective_south(month, day, latitude, H_bar, TiltAngle_deg_south, SolarTime, T_air):
    #omega_deg = Sunset_Angle(month, day, latitude)
    #omega_rad = math.radians(omega_deg)
    #HourAngle = Hour_angle(SolarTime)
    TiltAngle_rad_south = np.radians(TiltAngle_deg_south)
    I_bar_horizontal = Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime)
    if T_air > 0 :
        rho_g = 0.2
    elif T_air < -5 :
        rho_g = 0.7
    else: 
        rho_g = -0.1 * T_air + 0.2
    I_bar_g_south = I_bar_horizontal * rho_g * ( 1 - np.cos(TiltAngle_rad_south)) / 2
    return I_bar_g_south


# Calculate ground reflective solar radiation (north)
def Ground_reflective_north(month, day, latitude, H_bar, TiltAngle_deg_north, SolarTime, T_air):
    #omega_deg = Sunset_Angle(month, day, latitude)
    #omega_rad = math.radians(omega_deg)
    #HourAngle = Hour_angle(SolarTime)
    TiltAngle_rad_north = np.radians(TiltAngle_deg_north)
    I_bar_horizontal = Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime)
    if T_air > 0 :
        rho_g = 0.2
    elif T_air < -5 :
        rho_g = 0.7
    else: 
        rho_g = -0.1 * T_air + 0.2
    I_bar_g_north = I_bar_horizontal * rho_g * ( 1 - np.cos(TiltAngle_rad_north)) / 2
    return I_bar_g_north

# Calculate ground reflective solar radiation (east)
def Ground_reflective_east(month, day, latitude, H_bar, TiltAngle_deg_east, SolarTime, T_air):
    #omega_deg = Sunset_Angle(month, day, latitude)
    #omega_rad = math.radians(omega_deg)
    #HourAngle = Hour_angle(SolarTime)
    TiltAngle_rad_east = np.radians(TiltAngle_deg_east)
    I_bar_horizontal = Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime)
    if T_air > 0 :
        rho_g = 0.2
    elif T_air < -5 :
        rho_g = 0.7
    else: 
        rho_g = -0.1 * T_air + 0.2
    I_bar_g_east = I_bar_horizontal * rho_g * ( 1 - np.cos(TiltAngle_rad_east)) / 2
    return I_bar_g_east

# Calculate ground reflective solar radiation (west)
def Ground_reflective_west(month, day, latitude, H_bar, TiltAngle_deg_west, SolarTime, T_air):
    #omega_deg = Sunset_Angle(month, day, latitude)
    #omega_rad = math.radians(omega_deg)
    #HourAngle = Hour_angle(SolarTime)
    TiltAngle_rad_west = np.radians(TiltAngle_deg_west)
    I_bar_horizontal = Solar_radiation_on_horizontal_GHI(month, day, latitude, H_bar, SolarTime)
    if T_air > 0 :
        rho_g = 0.2
    elif T_air < -5 :
        rho_g = 0.7
    else: 
        rho_g = -0.1 * T_air + 0.2
    I_bar_g_west = I_bar_horizontal * rho_g * ( 1 - np.cos(TiltAngle_rad_west)) / 2
    return I_bar_g_west


















       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
