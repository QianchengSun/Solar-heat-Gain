"""
This is the python version code for 

Smart Wi-Fi physics-informed thermostat enabled to estimation of residential passive solar heat gain for any residence

"""

#%% programming about solar heat gain calculation based on LSTM-GRU

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
package_path = r"/Users/qianchengsun/PhD/github/Solar-heat-Gain"
sys.path.append(package_path)
import solar_calculation as sc

# set path to the prepared data frame
data_path = r"/Users/qianchengsun/PhD/solar_heat_gain/data/data_frame_house15.csv"
df_house_15 = pd.read_csv(data_path)
# apply solar function into df_house_15
# add new columns in exsit data frame
df_house_15["day_of_year"] = pd.NA
df_house_15["solar_time"] = pd.NA
df_house_15["GHI"] = pd.NA
"""
the order of the orientation is North, South, West, and East
"""
# create angles in different orientation
df_house_15["declination_angle"] = pd.NA
df_house_15["solar_altitude_angle"] = pd.NA
df_house_15["solar_azimuth_angle"] = pd.NA
# create solar incidence angle in different orientation
df_house_15["solar_incidence_angle_north"] = pd.NA
df_house_15["solar_incidence_angle_south"] = pd.NA
df_house_15["solar_incidence_angle_west"] = pd.NA
df_house_15["solar_incidence_angle_east"] = pd.NA
# create beam solar radiation on different orientation
df_house_15["beam_solar_radiation_north"] = pd.NA
df_house_15["beam_solar_radiation_south"] = pd.NA
df_house_15["beam_solar_radiation_west"] = pd.NA
df_house_15["beam_solar_radiation_east"] = pd.NA
# create diffuse solar radiation on different orientation
df_house_15["diffuse_solar_radiation_north"] = pd.NA
df_house_15["diffuse_solar_radiation_south"] = pd.NA
df_house_15["diffuse_solar_radiation_west"] = pd.NA
df_house_15["diffuse_solar_radiation_east"] = pd.NA
# create ground reflected solar radiation on different orientation
df_house_15["ground_reflected_solar_radiation_north"] = pd.NA
df_house_15["ground_reflected_solar_radiation_south"] = pd.NA
df_house_15["ground_reflected_solar_radiation_west"] = pd.NA
df_house_15["ground_reflected_solar_radiation_east"] = pd.NA
# surface azimuth angle define
df_house_15["surface_azimuth_angle_north"] = 180
df_house_15["surface_azimuth_angle_south"] = 0
df_house_15["surface_azimuth_angle_west"] = 90
df_house_15["surface_azimuth_angle_east"] = -90

# other input define
GMT = -5 # greenwhich time zone
Long_local = -84.2 # longitude for Dayton, OH
latitude = 39.7588 # latitude for Dayton, OH
"""
All the tilt angle define as 90 degree, is because in this case the tiltangle is represent for the residential wall
"""
# define tilt angle
TiltAngle_deg_north = 90
TiltAngle_deg_south = 90
TiltAngle_deg_west = 90
TiltAngle_deg_east = 90

# Surface Azimuth Angle define
SurfaceAzimuthAngle_south = 0
SurfaceAzimuthAngle_north = 180
SurfaceAzimuthAngle_west = 90
SurfaceAzimuthAngle_east = -90

# Turn the unit of temperature from F to C
T_air_outdoor = (df_house_15["Temp.Outdoor.Avg"] - 32) * 5 / 9
T_air_indoor = (df_house_15["Temp"] - 32) * 5 / 9
H_bar = df_house_15["All_sky_insolation_incident_on_horizontal_surface"]
# Apply the solar function to data frame to obtain the solar radiation
# combine hour and minutes together
df_house_15["new_hour"] = np.round(df_house_15["Hour"] + df_house_15["Minutes"] / 60, decimals = 3)
# calculate day_of_year
df_house_15["day_of_year"] = sc.day_of_year(month = df_house_15["Month"],
                                            day = df_house_15["Day"])
# declination angle
df_house_15["declination_angle"] = np.round(sc.Declination_Angle(month = df_house_15["Month"],
                                                                day = df_house_15["Day"]), decimals = 3)

# calculate solar time of the day
df_house_15["solar_time"] = np.round(sc.Solar_time(hour = df_house_15["new_hour"],
                                        GMT = GMT,
                                        long_local = Long_local,
                                        month = df_house_15["Month"],
                                        day = df_house_15["Day"]), decimals = 3)
# calculate the GHI (solar radiation on horizontal GHI)
df_house_15["GHI"] = np.round(sc.Solar_radiation_on_horizontal_GHI(month = df_house_15["Month"],
                                                        day = df_house_15["Day"],
                                                        latitude = latitude,
                                                        H_bar = H_bar, 
                                                        SolarTime = df_house_15["solar_time"]), decimals = 3)
# plt.plot(df_house_15["GHI"])
# caclulate solar altitude angle
df_house_15["solar_altitude_angle"] = sc.Solar_alitude_angle(latitude = latitude,
                                                            DeclinationAngle = df_house_15["declination_angle"],
                                                            SolarTime = df_house_15["solar_time"])

# calculate solar azimuth angle
"""
there is a warning in arccos(), that made lots of values are NAN in the dataframe
"""
df_house_15["solar_azimuth_angle"] = sc.Solar_azimuth_angle(latitude = latitude,
                                                            DeclinationAngle = df_house_15["declination_angle"],
                                                            SolarTime = df_house_15["solar_time"],
                                                            SolarAltitudeAngle = df_house_15["solar_altitude_angle"])

# calculate solar incidence angle 
"""
there is one point is different from Rstudio version in solar incidence angle calculation
""" 
# calculate solar incidence angle on North side
df_house_15["solar_incidence_angle_north"] = sc.Solar_incidence_angle_north(SurfaceAzimuthAngle_north = SurfaceAzimuthAngle_north,
                                                                            TiltAngle_deg_north = TiltAngle_deg_north,
                                                                            SolarAltitudeAngle_deg = df_house_15["solar_altitude_angle"],
                                                                            SolarAzimuthAngle = df_house_15["solar_azimuth_angle"])                                                                     
# calculate solar incidence angle on South side
df_house_15["solar_incidence_angle_south"] = sc.Solar_incidence_angle_south(SurfaceAzimuthAngle_south = SurfaceAzimuthAngle_south,
                                                                            TiltAngle_deg_south = TiltAngle_deg_south,
                                                                            SolarAltitudeAngle_deg = df_house_15["solar_altitude_angle"],
                                                                            SolarAzimuthAngle = df_house_15["solar_azimuth_angle"])                                                                       

# calculate solar incidence angle on West side
df_house_15["solar_incidence_angle_west"] = sc.Solar_incidence_angle_west(SurfaceAzimuthAngle_west = SurfaceAzimuthAngle_west,
                                                                            TiltAngle_deg_west = TiltAngle_deg_west,
                                                                            SolarAltitudeAngle_deg = df_house_15["solar_altitude_angle"],
                                                                            SolarAzimuthAngle = df_house_15["solar_azimuth_angle"])    
                                                                                                                                              
# calculate solar incidence angle on East side
df_house_15["solar_incidence_angle_east"] = sc.Solar_incidence_angle_east(SurfaceAzimuthAngle_east = SurfaceAzimuthAngle_east,
                                                                            TiltAngle_deg_east = TiltAngle_deg_east,
                                                                            SolarAltitudeAngle_deg = df_house_15["solar_altitude_angle"],
                                                                            SolarAzimuthAngle = df_house_15["solar_azimuth_angle"])    

# %%
"""
Beam solar radiation calculation has logic bug
"""
# calculate beam solar radiation
# calculate beam solar radiation on South side
df_house_15["beam_solar_radiation_south"] = sc.Beam_solar_radiation_south(month = df_house_15["Month"],
                                                                        day = df_house_15["Day"],
                                                                        latitude = latitude,
                                                                        H_bar = H_bar,
                                                                        SolarTime = df_house_15["solar_time"],
                                                                        SolarAltitudeAngle = df_house_15["solar_altitude_angle"],
                                                                        SolarIncidenceAngle_south = df_house_15["solar_incidence_angle_south"],
                                                                        TiltAngle_deg_south = TiltAngle_deg_south)
plt.plot(df_house_15["beam_solar_radiation_south"])
# %%
def H_bar_d_south(month, day, latitude, H_bar, TiltAngle_deg_south):
    if latitude < 0 :
        TiltAngle_deg_south = -TiltAngle_deg_south
    TiltAngle_south_rad = np.radians(TiltAngle_deg_south)

    H_bar_d_south_list = []
    Hbar_0 = np.array(sc.H_bar_0(month, day, latitude))
    omega_deg = np.array(sc.Sunset_Angle(month, day, latitude))

    for i in range(0, len(month)):
        k_t = H_bar[i] / Hbar_0[i]
        if omega_deg[i] < 81.4 :
            H_bar_d_south = np.max(0, ((1.391 - 3.56 * k_t + 4.189 * k_t ** 2 - 2.137 * k_t ** 3) * H_bar[i] * ( 1 + np.cos(TiltAngle_south_rad))/2) )
        else:
            H_bar_d_south = np.max(0, ((1.311 - 3.022 * k_t + 3.427 * k_t ** 2 - 1.821 * k_t ** 3) * H_bar[i] * ( 1 + np.cos(TiltAngle_south_rad))/2) )

        H_bar_d_south_list.append(H_bar_d_south)

    return H_bar_d_south_list
#%%
H_bar_d_south_value = H_bar_d_south(month = df_house_15["Month"],
                                    day = df_house_15["Day"],
                                    latitude = latitude,
                                    H_bar = H_bar,
                                    TiltAngle_deg_south = TiltAngle_deg_south)
# %%
sunset_angle = np.array(sc.Sunset_Angle(month = df_house_15["Month"], 
                                    day = df_house_15["Day"],
                                    latitude = latitude))
# %%
