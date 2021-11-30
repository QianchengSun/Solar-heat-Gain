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
df_house_15["declination_angle"] = np.round(sc.equation_of_time(month = df_house_15["Month"],
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
plt.plot(df_house_15["GHI"])
# caclulate solar altitude angle
df_house_15["solar_altitude_angle"] = sc.Solar_alitude_angle(latitude = latitude,
                                                            DeclinationAngle = df_house_15["declination_angle"],
                                                            SolarTime = df_house_15["solar_time"])
# %% calculate solar azimuth angle
df_house_15["solar_azimuth_angle"] = sc.Solar_azimuth_angle(latitude = latitude,
                                                            DeclinationAngle = df_house_15["declination_angle"],
                                                            SolarTime = df_house_15["solar_time"],
                                                            SolarAltitudeAngle = df_house_15["solar_altitude_angle"])

# %%
