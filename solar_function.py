#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description: 
    The function edit is based on Dr. Chiasson's RCL 561 Solar Energy Engineering class in University of Dayton.
    The function is used to calculate the solar radiation on an object for every minutes.
How to use :
    Here is an example code for how to use all the functions as a package in python.
Example code : 
import sys
package_path = r"/Users/qianchengsun/PhD/github/Solar-heat-Gain"
sys.path.append(package_path)
import solar_calculation as sc
Created on Tue May 18, 2021
Recent update on Nov 29, 2021
@author: qianchengsun
"""
"""
All the function is applying on DataFrame calculation, 
The internal of calculation will be build based on Numpy.
"""
import numpy as np
# Day of the year
def day_of_year(month, day):
    # turn columns into numpy array
    month = np.array(month)
    day = np.array(day)
    # create list to store the value of day_of_year
    n_list = []
    for i in range(0, len(month)):
        if month[i] == 1:
            n = day[i]
        elif month[i] == 2 :
            n = 31 + day[i]
        elif month[i] == 3 :
            n =  59 + day[i]
        elif month[i] == 4 :
            n = 90 + day[i] 
        elif month[i] == 5 :
            n = 120 + day[i] 
        elif month[i] == 6 :
            n = 151 + day[i] 
        elif month[i] == 7 :
            n = 181 + day[i] 
        elif month[i] == 8 :
            n = 212 + day[i]
        elif month[i] == 9 :
            n = 243 + day[i]
        elif month[i] == 10 :
            n = 273 + day[i]
        elif month[i] == 11 :
            n = 304 + day[i]
        elif month[i] == 12 :
            n = 304 + day[i]
        else:
            n = 334 + day[i]
            
        n_list.append(n)
    # the value return data type is list
    return n_list