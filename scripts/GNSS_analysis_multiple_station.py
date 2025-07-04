# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 09:08:46 2020

@author: zampa
"""

import inspect, os
from def_gnss_sql import gnss_sql
import numpy as np


# http://geodesy.unr.edu/NGLStationPages/gpsnetmap/GPSNetMap.html
# http://geodesy.unr.edu/NGLStationPages/steps.txt

# -----------------------------------------------------------------------------
# Here I create an external funtion, named  ''def_gps_sql.py'', imported as 'gps_sql'
# The function wants the following arguments:
# 1) The name of the input station file
# 2) A list with the bracking points dates e.g [12OCT25, 07JUN01, ...]
# 3) Name of the station e.g 'PORD'
# The outputs are: 
# 1) Figure with simple linear trend + bracking points, 
# 2) Figure with linear trend + annual component,
# 3) Dictionary where all data are stored.
# In the dictionary you have keys starting with 'l', refering to simple linear trend data
# and keys starting with 'a' refering to trend + annual component data
# r = residuals, i.e ar = (linear trend + annual component) residuals
# e = estimeted model values
# p = estimated parameters 
# e,n,u = east, north, up components
#Es: key 'ae_enu' = list with linear trend + annual component estimated values, for east north up componets

# System separator
s = os.sep 

# Home dir
home_dir = os.path.dirname( os.path.dirname( __file__ ) )


pord = gnss_sql( input_file_name = home_dir +s+ 'data' +s+ 'PORD.IGS08.tenv.txt', 
                 events = ['12OCT25'], 
                 output_fig = home_dir +s+ 'figs' +s+ 'PORD' )

codr = gnss_sql( input_file_name = home_dir +s+ 'data' +s+'CODR.IGS08.tenv.txt', 
                 events = [], 
                 output_fig =  home_dir +s+ 'figs' +s+ 'CODR')

# -----------------------------------------------------------------------------
