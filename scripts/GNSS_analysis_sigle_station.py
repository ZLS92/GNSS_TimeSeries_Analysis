# -*- coding: utf-8 -*-
"""
Exercise: GPS Station time series modelling 
Course: Space Geodesy and InSAR
Professor: Alessandra Borghi
Created on Tue Feb 4 11:23:45 2020
@author: Luigi Sante Zampa
"""

"""
Exercise: GPS Station time series modelling 
Course: Space Geodesy and InSAR
Professor: Alessandra Borghi
Created on Tue Feb 4 11:23:45 2020
@author: Luigi Sante Zampa
"""

###############################################################################
# INPUT PARAMETERS 
# name of your input file
# http://geodesy.unr.edu/NGLStationPages/gpsnetmap/GPSNetMap.html
# http://geodesy.unr.edu/gps_timeseries/tenv/IGS08/
# The input file MUST be in the same folder of this script(namescript.py)!
# Otherwise, the input MUST be preceded by the absolute path 
# (e.g. C:\Users\Public\Desktop\folder_name\file.txt)  

input_file_name = 'PORD.IGS08.tenv.txt' # name of the input file (same folder of the script)

# Known discontinuities e.g. --> events=['12OCT25', '15DEC30', '18MAR05', ...]
# http://geodesy.unr.edu/NGLStationPages/steps.txt

events = ['12OCT25'] # in this data only one disc.
###############################################################################

# -----------------------------------------------------------------------------
# Import inspect and os libraries (we need them to navigate through your computer folders)

import inspect, os

# Then, find the exercise folder absolute path 

folder_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

# os.chdir command will set YOUR exercise folder as main searching path for this script 
# i.e. each file you want to load in the script MUST be in the same folder 
# where this namescript.py file is saved!!
# Otherwise,they MUST be preceded by the absolute path 
# (e.g. C:\Users\Public\Desktop\folder_name\file.txt)  

os.chdir(folder_path)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Import python library for linear algebra and general math operations
import numpy as np
from numpy.linalg import inv # function to invert matrices 
# import library to read/use "any sort of" table/database...
import pandas as pd
# import python library for plotting data and results
import matplotlib.pyplot as plt
# -----------------------------------------------------------------------------

#%% Load your GPS station file

# -----------------------------------------------------------------------------
# Read_fwf is a funtion from pandas library (pd) used to import txt files
# It creates an object called DataFrame from which you can select your columns
# E.g. gps_file['____up(m)'] --> gives the column with the vertical component of GPS

gps_file = pd.read_fwf(input_file_name,header=None)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Now, based on what we now about columns headers, we can create a dictionary. 
# A dictionary allows to store objects/variables (e.g. our columns), 
# and call them later on in the code, using our own keyword (es: {'key_word': column}) 
# I am doing this because, sometimes, it is simpler to work with dictionaries (native python), 
# rather than with pandas DataFrame objects and I can also choose keywords names which are easy to 
# remember.
# I named the dictionary dat, which stands for data... but u can call it whatever u like

dat={
     'year_dec': gps_file[2].values,
     'e': gps_file[6].values,
     'n': gps_file[7].values,
     'u': gps_file[8].values,
     'se': gps_file[10].values,
     'sn': gps_file[11].values,
     'su': gps_file[12].values,
     'ymd': gps_file[1].tolist(),
}
# -----------------------------------------------------------------------------

#%% Compute LSQ

# -----------------------------------------------------------------------------
# First, we need to define input and output vectors
dat['enu'] = [dat['e'], dat['n'], dat['u']] # list with east north & up vectors 
dat['s_enu'] = [dat['se'], dat['sn'], dat['su']] # list with east north & up std
dat['p_enu'] = [] # empty list to store estimated parameters
dat['e_enu'] = [] # empty list to store estimated values 
dat['r_enu'] = [] # empty list to store residuals
dat['p_err_enu'] = [] # empty list to store parameters std errors
w = 2*np.pi # annual frequency i.e period T=1 (using decimal year as time variable) 
t = dat['year_dec'] # time vector i.e. decimal years
nume = len(events) # number of discontinuities
numi = -1 # counter that we need for indexing (python indexes start from 0!!)

# -----------------------------------------------------------------------------
# Then, we compute and iterate the LSQ algorithm on each signal component: e, n, u
for i in dat['enu']:
    numi += 1	
    Yo = i-np.mean(i)# observed values
    Yo = Yo.reshape((len(Yo),1)) # observed values shaped as column vector 
    A = np.ones([len(t), 6+nume]) # prepare design matrix
    A[:,1] = t # velocity 
    A[:,2] = np.sin(w*t) # annual periodicity 1
    A[:,3] = np.cos(w*t) # annual periodicity 2
    A[:,4] = np.sin(2*w*t) # semi-annual periodicity 1
    A[:,5] = np.cos(2*w*t) # semi-annual periodicity 2
    dat['event_index'] = [] # empty list to store discontinuities in decimal year format
    for e in range(0,nume): # loop to set the step function for each discontinuity
        indx = dat['ymd'].index(events[e])
        A[0:indx,6+e] = 0 # all 0 before the discontinuity, after remains 1 
        dat['event_index'].append(indx)
    # The inverse of a diagonal matrix (i.e. Q) can be simply obtained by replacing 
    # the diagonal elements with the reciprocals 
    IQ = np.asmatrix(np.diag(1/dat['s_enu'][numi]**2)) # Inverse of co-factor matrix Q  
    #-------------------------------------------------------------------------
    cb = inv(A.T.dot(IQ).dot(A)).dot(A.T).dot(IQ).dot(Yo) # Vector with estimated parameters
    #-------------------------------------------------------------------------
    # Parameters std error
    cb_err= np.sqrt(np.diag((((Yo-A.dot(cb)).T.dot(IQ).dot(Yo-A.dot(cb)))/
                             (len(Yo)-len(cb))).item()*(inv(A.T.dot(IQ).dot(A)))))
    #-------------------------------------------------------------------------	 
    # Now we can append the results for each component (e,n,u) to a list in the dictionary
    dat['p_enu'].append(cb) # list with parameters
    dat['e_enu'].append(np.asarray(A.dot(cb))) # list with estimated values
    dat['r_enu'].append(np.asarray(Yo-A.dot(cb))) # list with residuals 
    dat['p_err_enu'].append(cb_err) # list with parameters errors
# -----------------------------------------------------------------------------

#%% Plot Results

# -----------------------------------------------------------------------------

plt.figure(figsize=(22*0.393701,18*0.393701))
nplot = -1
y_labels = ['East (mm)', 'North (mm)', 'Up (mm)']
te = dat['year_dec']
mm = 1000 # factor m to mm conversion
staring_year = np.int(np.min(dat['year_dec']))-1
ending_year = np.int(np.max(dat['year_dec']))+2
x_label = np.arange(staring_year, ending_year, step=1)

for i in range(0, len(dat['enu'])):
    nplot += 2
    plt.subplot(3, 2, nplot)    
    ye = dat['e_enu'][i]
    s = plt.plot(te, ye*mm, 
    label=f'f(t)= A+Bt\n +Csin(wt)+Dcos(wt)\n +Esin(2wt)+Fcos(2wt)\n +Gstep(t-Te)')
    l = plt.scatter(te, (dat['enu'][i]-np.mean(dat['enu'][i]))*mm, 
                    s=1, c='r', label='Observed values')
    nume =0
    for ev in dat['event_index']:
        nume += 1
        plt.axvline(x=te[ev], linestyle='dashed', c='c',label='discontinuity (Te_'+str(nume)+')')      
    plt.xticks(x_label, rotation=45)
    if i==2:
        plt.xlabel('time')
    plt.ylabel(y_labels[i])
    if i==0:
        plt.legend(fontsize=7)
        plt.title('Observed data & Linear model')
    if i!=2:
        plt.tick_params(labelbottom=False)
    plt.grid('on', linestyle='--')     
    
    plt.subplot(3, 2, nplot+1)
    l = plt.scatter(te, dat['r_enu'][i]*mm, 
                    s=1, c='k', label='Residuals')
    plt.xticks(x_label, rotation=45)
    if i==2:
        plt.xlabel('time')
    if i==0:
       plt.title('Residuals')
    if i!=2:
        plt.tick_params(labelbottom=False)
    plt.grid('on', linestyle='--')   
    
plt.subplots_adjust(wspace=0.11, hspace=0.08,
                    top=0.95, bottom=0.1, right=0.98, left=0.09)
plt.savefig('All components')

# -----------------------------------------------------------------------------

