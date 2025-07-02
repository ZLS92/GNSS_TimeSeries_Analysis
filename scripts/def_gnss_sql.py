# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 09:04:26 2020

@author: zampa
"""

def gnss_sql(input_file_name, events, output_fig):
    
    import inspect, os
    folder_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    os.chdir(folder_path)
    import numpy as np
    from numpy.linalg import inv # to invert matrices 
    import pandas as pd
    import matplotlib.pyplot as plt
    
    #%% Load your GPS station file
    gps_file = pd.read_fwf(input_file_name,header=None)
    # -----------------------------------------------------------------------------
    
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
    
    #%%Linear Trend ---------------------------------------------------------------
    
    dat['s_enu'] = [dat['se'], dat['sn'], dat['su']] # list containing east north & up std_errors
    dat['enu'] = [dat['e'], dat['n'], dat['u']] # list containing east north & up vectors 
    dat['lp_enu'] = []
    dat['lr_enu'] = []
    dat['le_enu'] = []
    dat['lp_err_enu'] = []
    w = 2*np.pi # 1 year frequency for annual variations
    t = dat['year_dec'] # time vector
    numi = -1
    nume = len(events)
    
    # -----------------------------------------------------------------------------
    for i in dat['enu']:
        numi += 1
        Yo = i
        Yo = Yo.reshape((len(Yo),1)) # observed values shaped as column vector 
        A = np.ones([len(t), 2+nume]) # prepare designed matrix A 
        A[:,1] = t # velocity
        dat['event_index'] = []
        for e in range(0,nume):
            indx = dat['ymd'].index(events[e])
            A[0:indx,2+e] = 0
            dat['event_index'].append(indx)
        # Remember that the inverse of a digonal matrix is obtaine by replacing 
        # the digonal elements with the reciprocals 
        IQ = np.asmatrix(np.diag(1/dat['s_enu'][numi]**2)) # Inverse co-factor matrixr 
        cb = inv(A.T.dot(IQ).dot(A)).dot(A.T).dot(IQ).dot(Yo) # parameters vector        
        # Now we can append the results for each component (e,n,u) to lists in the dictionary
        dat['lp_enu'].append(cb) # list with parameters
        dat['le_enu'].append(np.asarray(A.dot(cb))) # list with parameters
        dat['lr_enu'].append(np.asarray(Yo-A.dot(cb))) # list with residuals
        dat['lp_err_enu'].append(np.sqrt(np.diag((((Yo-A.dot(cb)).T.dot(IQ).dot(Yo-A.dot(cb)))/
                                                   (len(Yo)-len(cb))).item()*(inv(A.T.dot(IQ).dot(A))))))           
    
    #%% Plot Figure Linear Trend ----------------------------------------------
    
    plt.figure(figsize=(22*0.393701,15*0.393701))
    nplot = -1
    y_labels = ['East (mm)', 'North (mm)', 'Up (mm)']
    te = dat['year_dec']
    mm = 1000 # factor m to mm conversion
    staring_year = int(np.min(dat['year_dec']))-1
    ending_year = int(np.max(dat['year_dec']))+2
    x_label = np.arange(staring_year, ending_year, step=1)
    
    for i in range(0, len(dat['enu'])):
        nplot += 2
        plt.subplot(3, 2, nplot)    
        ye = dat['le_enu'][i]
        plt.plot(te, ye*mm, 
                 label=f'f(t)= A + Bt')
        plt.scatter(te, (dat['enu'][i])*mm, 
                    s=1, c='r', label='Observed values')
        nume = 0
        for ev in dat['event_index']:
            nume += 1
            plt.axvline(x=te[ev], linestyle='dashed', c='c',label='discontinuity (Te'+str(nume)+')')     
        plt.xticks(x_label, rotation=45)
        if i==2:
            plt.xlabel('time')
        plt.ylabel(y_labels[i])
        if i==0:
            plt.legend(fontsize=8)
            plt.title('Observed data & Linear model')
        if i!=2:
            plt.tick_params(labelbottom=False)
        plt.grid('on', linestyle='--')
        
        plt.subplot(3, 2, nplot+1)
        plt.scatter(te, dat['lr_enu'][i]*mm, 
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
                        top=0.9, bottom=0.12, right=0.98, left=0.09)
    plt.suptitle( os.path.basename(output_fig))
    plt.savefig( output_fig +'_linear_trend.png' )    
    
    #%% Linear trend + annual component ---------------------------------------
    
    dat['s_enu'] = [dat['se'], dat['sn'], dat['su']] # list containing east north & up std_errors
    dat['ap_enu'] = []
    dat['ar_enu'] = []
    dat['ae_enu'] = []
    dat['ap_err_enu'] = []
    w = 2*np.pi # 1 year frequency for annual variations
    t = dat['year_dec'] # time vector
    numi = -1
    nume = len(events)
    
    # -------------------------------------------------------------------------
    for i in dat['enu']:
        numi += 1
        Yo = i-np.mean(i)
        Yo = Yo.reshape((len(Yo),1)) # observed values shaped as column vector 
        A = np.ones([len(t), 4+nume]) # prepare designed matrix A 
        A[:,1] = t # velocity
        A[:,2] = np.sin(w*t) # annual periodicity
        A[:,3] = np.cos(w*t) # annual periodicity
#        A[:,4] = np.sin(2*w*t) # semi-annual periodicity
#        A[:,5] = np.cos(2*w*t) # semi-annual periodicity
        dat['event_index'] = []
        for e in range(0,nume):
            indx = dat['ymd'].index(events[e])
            A[0:indx,4+e] = 0
            dat['event_index'].append(indx)
        # Remember that the inverse of a digonal matrix is obtaine by replacing 
        # the digonal elements with the reciprocals 
        IQ = np.asmatrix(np.diag(1/dat['s_enu'][numi]**2)) # Inverse co-factor matrixr 
        cb = inv(A.T.dot(IQ).dot(A)).dot(A.T).dot(IQ).dot(Yo) # parameters vector 
        # Now we can append the results for each component (e,n,u) to lists in the dictionary
        dat['ar_enu'].append(np.asarray(Yo-A.dot(cb))) # list with residuals 
        dat['ap_enu'].append(cb) # list with parameters
        dat['ae_enu'].append(np.asarray(A.dot(cb))) # list with parameters
        dat['ap_err_enu'].append(np.sqrt(np.diag((((Yo-A.dot(cb)).T.dot(IQ).dot(Yo-A.dot(cb)))/
                                                   (len(Yo)-len(cb))).item()*(inv(A.T.dot(IQ).dot(A))))))        
        
    #%% Plot trend + annual component ------------------------------------------
        
    plt.figure(figsize=(22*0.393701,15*0.393701))
    nplot = -1
    y_labels = ['East (mm)', 'North (mm)', 'Up (mm)']
    te = dat['year_dec']
    mm = 1000 # factor m to mm conversion
    staring_year = int(np.min(dat['year_dec']))-1
    ending_year = int(np.max(dat['year_dec']))+2
    x_label = np.arange(staring_year, ending_year, step=1)
    
    for i in range(0, len(dat['enu'])):
        nplot += 2
        plt.subplot(3, 2, nplot)    
        ye = dat['ae_enu'][i]
        plt.plot(te, ye*mm, 
                 label=f'f(t) = A + Bt + Csin(wt) + Dcos(wt)')
        plt.scatter(te, (dat['enu'][i]-np.mean(dat['enu'][i]))*mm, 
                    s=1, c='r', label='Observed values')  
        nume=0
        for ev in dat['event_index']:
            nume += 1
            plt.axvline(x=te[ev], linestyle='dashed', c='c',label='discontinuity (Te'+str(nume)+')')     
        plt.xticks(x_label, rotation=45)
        if i==2:
            plt.xlabel('time')
        plt.ylabel(y_labels[i])
        if i==0:
            plt.legend(fontsize=8)
            plt.title('Observed data & Linear model')
        if i!=2:
            plt.tick_params(labelbottom=False)
        plt.grid('on', linestyle='--')
        plt.text=()     
        
        
        plt.subplot(3, 2, nplot+1)
        plt.scatter(te, dat['ar_enu'][i]*mm, 
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
                        top=0.9, bottom=0.12, right=0.98, left=0.09)
    plt.suptitle( os.path.basename( output_fig ) )
    plt.savefig( output_fig+'_trend_annual.png' )   
    
    # -------------------------------------------------------------------------
    
    return dat