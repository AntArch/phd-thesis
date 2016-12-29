# -*- coding: utf-8 -*-
"""
Created on Mon May 20 20:08:46 2013
Stuff to pull out stuff for full waveform

@author: davstott
davstott@gmail.com
"""

import os
import numpy as np
from scipy import interpolate
from scipy import signal
import matplotlib.pyplot as plt

import csv

#******************* PATHS****************************
dirpath = os.path.dirname(os.path.abspath('...'))
print ' running in: ', dirpath
datapath = os.path.join(dirpath,'data')
print 'data located in :', datapath

outputpath = os.path.join(dirpath,'output')
os.chdir(outputpath)

header = str('x'+','+
             'y'+','+
             'z'+','+
             'intensity'+','+
             'peak_start'+','+
             'peak_end'+','+
             'peak_location'+','+
             'peak_width'+','+
             'max_intensity'+','+
             'peak_sum'+','+
             'shoulder_location'+','+
             'shoulder_intensity'+
             '\n')



#*******************functions*****************************

#a smoothing spline
def smoothing(waveform, kparam, sparam, weights):
    sm_x = np.arange(1,257,1)
    sm_weights = np.zeros(256)+weights
    sm_spline = interpolate.UnivariateSpline(sm_x, 
                                             waveform, 
                                             k=kparam,
                                             w=sm_weights,
                                             s=sparam)

    spline = sm_spline(sm_x)
    return spline  
    

#***************** Parameters*********************************
#for spline
kparam = 1.3
sparam = 191
weights = 1.4

#x values
x_vals = np.arange(1,257,1)



#find the data
os.chdir(datapath)
indir = os.listdir(datapath) 
#open the data in a for loop
for file in indir:
    data_list = []
    reading_count = 0
    with open(file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if len(row)>8:
                data_list.append(row)
                #print len(row)
                reading_count = reading_count+1
        print reading_count, len(data_list)
        
    #read in the csv file as a numpy array
    #ata = np.genfromtxt(file, delimiter =',', skiprows=20)
    data = np.asarray(data_list, dtype=np.float32)
    print 'shape',data.shape
    
    #print 'dodgy row-1', data[63861,:]
    #print 'dodgy row', data[63862,:]
    #print 'dodgy row+1', data[63863,:]
    
    
    reading_count = data.shape[0]
    
    #create the output file & HEADER    
    os.chdir(outputpath)
    '''fname1 = str(file[0:-4])
    output_fname = fname1+'_output.txt'
    print 'output file:', output_fname
    output = open(output_fname,'w')
    output.write(header)
    output.close()'''
    
    j = 0
    
    fname1 = str(file[0:-4])
    output_fname = fname1+'_output.txt'
    
    array = []
    
    
    
    #for each row pull out the parameters by slicig the array
    for row in data:
        try:        
            j = j+1
            #print 'File:',file,'Line',j,'of', reading_count, 'processed'
            #print row.shape
            #osgb eastng
            x = row[0]
            #osgb northing
            y = row[1]
            #osgb /newlyn elevation
            z = row[2]
            #intensity derived from LAS tools
            intensity = row[3]
            #number of returns identified by lastools
            #returns = row[4]
            #and lastly pull out the waveform
            waveform = row[9:]
            #print waveform
                        
            #smooth the waveform using a univariate spline using the parameters above
            smoothed2 = smoothing(waveform, kparam, sparam, weights)
                        
            #identify the peaks in the smoothed waveform
            #peaks = signal.find_peaks_cwt(smoothed2, np.arange(19,27,1))
            
            #first derivative of smoothed waveform
            diff = np.diff(smoothed2, n=1)
            
            #second derivative of smoothed waveform
            #diff2 = np.diff(smoothed2, n=2)
            
            #find the maximal value in waveform
            max_intensity = np.argmax(waveform)
            #print 'MAX', max_intensity
            
            #define the region of the returns
            diffreg = np.logical_or(diff>1.5,diff<-0.75)
            #get the x values for slicing the other arrays
            diffx = x_vals[1:]
            regx  = diffx[diffreg]
            #get the first value
            reg_l = regx[0]
            #get the last value
            reg_r = regx[-1]
            #print 'diffreg', reg_l, reg_r
            
            shoulder = np.argmin(diff[reg_l:reg_r])
            #print 'shoulder pos', shoulder
                    
                    
            peak_value = waveform[max_intensity]
            #print peak_value
                    
            peak_width = reg_r-reg_l
            #print peak_width
            peak_sum = np.sum(waveform[reg_l:reg_r])
            shoulder_pos = shoulder+reg_l
            shoulder_int = waveform[shoulder_pos]
            #print shoulder_pos
            #print peak_sum
            
            os.chdir(outputpath)
            
            odata = np.column_stack((x,
                                     y,
                                     z,
                                     intensity,
                                     reg_l,
                                     reg_r,
                                     max_intensity,
                                     peak_width,
                                     peak_value,
                                     peak_sum,
                                     shoulder_pos,
                                     shoulder_int))
            #print 'STACKED', odata.shape               
            #array.append(odata[0,:])
            
            array.append(odata)
            
        except:
            continue
            
        
    print 'Length of array', len(array)
    processed = np.asarray(array)
    print processed
    print processed.shape
    np.savetxt(output_fname, processed[:,0,:], delimiter=',')

    os.chdir(datapath)