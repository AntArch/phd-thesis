# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 08:48:28 2014

@author: david
"""

#*************** IMPORT DEPENDANCIES*******************************************
import numpy as np
#import spec_gdal4 as spg
from osgeo import gdal
import os
import csv
import h5py
import datetime
import numpy.ma as ma
#from StringIO import StringIO
#import shapely
#import r2py

from osgeo import gdal_array
from osgeo import gdalconst
from osgeo.gdalconst import *

from scipy.spatial import ConvexHull
from scipy.signal import find_peaks_cwt
from scipy import interpolate

import matplotlib.pyplot as plt

from shapely.geometry import LineString

################# Functions ###################################################
'''These here are functions that are not part of any specific class- these 
are used by the data import classes for functions such as smoothing'''

def smoothing(perc_out, block_start, block_end, kparam, weight, sparam):   
    #D
    sm_spline_block = perc_out[block_start:block_end,:]
    sm_x = sm_spline_block[:,0]
    sm_y = sm_spline_block[:,1]
    sm_len = sm_x.shape
      
    sm_weights = np.zeros(sm_len)+weight
   
    
    sm_spline = interpolate.UnivariateSpline(sm_x, 
                                             sm_y, 
                                             k=kparam,
                                             w=sm_weights,
                                             s=sparam)
                                             
    
    spline = sm_spline(sm_x)
    
    
    spline = np.column_stack((sm_x,spline))
    
    
    
    return spline
    
    
def interpolate_gaps(array1, array2):          
   array_end = array1.shape[0]-1
   array1_endx = array1[array_end, 0]
   #get the start point of the second array
   array2_start = array2[0,0] 
   #get the length of the area to be interpolated
   x_len = array2_start-array1_endx+1
   #generate x values to use for the array
   xvals = np.linspace(array1_endx, array2_start, num=x_len)
   #y val for the start of  the interpolated area
   yval_array1 = array1[array_end,1]
   # y val for the end of interpolated area    
   yval_array2 = array2[0,1]
   #stack the values into a new array
   xin = np.append(array1_endx, array2_start)
   yin = np.append(yval_array1, yval_array2) 
   
   #numpy.interp(x, xp, fp)
   gap_filling = np.interp(xvals, xin, yin)
   filled_x = np.column_stack((xvals, gap_filling))
   print filled_x.shape
   return filled_x

class absorption_feature():
    '''this class is used for the characterisation of spectral absortion features, 
    and their investigation using continuum removal'''
    def __init__(self, spectra, feat_start, feat_end, feat_centre): 
        self.wl = spectra[:,0]
        self.values = spectra[:,1]
        print ('CALL TO ABSORPTION FEATURE')
        # start of absorption feature
        self.feat_start = feat_start
        # end of absorption feature
        self.feat_end = feat_end
        # approximate 'centre' of feature
        self.feat_centre = feat_centre        
        
        #get the range of the data
        self.min_wl = self.wl[0]
        self.max_wl = self.wl[-1]
        
        print ('Absorption feature',self.feat_start,self.feat_end)

        #define feature name
        self.feat_name = str(self.feat_start)+'_'+str(self.feat_end)
         
        '''# if the feature is within the range of the sensor, do stuff
        if self.feat_start > self.min_wl and self.feat_end < self.max_wl:
            print 'can do stuff with this data'
            try:
                self.abs_feature()
                print ('Absorption feature analysis sussceful')
            except:
                print ('ERROR analysing absorption feature', self.feat_name)
                pass
        else:
            print ('Cannot define feature: Out of range')'''                  
    
    ########## Methods ##################################################
    def abs_feature(self):
        print ('Call to abs_feature made')
        # Meffod to calculate the end points of the absorption feature
        # Does this using the Qhull algorithim form scipy spatial
        
        #use the initial defintnion of the absorption feature as a staring point
        # get the indices for these
        
        cont_rem_stacked = None
        ft_def_stacked = None
        
        
        
        start_point = np.argmin(np.abs(self.wl-self.feat_start))
        end_point = np.argmin(np.abs(self.wl-self.feat_end))
        centre = np.argmin(np.abs(self.wl-self.feat_centre))
            
        #find the index minima of reflectance 
        minima = np.argmin(self.values[start_point:end_point])+start_point
        
        # if the minima = the start point then the start point is the minima 
        if minima == start_point:
            left = minima
        #if not then the left side of the feature is the maixima on the left of the minima 
        elif minima <= centre:
            left = start_point+np.argmax(self.values[start_point:centre])    
        else:
            left = start_point+np.argmax(self.values[start_point:minima])
            
        #right is the maxima on the right of the absorption feature
        if minima == end_point:
            right = minima     
        else:
            right = minima+np.argmax(self.values[minima:end_point])
        
        # use left and right to create a 2D array of points
        hull_in = np.column_stack((self.wl[left:right],self.values[left:right]))
        
        #determine the minima of the points
        hull_min = minima-left
        
        if hull_min <= 0:
            hull_min=0
        
        #find the wavelength at minima
        hull_min_wl = hull_in[hull_min,0]
        
        # define the wavelength ranges we'll use to select simplices
        ft_left_wl = hull_min_wl-((hull_min_wl-hull_in[0,0])/2)
        ft_right_wl = hull_min_wl+((hull_in[-1,0]-hull_min_wl)/2)
        
        #use scipy.spatial convex hull to determine the convex hull of the points
        hull = ConvexHull(hull_in)
        
        # get the simplex tuples from the convex hull
        simplexes = hull.simplices
        # create an empty list to store simplices potentially related to our feature
        feat_pos = []
        #iterate through the simplices
        for simplex in simplexes:
            #extract vertices from simplices
            vertex1 = simplex[0]
            vertex2 = simplex[1]
            
            #print 'VERT!',hull_in[vertex1,0],hull_in[vertex2,0]
            
            ''' We're only interested in the upper hull. Qhull moves counter-
            clockwise. Therefore we're only interested in those points where 
            vertex 1 is greater than vertex 2'''
            '''The above may be total bollocks'''
            
            if not vertex1 < vertex2:
                '''We then use the wavelength ranges to determine which simplices
                relate to our absorption feature'''
                if hull_in[vertex2,0] <= ft_left_wl and \
                   hull_in[vertex2,0] >= self.wl[left] and \
                   hull_in[vertex1,0] >= ft_right_wl and \
                   hull_in[vertex1,0] <= self.wl[right]:
                    # append the vertices to the list
                    print hull_in[vertex2,0]
                    print hull_in[vertex1,0]
                    feat_pos.append((vertex2,vertex1))
                    print 'feat_pos length:',len(feat_pos), type(feat_pos)
                    #print feat_pos[0],feat_pos[1]
                
                else:
                    continue
        
        '''We only want one feature here. If there's more than one or less
        than one we're not interested as we're probably not dealing with
        vegetation'''
        # If there's less than one feature...
        if len(feat_pos) < 1:
            print ('Absorption feature cannot be defined:less than one feature')
            ft_def_stacked = None 
            ft_def_hdr = None 
            cont_rem_stacked = None

            
            
            
            
        elif len(feat_pos) == 1:
            feat_pos=feat_pos[0]
            print '£££££',feat_pos, type(feat_pos)
        
        else:
            #if theres more than one fid the widest one. this is not optimal.
            if len(feat_pos) >1:
                feat_width = []
                
                for pair in feat_pos:
                    feat_width.append(pair[1]-pair[0])
                    print 'feat width:', feat_width
                #feat_width = np.asarray(feat_width)
                print feat_width
                f_max = feat_width.index(max(feat_width))
                print f_max
                feat_pos = feat_pos[f_max]
                print type(feat_pos)
                print 'blaaarg',len(feat_pos)
                print feat_pos[0]
                print feat_pos[1]
            
       
        
        if not feat_pos==None:
            feat_pos = feat_pos[0], feat_pos[1]
            print 'DOES MY FEAT_POS CONVERSION WORK?', feat_pos
                           
            print ('Analysing absorption feature')
            #slice
            feature = hull_in[feat_pos[0]:feat_pos[1],:]
            print ('Feature shape',feature.shape,'start:',feature[0,0],'end:',feature[-1,0])
            #get the minima in the slice
            minima_pos = np.argmin(feature[:,1])
            #continuum removal
            contrem = self.continuum_removal(feature,minima_pos)
    
            # set up single value outputs
            # start of feature
            refined_start = feature[0,0]
            # end of feature
            refined_end = feature[-1,0]
            # wavelength at minima
            minima_WL = feature[minima_pos,0]
            # reflectance at minima
            minima_R = feature[minima_pos,1]
            # area of absorption feature
            feat_area = contrem[4]
            # two band normalised index of minima and start of feature
            left_tbvi = (refined_start-minima_R)/(refined_start+minima_R)
            # two band normalised index of minima and right of feature
            right_tbvi = (refined_end-minima_R)/(refined_end+minima_R)
            # gradient of the continuum line
            cont_gradient = np.mean(np.gradient(contrem[0]))
            # area of continuum removed absorption feature
            cont_rem_area = contrem[3]
            # maxima of continuum removed absorption feature
            cont_rem_maxima = np.max(contrem[1])
            # wavelength of maxima of continuum removed absorption feature
            cont_rem_maxima_wl = feature[np.argmax(contrem[1]),0]
            #area of left part of continuum removed feature
            cont_area_l = contrem[5]
            if cont_area_l == None:
                cont_area_l=0
            #are aof right part of continuum removed feature
            cont_area_r = contrem[6]

            #stack these into a lovely array
            ft_def_stacked = np.column_stack((refined_start,
                                              refined_end,
                                              minima_WL,
                                              minima_R,
                                              feat_area,
                                              left_tbvi,
                                              right_tbvi,
                                              cont_gradient,
                                              cont_rem_area,
                                              cont_rem_maxima,
                                              cont_rem_maxima_wl,
                                              cont_area_l,
                                              cont_area_r))
                                              
            ft_def_hdr = str('"Refined start",'+
                               '"Refined end",'+
                               '"Minima Wavelenght",'+
                               '"Minima Reflectance",'+
                               '"Feature Area",'+
                               '"Left TBVI",'+
                               '"Right TBVI",'+
                               '"Continuum Gradient",'+
                               '"Continuum Removed Area",'+
                               '"Continuum Removed Maxima",'+
                               '"Continuum Removed Maxima WL",'+
                               '"Continuum Removed Area Left",'+
                               '"Continuum Removed Area Right",')
                                              
            print ft_def_stacked.shape #save the stacked outputs as hdf
                
                
                                                  
            # stack the 2d continuum removed outputs
            cont_rem_stacked = np.column_stack((feature[:,0],
                                                feature[:,1],
                                                contrem[0],
                                                contrem[1],
                                                contrem[2]))
                                                
            print ('CREM', cont_rem_stacked.shape)
                                            
        return ft_def_stacked, ft_def_hdr, cont_rem_stacked
            
            
            
    def continuum_removal(self,feature,minima):  
        #method to perform continuum r=<emoval
        #pull out endmenmbers
        end_memb = np.vstack((feature[0,:],feature[-1,:]))
        #interpolate between the endmembers using x intervals    
        continuum_line = np.interp(feature[:,0], end_memb[:,0], end_memb[:,1])
        #continuum removal
        continuum_removed = continuum_line/feature[:,1]

        #stack into coord pairs so we can measure the area of the feature
        ft_coords = np.vstack((feature,
                               np.column_stack((feature[:,0],continuum_line))))
        #get the area
        area = self.area(ft_coords)
                
        #get the area of the continuum removed feature
        cont_rem_2d = np.column_stack((feature[:,0],continuum_removed))        
        cont_r_area = self.area(cont_rem_2d)

        #band-normalised by area continuum removal
        cont_BNA = (1-(feature[:,1]/continuum_line))/area
        
        #continuum removed area on left of minima
        cont_area_left = self.area(cont_rem_2d[0:minima,:])

        #continuum removed area on right of minima
        cont_area_right = self.area(cont_rem_2d[minima:,:])
        
        return (continuum_line, 
                continuum_removed, 
                cont_BNA, 
                cont_r_area,  
                area,
                cont_area_left,
                cont_area_right)
         
    #define area of 2d polygon- using shoelace formula
    def area(self, coords2d):
        #setup counter
        total = 0.0
        #get the number of coorsinate pairs
        N = coords2d.shape[0]
        #iterate through these
        for i in range(N):
            #define the first coordinate pair
            vertex1 = coords2d[i]
            #do the second
            vertex2 = coords2d[(i+1) % N]
            #append the first & second distance to the toatal
            total += vertex1[0]*vertex2[1] - vertex1[1]*vertex2[0]
            #return area
            return abs(total/2)  
    
class Indices():
    #class that does vegetation indices
    def __init__(self,spectra):
        self.wl = spectra[:,0]
        self.values = spectra[:,1]
        self.range = (np.min(self.wl),np.max(self.wl))
        '''So, the init method here checks the range of the sensor and runs
        the appropriate indices within that range, and saves them as hdf5. 
        The indices are all defined as methods of this class'''
    def visnir(self):
    # Sensor range VIS-NIR
        if self.range[0] >= 350 and \
           self.range[0] <= 500 and \
           self.range[1] >= 900:
               vis_nir = np.column_stack((self.sr700_800(),
                                          self.ndvi694_760(),
                                          self.ndvi695_805(),
                                          self.ndvi700_800(),
                                          self.ndvi705_750(),
                                          self.rdvi(),
                                          self.savi(),
                                          self.msavi2(),
                                          self.msr(),
                                          self.msrvi(),
                                          self.mdvi(),
                                          self.tvi(),
                                          self.mtvi(),
                                          self.mtvi2(),
                                          self.vog1vi(),
                                          self.vog2(),
                                          self.prsi(),
                                          self.privi(),
                                          self.sipi(),
                                          self.mcari(),
                                          self.mcari1(),
                                          self.mcari2(),
                                          self.npci(),
                                          self.npqi(),
                                          self.cri1(),
                                          self.cri2(),
                                          self.ari1(),
                                          self.ari2(),
                                          self.wbi()))
                                          
                                          
               vis_nir_hdr=str('"sr700_800",'+
                               '"ndvi694_760",'+
                               '"ndvi695_805",'+
                               '"ndvi700_800",'+
                               '"ndvi705_750",'+
                               '"rdvi",'+
                               '"savi",'+
                               '"msavi2",'+
                               '"msr",'+
                               '"msrvi",'+
                               '"mdvi",'+
                               '"tvi",'+
                               '"mtvi",'+
                               '"mtvi2",'+
                               '"vog1vi",'+
                               '"vog2",'+
                               '"prsi"'+
                               '"privi",'+
                               '"sipi",'+
                               '"mcari",'+
                               '"mcari1",'+
                               '"mcari2",'+
                               '"npci",'+
                               '"npqi",'+
                               '"cri1",'+
                               '"cri2",'+
                               '"ari1",'+
                               '"ari2",'+
                               '"wbi"')
        else:
            vis_nir = None
            vis_nir_hdr = None
        
        return vis_nir,vis_nir_hdr
            #Range NIR-SWIR
        
    def nir_swir(self):
        if self.range[0] <= 900 and self.range[1] >=2000:
            nir_swir = np.column_stack((self.ndwi(),
                                        self.msi(),
                                        self.ndii()))
            nir_swir_hdr = str('"ndwi",'+
                               '"msi",'+
                               '"ndii"')                                
            
        else:
            #continue
            print 'not nir-swir'
            nir_swir=None
            nir_swir_hdr=None            
            
            
        return nir_swir, nir_swir_hdr
        #range SWIR
    def swir(self):
        if self.range[1] >=2000:
            swir = np.column_stack((self.ndni(),
                                    self.ndli()))
            
            swir_hdr=str('"ndni",'+
                         '"ndli"')
        else:
            print 'swir-nir'
            swir = None
            swir_hdr = None
            #continue
        return swir,swir_hdr
        
        
    #||||||||||||||||||||| Methods |||||||||||||||||||||||||||||||||||||||||||||||        
    # function to run every permutation of the NDVI type index across the Red / IR
    # ...... VIS / NIR methods ....
    def multi_tbvi (self, red_start=650, red_end=750, ir_start=700, ir_end=850):
        # get the indicies of the regions we're going to use.
        # we've added default values here, but they can happily be overidden
        #start of red
        red_l =np.argmin(np.abs(self.wl-red_start))
        #end of red
        red_r = np.argmin(np.abs(self.wl-red_end))
        #start of ir
        ir_l = np.argmin(np.abs(self.wl-ir_start))
        #end of ir
        ir_r = np.argmin(np.abs(self.wl-ir_end))
        
        #slice
        left = self.values[red_l:red_r]
        right = self.values[ir_l:ir_r]
        
        #set up output
        values = np.empty(3)
        
        #set up counter
        l = 0    
        #loop throught the values in the red
        for lvalue in left:
            l_wl = self.wl[l+red_l]
            r = 0
            l = l+1
            #then calculate the index with each wl in the NIR
            for rvalue in right:
                value = (rvalue-lvalue)/(rvalue+lvalue)
                r_wl = self.wl[r+ir_l]
                out = np.column_stack((l_wl,r_wl,value))
                values = np.vstack((values, out))
                out = None
                r = r+1
        return values[1:,:]
        
    def sr700_800 (self, x=700, y=800):
        index = self.values[np.argmin(np.abs(self.wl-x))]/self.values[np.argmin(np.abs(self.wl-y))]
        return index
    
    def ndvi705_750 (self, x=705, y=750):
        index = (self.values[np.argmin(np.abs(self.wl-y))]-self.values[np.argmin(np.abs(self.wl-x))])/\
                (self.values[np.argmin(np.abs(self.wl-y))]+self.values[np.argmin(np.abs(self.wl-x))])
        return index    
    
    def ndvi700_800 (self, x=700, y=800):
        index = (self.values[np.argmin(np.abs(self.wl-y))]-self.values[np.argmin(np.abs(self.wl-x))])/\
                (self.values[np.argmin(np.abs(self.wl-y))]+self.values[np.argmin(np.abs(self.wl-x))])
        return index
        
    def ndvi694_760 (self, x=694, y=760):
        index = (self.values[np.argmin(np.abs(self.wl-y))]-self.values[np.argmin(np.abs(self.wl-x))])/\
                (self.values[np.argmin(np.abs(self.wl-y))]+self.values[np.argmin(np.abs(self.wl-x))])
        return index
    
    def ndvi695_805 (self, x=695, y=805):
        index = (self.values[np.argmin(np.abs(self.wl-y))]-self.values[np.argmin(np.abs(self.wl-x))])/\
                (self.values[np.argmin(np.abs(self.wl-y))]+self.values[np.argmin(np.abs(self.wl-x))])
        return index
        
    def npci (self, x=430, y=680):
        index = (self.values[np.argmin(np.abs(self.wl-y))]-self.values[np.argmin(np.abs(self.wl-x))])/\
                (self.values[np.argmin(np.abs(self.wl-y))]+self.values[np.argmin(np.abs(self.wl-x))])
        return index
        
    def npqi (self, x=415, y=435):
        index = (self.values[np.argmin(np.abs(self.wl-y))]-self.values[np.argmin(np.abs(self.wl-x))])/\
                (self.values[np.argmin(np.abs(self.wl-y))]+self.values[np.argmin(np.abs(self.wl-x))])
        return index
    
        #mSRvi
    #= (750-445)/(705+445)
    def msrvi (self):   
        x = 750
        y = 445
        z = 705    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        msrvi_val = (x_val-y_val)/(z_val+y_val)
        return msrvi_val
        
    #Vogelmann Red Edge 1
    #740/720
    def vog1vi (self):   
        x = 740
        y = 720    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        vog1vi_val = (x_val/y_val)
        return vog1vi_val
        
    #Vogelmann Red Edge 2
    #= (734-747)/(715+726)
    def vog2 (self):   
        v = 734
        x = 747
        y = 715
        z = 726    
        v_val = self.values[np.argmin(np.abs(self.wl-v))]
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        vog2_val = (v_val-x_val)/(y_val+z_val)
        return vog2_val
        
    #PRI
    # (531-570)/(531+570)
    def privi (self):   
        x = 531
        y = 570    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        privi_val = (x_val-y_val)/(x_val+y_val)
        return privi_val
                
    #SIPI
    #(800-445)/(800-680)
    def sipi (self):   
        x = 800
        y = 445
        z = 680    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        sipi_val = (x_val-y_val)/(x_val+z_val)
        return sipi_val

    #Water band index
    # WBI = 900/700
    def wbi (self):   
        x = 900
        y = 700    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        wbi_val = (x_val/y_val)
        return wbi_val
        
    #mNDVI
    #= (750-705)/((750+705)-(445))
    def mdvi (self):
        x = 750
        y = 705
        z = 445    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        mdvi_val = (x_val-y_val)/((x_val+y_val)-z_val)
        return mdvi_val
        
    #Carotenid Reflectance Index
    #CRI1 = (1/510)-(1/550)
    def cri1 (self):   
        x = 510
        y = 550    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        cri1_val = (1/x_val)-(1/y_val)
        return cri1_val
        
    #CRI2 = (1/510)-(1/700)
    def cri2 (self):   
        x = 510
        y = 700    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        cri2_val = (1/x_val)-(1/y_val)
        return cri2_val
        
    #Anthocyanin
    #ARI1 = (1/550)-(1/700)
    def ari1 (self):   
        x = 550
        y = 700   
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        ari1_val = (1/x_val)-(1/y_val)
        return ari1_val
        
    #ARI2 = 800*((1/550)-(1/700)_))
    def ari2 (self):
        x = 510
        y = 700   
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        ari2_val = 800*((1/x_val)-(1/y_val))
        return ari2_val
        
    #MSR
    #=((800/670)-1)/SQRT(800+670)
    def msr (self):   
        x = 800
        y = 670   
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        msr_val = ((x_val/y_val)-1)/(np.sqrt(x_val+y_val))
        return msr_val
        
    #SAVI
    #= (1+l)(800-670)/(800+670+l)
    def savi (self, l=0.5):   
        x = 800
        y = 670
        l = 0.5
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        savi_val = ((1+l)*(x_val-y_val))/(x_val+y_val+l)
        return savi_val
        
    #MSAVI
    #=1/2(sqrt(2*800)+1)-SQRT(((2*800+1)sqr)-8*(800-670)
    def msavi2 (self):   
        x = 800
        y = 670
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        msavi2_top1 = (2*x_val+1)
        msavi2_top2 = (np.sqrt(np.square(2*x_val+1)-(8*(x_val-y_val))))
        msavi2_top = msavi2_top1-msavi2_top2
        msavi2_val = msavi2_top/2
        return msavi2_val
        
    #Modified clhoropyll absorption indec
    #MCARI = ((700-670)-0.2*(700-550))*(700/670)
    def mcari (self):   
        x = 700
        y = 670
        z = 550
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        mcari_val = (x_val-y_val)-(0.2*(x_val-z_val)*(x_val/y_val))    
        return mcari_val
        
    #Triangular vegetation index
    #TVI 0.5*(120*(750-550))-(200*(670-550))
    def tvi (self):   
        x = 750
        y = 550
        z = 670    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        tvi_val = 0.5*((120*(x_val-y_val))-(200*(z_val+y_val)))
        return tvi_val
        
    #MCAsavRI1 = 1.2*(2.5*(800-67-)-(1.3*800-550)
    def mcari1 (self):   
        x = 800
        y = 670
        z = 550    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        mcari1_val = (1.2*((2.5*(x_val-y_val)))-(1.3*(x_val+z_val)))
        return mcari1_val
        
    #MTVI1
    #=1.2*((1.2*(800-550))-(2.5(670-550)))
    def mtvi (self):   
        x = 800
        y = 550
        z = 670    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        mtvi_val = (1.2*(12*(x_val-y_val)))-(2.5*(z_val-y_val))
        return mtvi_val

    def mcari2 (self):   
        x = 800
        y = 670
        z = 550
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        mcari2_top = (1.5*(2.5*(x_val-y_val)))-(1.3*(x_val-z_val))
        mcari2_btm = np.sqrt((np.square(2*x_val)+1)-((6*x_val)-(5*(np.sqrt(y_val))))-0.5)
        mcari2_val = mcari2_top/mcari2_btm
        return mcari2_val
        
    #MTVI2=(1.5*(2.5(800-670)-2.5*(800-550))/sqrt((2*800+1s)sq)-((6*800)-(5*sqrt670))-0.5
    def mtvi2 (self):   
        x = 800
        y = 670
        z = 550    
    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        mtvi2_top = (1.5*(2.5*(x_val-z_val)))-(1.3*(x_val-z_val))
        mtvi2_btm = np.sqrt((np.square(2*x_val)+1)-((6*x_val)-(5*(np.sqrt(y_val))))-0.5)
        mtvi2_val = mtvi2_top/mtvi2_btm
        return mtvi2_val
        
    #Renormalised DVI
    #RDVI = (800-670)/sqrt(800+670)
    def rdvi (self):
        x = 800
        y = 670
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        rdvi_val = (x_val-y_val)/np.sqrt(x_val+y_val)
        return rdvi_val
        
    #Plant senescance reflectance index
    #PRSI = (680-500)/750
    def prsi (self):   
        x = 680
        y = 500
        z = 750
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        prsi_val = (x_val-y_val)/z_val
        return prsi_val
        
    #||||||||||||||||||||||| SWIR methods ||||||||||||||||||||||||||||||||||||

    #Cellulose Absorption Index
    #CAI =0.5*(2000-2200)/2100
    def cai (self):   
        x = 2000
        y = 2200
        z = 2100
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        z_val = self.values[np.argmin(np.abs(self.wl-z))]
        cai_val = 0.5*(x_val-y_val)-z_val
        return cai_val
        
    #Normalized Lignin Difference
    #NDLI = (log(1/1754)-log(1/1680))/(log(1/1754)+log(1/1680))
    def ndli (self):   
        x = 1754
        y = 2680
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))] 
        ndli_val = (np.log(1/x_val)-np.log(1/y_val))/(np.log(1/x_val)+np.log(1/y_val))
        return ndli_val
        
    #Canopy N
    #NDNI =(log(1/1510)-log(1/1680))/(log(1/1510)+log(1/1680))
    def ndni (self):   
        x = 1510
        y = 1680
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        ndni_val = (np.log(1/x_val)-np.log(1/y_val))/(np.log(1/x_val)+np.log(1/y_val))
        return ndni_val
        
    #|||||||||||||||||||||| Full spectrum (VIS-SWIR)||||||||||||||||||||||||||||
    #Normalised Difference IR index
    #NDII = (819-1649)/(819+1649)#NDII = (819-1649)/(819+1649)
    def ndii (self):   
        x = 819
        y = 1649    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        ndii_val = (x_val-y_val)/(x_val+y_val)
        return ndii_val
    
    #Moisture Stress Index
    #MSI = 1599/819http://askubuntu.com/questions/89826/what-is-tumblerd
    def msi (self):   
        x = 1599
        y = 810    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        msi_val = (x_val/y_val)
        return msi_val
        
    #NDWI
    #(857-1241)/(857+1241)
    def ndwi (self):   
        x = 857
        y = 1241    
        x_val = self.values[np.argmin(np.abs(self.wl-x))]
        y_val = self.values[np.argmin(np.abs(self.wl-y))]
        ndwi_val = (x_val-y_val)/(x_val+y_val)
        return ndwi_val

class red_edge():
    '''Class to derive red edge position using a number of different methods'''
    def __init__(self, spectra):
        self.wl = spectra[:,0]
        self.values = spectra[:,1]
        self.range = (np.min(self.wl),np.max(self.wl))
        '''Again, the mehtod that initialises this class uses the range of the 
        sensor to check to see if it falls within the red-edge reigion. If so,
        it will derive the red edge using the differnet methods and save these 
        as seprate hdf5 datasets in the appropriate group''' 
        if self.range[0] <= 670 and self.range[1] >=750:
            self.redge_vals = np.column_stack((self.redge_linear(),
                                               self.redge_lagrange(),
                                               self.redge_linear_extrapolation()))
            print self.redge_vals
            print self.redge_linear,self.redge_lagrange,self.redge_linear_extrapolation
            
            self.redge_hdr = str('"linear",'+
                                 '"lagrange",'+
                                 '"extrapolated"')
            
        else:
            print ('red_edge out of range')
            self.redge_vals = None
            self.redge_hdr = None
            
            
    ##################### METHODS #########################################
    #linear- defined by clevers et al 1994:
    def redge_linear(self):
        r670 = self.values[np.argmin(np.abs(self.wl-670))]
        r780 = self.values[np.argmin(np.abs(self.wl-780))]
        r700 = self.values[np.argmin(np.abs(self.wl-700))]
        r740 = self.values[np.argmin(np.abs(self.wl-740))]
        r_edge = (r670+r780)/2
        lin_rep =700+40*((r_edge-r700)/(r740-r700))
        print 'REDGE_LINEAR',lin_rep
        return lin_rep
        
    #Lagrangian method, after Dawson & Curran 1998
    def redge_lagrange(self):        
        #select the red edge region of the first derviative and associate this
        #with wavelength
        x = 680
        y = 730   
        first_diff = np.diff(self.values, 1) 
        spec_in = np.column_stack((self.wl[1:], first_diff))
        l680 = np.argmin(np.abs(spec_in[:,0]-x))
        r680 = spec_in[l680,0]
        l730 = np.argmin(np.abs(spec_in[:,0]-y))
        r730 = spec_in[l730,0]
        redge_region_sel = np.where(np.logical_and(spec_in[:,0]>r680-1, 
                                                   spec_in[:,0]<r730+1))
        redge_region = spec_in[redge_region_sel]        
        #find the maximum first derivative, return index
        dif_max = np.argmax(redge_region[:,1], axis=0)
        #find band with the max derivative -1, return index     
        dif_max_less = (np.argmax(redge_region[:,1], axis=0))-1
        #find band with the max derivative +1, return index
        dif_max_more = (np.argmax(redge_region[:,1], axis=0))+1 
        
        
        if dif_max_more >= redge_region.shape[0]:
            dif_max_more = redge_region.shape[0]-1
        
        #use these indeces to slice the array    
        rmax = redge_region[dif_max]    
        rmax_less =redge_region[dif_max_less]
        rmax_more =redge_region[dif_max_more]
        #lagrangian interpolation with three points
        #this has been expanded to make the syntax easier    
        a = rmax_less[1]/(rmax_less[0]-rmax[0])*(rmax_less[0]-rmax_more[0])
        b = rmax[1]/(rmax[0]-rmax_less[0])*(rmax[0]-rmax_more[0])
        c = rmax_more[1]/(rmax_more[0]-rmax_less[0])*(rmax_more[0]-rmax[0])
        d = a*(rmax[0]+rmax_more[0])
        e = b*(rmax_less[0]+rmax_more[0])
        f = c*(rmax_less[0]+rmax[0])    
        lg_rep = (d+e+f)/(2*(a+b+c))
        print 'Lagrangian', lg_rep
        return lg_rep
    
    
    
    
    
    #Linear extrapolation- after Cho & Skidmore 2006, Cho et al 2007
    def redge_linear_extrapolation(self):
        diff = np.diff(self.values)
        d680 = diff[np.argmin(np.abs(self.wl-680+1))]
        d694 = diff[np.argmin(np.abs(self.wl-694+1))]
        
        d724 = diff[np.argmin(np.abs(self.wl-724+1))]
        d760 = diff[np.argmin(np.abs(self.wl-760+1))]
        
        red_slope = ((d694-d680)/(694-680))
        ir_slope = ((d760-d724)/(760-724))
        
        red_inter = d680-(red_slope*680)
        
        ir_inter = d724-(ir_slope*724)
        
        wl = (ir_inter-red_inter)/(ir_slope-red_slope)

        print '^!!!!!!!!! Linear:',wl
        return np.abs(wl)

        

        
        

class fluorescence():
    '''this class is inteded to look for evidence of photosynthetic flourescence
    currently this is limited to simple reflectance indices. This should be
    expanded to take in other more complex methods to invesitgae fluorescence'''
    
    def __init__(self, spectra):
        self.wl = spectra[:,0]
        self.values = spectra[:,1]
        self.range = (np.min(self.wl),np.max(self.wl))
        print 'call to fluor'
        '''The init method checks the range to establish if it overlaps with 
        region of chlorophyll flourescence. If so it will will perform the 
        analysis methods and output to hdf5'''
        

    def wl_selector(self, x):
        '''this method finds the index of the wavelength closest to that 
        specified for reflectance'''
        value = self.values[np.argmin(np.abs(self.wl-x))]
        return value
    
    def d_wl_selector(self, x):
        '''this method finds the index of the wavelength closest to that 
        specified for the first derivative'''
        diff = np.diff(self.values)
        value = diff[np.argmin(np.abs(self.wl-x))+1]
        return value
        
    def wl_max_d(self):
        '''method to extract wavelength of the maxima of the first derivative
        and return this'''
        start = np.argmin(np.abs(self.wl-650))
        end = np.argmin(np.abs(self.wl-760))
        diff = np.diff(self.values[start:end])
        maxdiff = np.argmax(diff)
        maxdiffwl = self.wl[maxdiff+start+1]
        return maxdiffwl, diff[maxdiff]
        
    def simple_ratios(self):
        ''' This method runs flourescence indices ratios and returns them as a
        stacked numpy array'''
        #r680/r630
        r680r630 = self.wl_selector(680)/self.wl_selector(630)
        print r680r630
        #r685/r630
        r685r630 = self.wl_selector(685)/self.wl_selector(630)
        print r685r630
        #r685/r655
        r685r655 = self.wl_selector(685)/self.wl_selector(655)
        print r685r655
        #r687/r630
        r687r630 = self.wl_selector(687)/self.wl_selector(630)
        print r687r630
        #r690/r630
        r690r630 = self.wl_selector(690)/self.wl_selector(630)
        print r690r630
        #r750/r800
        r750r800 = self.wl_selector(750)/self.wl_selector(800)
        print r750r800
        
        #sq(r685)/(r675-r690)
        sqr685 = np.square(self.wl_selector(685))/(self.wl_selector(675)-self.wl_selector(690))
        print sqr685
        
        #(r675-r690)/sq(r683) Zarco-Tejada 2000
        r675r690divsq683 = (self.wl_selector(675)-self.wl_selector(690))/np.square(self.wl_selector(683))
        print r675r690divsq683
         
        #d705/d722
        d705d722 = self.d_wl_selector(705)/self.d_wl_selector(722)
        print d705d722
        
        #d730/d706
        d730d706 = self.d_wl_selector(730)/self.d_wl_selector(706)
        print d730d706
        
        #(d688-d710)/sq(d697)
        d686d710sq697 = (self.d_wl_selector(688)-self.d_wl_selector(710))\
                        /np.square(self.d_wl_selector(697))
        print d686d710sq697
        
        #wl at max d / d720
        maxdd720 = self.wl_max_d()[1]/self.d_wl_selector(720)
        print maxdd720
        
        #wl at max d / d703
        maxdd703 = self.wl_max_d()[1]/self.d_wl_selector(703)
        print maxdd703
        
        #wl at max d / d(max d+12)
        print self.wl_max_d()[0]
        maxd12 = self.wl_max_d()[1]/self.d_wl_selector(self.wl_max_d()[0]+12)
        print maxd12
        
        combined = np.vstack((r680r630,
                              r685r630,
                              r685r655,
                              r687r630,
                              r690r630,
                              r750r800,
                              sqr685,
                              r675r690divsq683,
                              d705d722,
                              d730d706,
                              d686d710sq697,
                              maxdd720,
                              maxdd703,
                              maxd12))
                              
                              
        fluo_hdr = str('"r680r630",'+
                       '"r685r630",'+
                       '"r685r655",'+
                       '"r687r630",'+
                       '"r690r630",'+
                       '"r750r800",'+
                       '"sqr685",'+
                       '"r675r690divsq683",'+
                       '"d705d722",'+
                       '"d730d706",'+
                       '"d686d710sq697",'+
                       '"maxdd720",'+
                       '"maxdd703",'+
                       '"maxd12"')
        
        return combined, fluo_hdr
        
    def dual_peak(self):
        '''This fuction loogs for a dual peak in the red-edge region. If it's
        there it measures the depth of the feature between the two peaks. 
        UNTESTED'''
        start = self.wl_selector(640)
        end = self.wl_selector(740)
        
        d1_region = np.diff(self.values[start:end])
        #d2_region = np.diff(self.values[start:end], n=2)
        

        
        peak_finder = find_peaks_cwt(d1_region, np.arange(3,10))
        
            
        peak_wl = wavelengths[peak_finder]
        
        fluor_peaks = []

        for peak in peak_finder:
            if peak_wl[peak] == self.wl[self.wl_selector(668)]:
                print ('found flourescence peak at 668nm')
                fluor_peaks.append(peak)
                
            elif peak_wl[peak] == self.wl[self.wl_selector(735)]:
                print ('found flourescence peak at 735nm')
                fluor_peaks.append[peak]
                
            else:
                print ('unknown peak')
                
        '''if len(fluor_peaks) == 2:
            something = 'something'''
            
class load_asd():
                       
    def __init__(self, indir, output_dir):
        data_list = os.listdir(indir)
        print data_list
        
        #output_dir = os.path.join(indir,'output')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        for directory in data_list:
            parent = os.path.join(indir, directory)
            spectra_dir = os.path.join(parent, 'raw_spectra')
            reading_info_dir = os.path.join(parent, 'reading_info')
            
            sensor_name = 'ASD FieldSpec Pro'
            sensor_type = 'SPR'
            sensor_units = 'nm'
            sensor_range = [350,2500]
            
            os.chdir(reading_info_dir)
            
            reading_info_file = open('reading_atributes.txt','rb')
            
            reading_info = csv.DictReader(reading_info_file)         
            
            reading_info_array = np.empty(12)
            
            readings_list = [row for row in reading_info]
            
            
            for reading in readings_list[:]:
                
                reading_filename = str(reading['reading_id']+'.txt')
                reading_info_line = np.column_stack((reading['reading_id'],
                                                     reading['dartField'],
                                                     reading['transect'],
                                                     reading['transectPosition'],
                                                     reading['reading_type'],
                                                     reading['reading_coord_osgb_x'],
                                                     reading['reading_coord_osgb_y'],
                                                     reading['dateOfAcquisition'],
                                                     reading['timeOfAcquisition'],
                                                     reading['instrument_number'],
                                                     reading['dark_current'],
                                                     reading['white_ref']))
                #print reading_info_line
                
                if reading['reading_type']== 'REF':
                    reading_info_array = np.vstack((reading_info_array,reading_info_line))
                    #print reading_info_array
                    
                    print'*********** Loading File', reading_filename, '***********'
                    os.chdir(spectra_dir)
                    spec = np.genfromtxt(reading_filename, 
                                         delimiter=', ',
                                         skiprows=30)
                                         
                    spec = np.column_stack((spec[:,0],spec[:,1]*100))
                              
                    nir_start = 0
                    nir_end = 990
                    nir_weight = 3.5
                    nir_k = 4.9
                    nir_s =45
                    
                    swir1_start = 1080
                    swir1_end = 1438
                    swir1_weight = 8.5
                    swir1_k = 3.5
                    swir1_s = 35
                    
                    swir2_start = 1622
                    swir2_end = 2149
                    swir2_weight = 1.2
                    swir2_s = 92
                    swir2_k = 2.8
                    
                    #smoothing(perc_out, block_start, block_end, kparam, weight, sparam)
                    nir_smoothed = smoothing(spec, nir_start, nir_end, nir_k, nir_weight, nir_s)
                    swir1_smoothed = smoothing(spec, swir1_start, swir1_end, swir1_k, swir1_weight, swir1_s)
                    swir2_smoothed = smoothing(spec, swir2_start, swir2_end, swir2_k, swir2_weight, swir2_s)
                    
                    print 'Smoothed array shape', nir_smoothed.shape,swir1_smoothed.shape,swir2_smoothed.shape
                    
                    
                    
                    nir_swir_gap = interpolate_gaps(nir_smoothed,swir1_smoothed)
                    swir2_gap = interpolate_gaps(swir1_smoothed,swir2_smoothed)
                    
                 
                    spec_smoothed = np.vstack((nir_smoothed,
                                               nir_swir_gap,
                                               swir1_smoothed,
                                               swir2_gap,
                                               swir2_smoothed))
                                      
                    print 'Spec SHAPE:', spec.shape
                    
                    survey_dir = os.path.join(output_dir, directory)
                    if not os.path.exists(survey_dir):
                        os.mkdir(survey_dir)
                        
                    os.chdir(survey_dir)
                        
                    try:
                        abs470 = absorption_feature(spec_smoothed,400,518,484)
                        print abs470.abs_feature()[0]
                        abs470_ftdef = abs470.abs_feature()[0]
                        print abs470_ftdef
                        abs470_crem = abs470.abs_feature()[2]
                        if not abs470_ftdef == None:
                            np.savetxt(reading_filename[0:-4]+'_abs470_ftdef.txt',
                                       abs470_ftdef,
                                       header=abs470.abs_feature()[1],
                                       delimiter=',')
                            np.savetxt(reading_filename[0:-4]+'_abs470_crem.txt',
                                       abs470_crem,
                                       delimiter=',')
                    except:
                        pass
                    try:
                        abs670 = absorption_feature(spec_smoothed,548,800,670)
                        abs670_ftdef = abs670.abs_feature()[0]
                        abs670_crem = abs670.abs_feature()[2]
                        if not abs670_ftdef == None:
                            np.savetxt(reading_filename[0:-4]+'_abs670_ftdef.txt',
                                       abs670_ftdef,
                                       header=abs670.abs_feature()[1],
                                       delimiter=',')
                            np.savetxt(reading_filename[0:-4]+'_abs670_crem.txt',
                                       abs670_crem,
                                       delimiter=',')
                    except:
                        pass
                    
                    try:
                        abs970 = absorption_feature(spec_smoothed,880,1115,970)
                        abs970_ftdef = abs970.abs_feature()[0]
                        abs970_crem = abs970.abs_feature()[2]
                        if not abs1200_ftdef == None:
                            np.savetxt(reading_filename[0:-4]+'_abs970_ftdef.txt',
                                       abs970_ftdef,
                                       header=abs970.abs_feature()[1],
                                       delimiter=',')
                            np.savetxt(reading_filename[0:-4]+'_abs970_crem.txt',
                                       abs970_crem,
                                       delimiter=',')
                    except:
                        pass
                    
                    try:
                        abs1200 = absorption_feature(spec_smoothed,1080,1300,1190)
                        abs1200_ftdef = abs1200.abs_feature()[0]
                        abs1200_crem = abs1200.abs_feature()[2]
                        if not abs1200_ftdef == None:
                            np.savetxt(reading_filename[0:-4]+'_abs1200_ftdef.txt',
                                       abs1200_ftdef,
                                       header=abs1200.abs_feature()[1],
                                       delimiter=',')
                            np.savetxt(reading_filename[0:-4]+'_abs1200_crem.txt',
                                       abs1200_crem,
                                       delimiter=',')
                    except:
                        pass
                    try:               
                        abs1730 = absorption_feature(spec_smoothed,1630,1790,1708)
                        abs1730_ftdef = abs1730.abs_feature()[0]
                        abs1730_crem = abs1730.abs_feature()[2]
                        if not abs1730_ftdef == None:
                            np.savetxt(reading_filename[0:-4]+'_abs1730_ftdef.txt',
                                       abs1730_ftdef,
                                       header=abs1730.abs_feature()[1],
                                       delimiter=',')
                            np.savetxt(reading_filename[0:-4]+'_abs1730_crem.txt',
                                       abs1730_crem,
                                       delimiter=',')
                    except:
                        pass
                    
                    print spec_smoothed.shape
                    
                    try:
                        abs2100 = absorption_feature(spec_smoothed,2001,2196,2188)
                        abs2100_ftdef = abs2100.abs_feature()[0]
                        abs2100_crem = abs2100.abs_feature()[2]
                        if not abs2100_ftdef == None:
                            np.savetxt(reading_filename[0:-4]+'_abs2100_ftdef.txt',
                                       abs2100_ftdet,
                                       header=abs2100.abs_feature()[1],
                                       delimiter=',')
                            np.savetxt(reading_filename[0:-4]+'_abs2100_crem.txt',
                                       abs2100_crem,
                                       delimiter=',')
                    except:
                        pass
                    
                    veg_indices = Indices(spec_smoothed)
                    indices = np.column_stack((veg_indices.visnir()[0],
                                               veg_indices.nir_swir()[0],
                                               veg_indices.swir()[0]))
                    print veg_indices.visnir()[1],veg_indices.nir_swir()[1],veg_indices.swir()[1]
                    hdr = str(veg_indices.visnir()[1]+','+veg_indices.nir_swir()[1]+','+veg_indices.swir()[1])
                    
                    np.savetxt(reading_filename[0:-4]+'_indices.txt',
                               indices,
                               header=hdr,
                               delimiter=',')
                    
                    mtbvi = veg_indices.multi_tbvi()
                    np.savetxt(reading_filename[0:-4]+'_mtbvi.txt',
                               mtbvi,
                               delimiter=',')
                               
                    redge = red_edge(spec_smoothed)
                    print redge.redge_vals.shape
                    print redge.redge_vals
                    np.savetxt(reading_filename[0:-4]+'_redge.txt',
                               redge.redge_vals,
                               delimiter=',')
                    
                    fluo = fluorescence(spec_smoothed)
                    np.savetxt(reading_filename[0:-4]+'_flou.txt',
                               np.transpose(fluo.simple_ratios()[0]),
                               header = fluo.simple_ratios()[1],
                               delimiter=',')
                    
                     
                    np.savetxt(reading_filename[0:-4]+'_spec.txt',
                               spec_smoothed,
                               delimiter=',')

             
if __name__ == "__main__": 
    dir_path = os.path.dirname(os.path.abspath('...'))
    input_data = os.path.join(dir_path, 'processed')
    output_path = os.path.join(dir_path, 'output_ascii')
    loader = load_asd(input_data, output_path)
    #spectra = loader.spec_smoothed

    
             
