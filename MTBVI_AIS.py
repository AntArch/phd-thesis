# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 20:50:47 2016

@author: dav
"""

import numpy as np
#import spec_gdal4 as spg
from osgeo import gdal
import os
import csv
#import h5py
import datetime
import numpy.ma as ma
#from StringIO import StringIO
#import shapely
#import r2py

from osgeo import gdal_array
from osgeo import gdalconst
from osgeo.gdalconst import *


from osgeo import ogr
from osgeo import osr



import matplotlib.pyplot as plt

class MtbviImage():
    def __init__(self,
                 image,
                 wls):
                     
        self.wl = wls
        
        self.image = image
        
        

    def multi_tbvi (self, red_start=630, red_end=750, ir_start=700, ir_end=850):
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
        left = self.image[:,:,red_l:red_r]
        right = self.image[:,:,ir_l:ir_r]
        
        out = None
                              
        #set up counter
        l = 0    
        #loop throught the values in the red
        for lband in range(left.shape[2]):
            l = left[:,:,lband]
            
            #then calculate the index with each wl in the NIR
            for rband in range(right.shape[2]):
                r = right[:,:,rband]
                mtbvi = (r-l)/(r+l)

                if lband == 0 and rband ==0:
                    out = mtbvi
                else:
                    out = np.dstack((out, mtbvi))
                    
        return out, self.wl[red_l:red_r],self.wl[ir_l:ir_r]
            
class LoadImage():
    def __init__(self,infile):
        
        # open the dataset
        self.image_name = infile[:-4]
        self.dataset = gdal.Open(infile) #GA_ReadOnly)
        # if there's nothign there print error
        #self.stacked = None
        if self.dataset is None: 
            print('BORK: Could not load file: %s' %(infile))
       # otherwise do stuff
        else:
            #get the bit depth of the source image
           
            '''try:
                pillow_image = Image.open(infile)
                self.bit_depth = pillow_image.bits()
                pillow_image.close()
            except:
                print ('Cant get the bit-depth of the image with pillow')'''
                
            #get the format
            self.driver = self.dataset.GetDriver().ShortName
            #get the x dimension
            self.xsize = self.dataset.RasterXSize
            #get the y dimension
            self.ysize = self.dataset.RasterYSize
            #get the projection
            self.proj = self.dataset.GetProjection()
            #get the number of bands
            bands = self.dataset.RasterCount
            print('BANDS:',bands)
            #get the geotransform Returns a list object. This is standard GDAL ordering:
                #spatial[0] = top left x
                #spatial[1] = w-e pixel size
                #spatial[2] = rotation (should be 0)
                #spatial[3] = top left y
                #spatial[4] = rotation (should be 0)
                #spatial[5] = n-s pixel size
            self.spatial = self.dataset.GetGeoTransform()

            #print some stuff to console to show  we're paying attention
            print('Found raster in %s format. Raster has %s bands' %(self.driver,bands))
            print('Projected as %s' %(self.proj))
            print('Dimensions: %s x %s' %(self.xsize,self.ysize))
            
            #instantiate a counter
            count = 1
            
            #OK. This is the bit that catually loads the bands in in a while loop
            # Loop through bands as long as count is equal to or less than total
            while (count<=bands):
                print('BANDS less than COUNT')
                #show that your computer's fans are whining for a reason
                print('Loading band: %s of %s' %(count,bands))
                #get the band
                band = self.dataset.GetRasterBand(count)
                # load this as a numpy array
                
                #mask the no data values
                data_array = band.ReadAsArray()
                data_array = ma.masked_where(data_array == 0, data_array)
                data_array = data_array.filled(-999)
                data_array = data_array.astype(np.float32, copy=False)
                # close the band object
                band = None
                #this bit stacks the bands into a combined numpy array
                #if it's the first band copy the array directly to the combined one
                if count == 1:
                    self.stacked = data_array
                #else combine these 
                else:
                    self.stacked = np.dstack((self.stacked,data_array))            
                # increment the counter
                count = count+1
                
            #self.coords_matrix = self.coords()
            #print self.coords_matrix.shape
            #print self.coords_matrix
                
def writeimage(outpath, 
               outname,
               image,
               spatial):

    data_out = image
    print('ROWS,COLS',image.shape)
    print('Call to write image')
    os.chdir(outpath)
    print('OUTPATH',outpath)
    print('OUTNAME',outname)
    
    #load the driver for the format of choice
    driver = gdal.GetDriverByName("Gtiff") 
    #create an empty output file
    #get the number of bands we'll need:
    
    bands = image.shape[2]
    print('BANDS OUT', bands)
    
    #file name, x columns, y columns, bands, dtype
    out = driver.Create(outname, image.shape[1], image.shape[0], bands, gdal.GDT_Float32)
    #define the location using coords of top-left corner
    # minimum x, e-w pixel size, rotation, maximum y, n-s pixel size, rotation
    out.SetGeoTransform(spatial)

    srs = osr.SpatialReference()
    #get the coodrinate system using the ESPG code
    srs.SetWellKnownGeogCS("EPSG:27700")
    #set pstackedstackedstackedtojection of output file 
    out.SetProjection(srs.ExportToWkt())
    
    band = 1
    
    if bands == 1:
        out.GetRasterBand(band).WriteArray(data_out)
        #set the no data value
        out.GetRasterBand(band).SetNoDataValue(-999)
        #apend the statistics to dataset
        out.GetRasterBand(band).GetStatistics(0,1)
        
        print('Saving %s/%s' % (band,bands))
        
    else:
        while (band<=bands):
            data = data_out[:,:,band-1]
            #write values to empty array
            out.GetRasterBand(band).WriteArray( data )    
            #set the no data value
            out.GetRasterBand(band).SetNoDataValue(-999)
            #apend the statistics to dataset
            out.GetRasterBand(band).GetStatistics(0,1)  
            print('Saving %s/%s' % (band,bands))
            band = band+1  
    out = None     
    
    print('Processing of %s complete' % (outname))       
    
    return outname
    
def dirhandler(dir, outdir):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    image_dir = os.path.join(dir,'image')
    
    wavelengths = np.genfromtxt(os.path.join(dir,'wavelengths','wavelengths.txt'))
    
    for image in os.listdir(image_dir):
        if image[-3:] == 'tif':
            image_in = LoadImage(os.path.join(image_dir,image))
            m = MtbviImage(image_in.stacked, wavelengths)
            
            mtbvi = m.multi_tbvi()
            
            writeimage(outdir, 
                       image[:-4]+'_mtbvi.tif',
                       mtbvi[0],
                       image_in.spatial)
                       
            np.savetxt(os.path.join(outdir,image[:-4]+'_lWL.txt'),
                       mtbvi[1],
                       delimiter=',')
            
            np.savetxt(os.path.join(outdir,image[:-4]+'_rWL.txt'),
                       mtbvi[2],
                       delimiter=',')
                              

if __name__ == '__main__':
    
    #rootdir = '/home/dav/data/temp/test/done'
    
    rootdir='/home/dav/data/temp/optimise_tmp/mtbvi'
    
    outdir= '/home/dav/data/temp/test/optimise_tmp/mtbvi_done'
    
    for dir in os.listdir(rootdir):
        if dir == 'subsets':
            subsets_dir = os.path.join(rootdir,'subsets')
            
            print(subsets_dir)
            
            for d1 in os.listdir(subsets_dir):
                d = os.path.join(subsets_dir, d1)
                
                for d2 in os.listdir(d):
                    ot = os.path.join(outdir, 'subsets', d2, 'output')
                    
                    dirhandler(os.path.join(d,d2),ot)
                    
        '''else:
            ot = os.path.join(outdir, dir, 'output')
            dirhandler(os.path.join(rootdir,dir),ot)'''
            
            
        
        
        
        
        
    