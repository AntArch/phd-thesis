# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 23:33:37 2016

@author: dav
"""

import numpy as np
from osgeo import gdal
from osgeo import osr


import numpy.ma as ma
import os

import matplotlib.pyplot as plt


#import Image

import fiona
from shapely.geometry import shape
from shapely.geometry import asPoint
from scipy import stats


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
    
    print (image.shape)
    if len(image.shape)>2:
        bands = image.shape[2]
    else:
        bands =1
    print('BANDS OUT', bands)
    
    #file name, x columns, y columns, bands, dtype
    out = driver.Create(outname, image.shape[1], image.shape[0], bands, gdal.GDT_Float32)
    #define the location using coords of top-left corner
    # minimum x, e-w pixel size, rotation, maximum y, n-s pixel size, rotation
    out.SetGeoTransform(spatial)

    srs = osr.SpatialReference()
    #get the coodrinate system using the ESPG code
    srs.SetWellKnownGeogCS("EPSG:4277")
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
    
#Class to load the raster image and create an empty raster of the same dimensions
#LoadImage: loads an image using gdal
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
        
def mtbvi_stats(image, classes, lwls, rwls):
    out= np.zeros((image.shape[2],2))
    
    cols = lwls.shape[0]
    
    rows = rwls.shape[0]
    
    for i in range(image.shape[2]):
        values=image[:,:,i]
        
        arc = values[np.where(classes==1)]
        bac = values[np.where(classes==2)]
        
        #print (arc)
        #print(bac)
        
        t = stats.ttest_ind(arc,bac, equal_var=False)
        
        out[i,0]=t[0]
        out[i,1]=t[1]
        
    mt_grid = np.reshape(out[:,0],(cols,rows))
    
    finite = np.isfinite(out[:,0])
    
    best = np.argmax(np.abs(out[finite,0]))
    
    lwl_otsh = np.repeat(lwls,out.shape[0]/cols)[finite]
    
    rwl_otsh = np.tile(rwls,out.shape[0]/rows)[finite]
    
    
    best = np.array([best,lwl_otsh[best],rwl_otsh[best],out[best,0],out[best,1]])
    
    '''best = np.argmax(np.abs(out[:,0]))
    
    lwl_otsh = np.repeat(lwls,out.shape[0]/cols)
    
    rwl_otsh = np.tile(rwls,out.shape[0]/rows)
    
    
    best = np.array([best,lwl_otsh[best],rwl_otsh[best],out[best,0],out[best,1]])'''
    print (best)
    
    
        
    return mt_grid, best
    
if __name__ == '__main__':
    
    root_dir = '/home/dav/data/temp/test/mtbvi/subsets'
    
    for dir in os.listdir(root_dir):
        
        d = os.path.join(root_dir,dir)
        
        outdir = os.path.join(d,'results')
        
        classes  = LoadImage(os.path.join(d,'grid_output','mask.tif')).stacked
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        image_dir = os.path.join(d,'output')
        
        lwl = np.genfromtxt(os.path.join(image_dir,dir+'_lWL.txt'))
        rwl = np.genfromtxt(os.path.join(image_dir,dir+'_rWL.txt'))
        
        corr_dir = '/home/dav/data/temp/test/correlation'
        
        for i in os.listdir(image_dir):
            if i[-3:]=='tif':
                im = LoadImage(os.path.join(image_dir,i))
                
                image = im.stacked
                
                spatial = im.spatial
                
                stats_ot = mtbvi_stats(image,classes,lwl,rwl)
                
                plt.imshow(np.rot90(stats_ot[0],k=1), 
                           extent=[lwl[0],lwl[-1],rwl[0],rwl[-1]],
                           interpolation='none',
                           cmap='gnuplot2')
                plt.colorbar()
                plt.xlabel('RED $\lambda  (nm)$')
                plt.ylabel('NIR  $\lambda (nm)$')
                
                i_name =dir.split('_')[0]
                i_date =dir.split('_')[1]
                
                if i_name == 'ddt1':
                    tname = 'Didington Transect 1 %s' %(i_date)
                if i_name == 'ddch':
                    tname = 'Didington Clay Field %s' %(i_date)
                if i_name == 'hhcc':
                    tname = 'Harnhill Cherry Copse %s' %(i_date)
                if i_name == 'hhqf':
                    tname = 'Harnhill Quarry Field %s' %(i_date)
                
                
                
                plt.title(tname+':\n NDVI variable wavelength: Contrast')
                plt.tight_layout()
            
                plt.savefig(os.path.join(outdir,'_'+dir+'_mtbvi.png'))
                
                plt.show()
                
                plt.close()
                
                np.savetxt(os.path.join(outdir,'_'+dir+'_bestmtbvi.txt'),stats_ot[1])
                
                '''writeimage(os.path.join(corr_dir,dir),
                           str(dir+'_'+str(int(stats_ot[1][1]))+'_'+str(int(stats_ot[1][2]))+'_mtbvi.tif'),
                           image[:,:,int(stats_ot[1][0])],
                           spatial)'''
                
                
                
                
            
        
        
        
        
        
        
    
    

        