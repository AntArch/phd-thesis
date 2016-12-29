# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 14:28:07 2013

@author: david
"""

'''SpecRGB is a wee program written for the analysis of vegetation properties 
in RGB aerial photography.. This was originally witten for the detection of 
archaeological cropamrks, but there will be other uses.
It will accept RGB images in any GDAL compliant format (e.g. *.tiff, *.jpg etc). 
The software preserves the spatial resolution of the source image, and if 
there's any spatial information attached it will preserve that too. 
It outputs the three-band original as a multi-band geotiff containing the 
derived data.

To use this you'll need Python with Numpy, Matplotlib, Pillows and GDAL installed '''

################# DEPENDANCIES ###############################################
import datetime
import numpy as np
from osgeo import gdal
from osgeo import osr

#import colorsys

import numpy.ma as ma
import os
from matplotlib import colors

from scipy.ndimage import filters

from PIL import Image


################# Classess ###################################################

#LoadImage: loads an image using gdal
class LoadImage():
    def __init__(self,infile):
        # open the dataset
        self.dataset = gdal.Open(infile) #GA_ReadOnly)
        # if there's nothign there print error
        if self.dataset is None: 
            print('BORK: Could not load file: %s' %(infile))
       # otherwise do stuff
        else:
            #get the bit depth of the source image
           
            try:
                pillow_image = Image.open(infile)
                self.bit_depth = pillow_image.bits()
                pillow_image.close()
            except:
                print ('Cant get the bit-depth of the image with pillow')
                
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

# Image RGB- the class that actually does the analysis
class ImageRGB(LoadImage):
    def __init__(self,peak_wl_r=680, peak_wl_g=550, peak_wl_b=440):
        ## bands from the source image
        print('call to image rgb')
        self.red = self.stacked[:,:,0]
        print('red', self.red.shape)
        self.green = self.stacked[:,:,1]
        print('green',self.green.shape)
        self.blue = self.stacked[:,:,2]
        print('blue',self.blue.shape)
        
        #wavelengths to use for calculating slope
        self.peak_wl_r = peak_wl_r
        self.peak_wl_g = peak_wl_g
        self.peak_wl_b = peak_wl_b
     
        #calculate  HSI
        '''self.intensity = self.red+self.blue+self.green
        print 'Intensity',self.intensity.shape
        self.hue = (self.green-self.blue)/(self.intensity-(3*self.blue))
        print  'Hue',self.hue.shape
        self.saturation = (self.intensity-(3*self.blue))/self.intensity
        print 'Saturation',self.saturation.shape'''
        
        '''self.hue = self.hsv()[:,:,0]
        self.saturation = self.hsv()[:,:,1]
        self.value = self.hsv()[:,:,2]'''
        
        self.stacked_out = np.dstack((self.red,
                                      self.green,
                                      self.blue,
                                      self.hue(),
                                      self.saturation(),
                                      self.intensity(),
                                      self.simpleratioGB(),
                                      self.simpleratioGR(),
                                      self.simpleratioRB(),
                                      self.indexNDGBI(),
                                      self.indexNDGRI(),
                                      self.indexNDRBI(),
                                      self.vari(),
                                      self.green_blue_slope(),
                                      self.green_red_slope(),
                                      self.red_blue_slope(),
                                      self.green_peak_slope(),
                                      self.red_transformation(), 
                                      self.green_transformation(),
                                      self.blue_transformation(),
                                      self.exg(),
                                      self.exr(),
                                      self.exgexr()))
        
        print('STACKED', self.stacked_out.shape)
        
    ####### Meffods ##########################################################
    
    def hsv(self):
        max_val = np.power(2,self.bit_depth)
        hsv = colors.rgb_to_hsv(self.stacked[:,:,0:3]/max_val)
        return hsv
     
    def intensity(self):
        intensity = (self.red+self.blue+self.green)/3
        print('Intensity',intensity.shape)
        return intensity
        
    def saturation(self):
        if self.intensity == 0:
            saturation = 0    
        else:
            m = np.min([self.red,self.green,self.blue])
            saturation = (1-m)/self.intensity()     
        print('saturation',saturation.shape)
        return saturation
        
        
    def hue(self):
        #Derived using Preucil colour circle
        hue = np.empty_like(self.red)
        rows=self.red.shape[0]
        columns = self.red.shape[1]
        print('call to hue method')
        print(hue.shape)
        for row in range(rows):
            column = 0
            for column in range(columns):
                if self.red[row,column] >= self.green[row,column]:
                    #r>g>b red-yellow
                    if self.green[row,column] >= self.blue[row,column]:
                        hue[row,column] = 60*(self.green[row,column]-self.blue[row,column])/(self.red[row,column]-self.blue[row,column])
                    #r>b>g magenta-red
                    elif self.blue[row,column] >= self.green[row,column]:
                        hue[row,column] = 60*(6-(self.blue[row,column]-self.green[row,column])/(self.red[row,column]-self.green[row,column]))
                    # b>r>g blue-magenta
                    elif self.blue[row,column] >= self.red[row,column]:
                        hue[row,column] = 60*(4+(self.red[row,column]-self.green[row,column])/(self.blue[row,column]-self.green[row,column]))
                else:
                    #g>r>b yellow-green
                    if self.red[row,column] >= self.blue[row,column]:
                        hue[row,column] = 60*(2-(self.red[row,column]-self.blue[row,column])/(self.green[row,column]-self.blue[row,column]))
                    #g>b>r green-cyan
                    elif self.green[row,column] >= self.blue[row,column]:
                        hue[row,column] = 60*(2+(self.blue[row,column]-self.red[row,column])/(self.green[row,column]-self.red[row,column]))
                    #b>g>r cyan-blue
                    elif self.blue[row,column] >= self.green[row,column]:
                        hue[row,column] = 60*(4-(self.green[row,column]-self.red[row,column])/(self.blue[row,column]-self.red[row,column]))
                column = column+1
       
        # deprecated methods that didn't work. delete at some point
        '''
        if self.blue > self.green:
            print 'H=360'
            h = 360
        else:
            h = 180
            print 'H=180
        
        #hue = np.tan(np.sqrt(3*(self.green-self.blue))/(2*(self.red,self.green,self.blue)))    
       
        hue = np.arctan2(np.sqrt(3*(self.green-self.blue)),
                         2*(self.red,self.green,self.blue))
            
        top = 0.5*((self.red-self.green)+(self.red-self.blue))
        bottom = np.power((np.square(self.red-self.green)+((self.red+self.blue)*(self.green-self.blue))),0.5)
        hue = 180-np.arccos(top/bottom) 
        
        print 'Hue',hue.shape
        hue = hue[1,:,:]
        print hue
        print 'HUEY',hue.shape'''
        return hue
   
   # simple ratios      
    def simpleratioGB(self):
        srGB = self.green/self.blue
        print('sRGB',srGB.shape)
        return srGB
        
    def simpleratioRB(self):
        srRB = self.red/self.blue
        print('sRB',srRB.shape)
        return srRB
        
    def simpleratioGR(self):
        srGR = self.green/self.red
        print('srGR',srGR.shape)
        return srGR
    
    #normalised differnece indices    
    def indexNDGRI(self):
        ndgri = (self.green-self.red)/(self.green+self.red)
        print('NDGRI',ndgri.shape)
        return ndgri
        
    def indexNDGBI(self):
        ndgbi = (self.green-self.blue)/(self.green+self.blue)
        print('NDGBI',ndgbi.shape)
        return ndgbi
        
    def indexNDRBI(self):
        ndrbi = (self.red-self.blue)/(self.red+self.blue)
        print('NDRBI',ndrbi.shape)
        return ndrbi
        
        
    def vari(self):
        vari = (self.green-self.red)/(self.green+self.red-self.blue)
        print('VARI',vari.shape)
        return vari
    
    #  Band transformations 
    def red_transformation(self):
        redx = self.red/self.intensity()
        print('REDX',redx.shape)
        return redx
        
    def green_transformation(self):
        greenx = self.green/self.intensity()
        print('GREENX',greenx.shape)
        return greenx
        
    def blue_transformation(self):
        bluex = self.blue/self.intensity()
        print('BLUEX',bluex.shape)
        return bluex
    
    #Slope between peaks and troughs
    def green_blue_slope(self):
        gbs = (self.green-self.blue)/(self.peak_wl_g-self.peak_wl_b)
        print('GBS',gbs.shape)
        return gbs
        
    def green_red_slope(self):
        grs = (self.green-self.blue)/(self.peak_wl_g-self.peak_wl_r)
        print('GRS',grs.shape)
        return grs
        
    def red_blue_slope(self):
        rbs = (self.red-self.blue)/(self.peak_wl_r-self.peak_wl_b)
        print('RBS',rbs.shape)
        return rbs
    
    #one of my own- angle between green peak and the troughs on either side
    def green_peak_slope(self):
        gps = np.arctan((self.green_blue_slope()-self.green_red_slope())/(1+self.green_blue_slope()*self.green_red_slope()))
        
        '''gps = 180-(np.arctan(self.green_blue_slope())+np.arctan(np.abs(self.green_red_slope())))'''
        print('GPS',gps.shape)
        return gps
        
    #Excess red and green
    def exg(self):
        print('call to exg')
        exg = (self.green_transformation()*2)-self.red_transformation()-self.blue_transformation()
        print('EXG',exg.shape)
        return exg
        
    def exr(self):
        print('call to exr')
        exr = (self.red_transformation()*1.4)-self.green_transformation()-self.blue_transformation()
        print('EXR',exr.shape)
        return exr
        
    def exgexr(self):
        print('call to exgexr')
        exgexr = self.exg()-self.exr()
        print('EXGEXR',exgexr.shape)
        return exgexr
            
class FilterImage(LoadImage, ImageRGB):
    def __init__(self, target_size=1):
        '''Work out how much we need to filter. The target_size parameter is the
        predicted size of hte features of interest'''
        
        print('Call to filter image')
        
        #Work out approximate spatial resolution
        spatial_resolution = (self.spatial[1]+np.abs(self.spatial[5]))/2
        print('spatial resolution',spatial_resolution)
        
        
        #****** We need to work out the filter size ***************************
        #how many pixels fit in our target size in the x dimension?
        self.x_num = int(target_size/self.spatial[1])
        print('xnum',self.x_num)
        #ditto y
        self.y_num = int(target_size/np.abs(self.spatial[5]))
        print('ynum',self.y_num)
        
        #if the spatial resoulution is smaller than the target size, filter it
        if spatial_resolution < target_size:
            '''we don't want a 1x1 or 1x2 filter. That would be silly, so if 
            the x or y pixel radius is too small then reassign it to be 3'''
            if self.x_num < 3:
                self.x_num=3
            if self.y_num < 3:
                self.y_num=3
            '''check to see if these are odd or even. We want them to be odd 
            (cos filter)- If they're even increment to the next odd number'''
            if self.x_num % 2 == 0:
                self.x_num = self.x_num+1
            if self.y_num % 2 == 0:
                self.y_num = self.y_num+1
            
            #filter the arrays!
            self.filtered_g = self.gaussian_filter()
            self.filtered_m = self.median_filter()
        
        else:
            print ('no need to filter')
            self.filtered_g = None
            self.filtered_m = None
            
    '''def gaussian_filter(self):
        print 'call to gaussian'
        filtered = filters.gaussian_filter(self.stacked_out, 3)
        print 'gshape', filtered.shape
        return filtered'''
        
    '''def median_filter(self):
        print 'call to median'
        print self.x_num, self.y_num
        filtered = filters.median_filter(self.stacked_out, 
                                         size=self.x_num)
        print 'mshape',filtered.shape
        return filtered'''
    
    def gaussian_filter(self):
        bands = self.stacked_out.shape[2]
        
        band = 1
        
        while band <= bands:
            if band == 1:
                #print 'Filtering 1st band'
                filtered = filters.gaussian_filter(self.stacked_out[:,:,band-1],
                                                   sigma=3)
                print('gaussian filtered band:',band)  
                band = band+1
                
            else:
                filtered = np.dstack((filtered,
                                      filters.gaussian_filter(self.stacked_out[:,:,band-1],
                                                              sigma=3)))
                print('gaussian filtered band:',band)                                             
                band = band+1
                
                                      
        return filtered
    
    
    def median_filter(self):
        bands = self.stacked_out.shape[2]
        
        band = 1
        
        while band <= bands:
            if band == 1:
                filtered = filters.median_filter(self.stacked_out[:,:,band-1],
                                                 size=(self.x_num,self.y_num))
                print('median filtered band:',band)  
                band = band+1
            else:
                filtered = np.dstack((filtered,
                                      filters.median_filter(self.stacked_out[:,:,band-1],
                                                            size=(self.x_num,self.y_num))))
                print('filtered band:',band)                                           
                band=band+1
                
                
                                      
        return filtered
                
        
        
        

class WriteImage(LoadImage, ImageRGB, FilterImage):
    def __init__(self,image,outpath,outname):
        
        #see if the image has been filtered
        '''if not self.filtered_g == None:
            #if it has then stack this with the unfiltered bands
            data_out = np.dstack((self.stacked_out,
                                  self.filtered_g,
                                  self.filtered_m))
                                  
        else:
            #otherwise carry on
            data_out = self.stacked_out'''
            
        data_out = image
        

        print('Call to write image')
        os.chdir(outpath)
        self.outpath = outpath
        print('OUTPATH', self.outpath)
        self.outname = outname
        print('OUTNAME',self.outname)
        #load the driver for the format of choice
        driver = gdal.GetDriverByName(self.driver) 
        #create an empty output file
        #get the number of bands we'll need:
        
        bands = data_out.shape[2]
        print('BANDS OUT', bands)
        
        #file name, x columns, y columns, bands, dtype
        out = driver.Create(outname, self.xsize, self.ysize, bands, gdal.GDT_Float32)
        #define the location using coords of top-left corner
        # minimum x, e-w pixel size, rotation, maximum y, n-s pixel size, rotation
        out.SetGeoTransform(self.spatial)

        srs = osr.SpatialReference()
        #get the coodrinate system using the ESPG code
        srs.SetWellKnownGeogCS(self.proj)
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
        
        print('Processing of %s complete' % (self.outname))
        
    def metacrap(self):
        os.chdir(self.outpath)
        meta_name = '%s_metadata.txt' % (self.outname[0:-3])
        
        version = 1
        
        date_time = datetime.datetime.now()
        
        blurb = '#RGB image processed using the SpecRGB Python software \
                 version %s at %s. For more information on the software email\
                 davstott@gmail.com. \n #BAND LIST: The image associated with \
                 this file is a multi-band image. The following information \
                 describes what the bands are:  \n\n "band_no",band_name"' \
                 % (version, date_time)
        
        band_names = ['0red',
                      '1green',
                      '2blue',
                      '3hue',
                      '4saturatrion',
                      '5intensity',
                      '6simpleratioGB',
                      '7simpleratioGR',
                      '8simpleratioRB',
                      '9indexNDGBI',
                      '10indexNDGRI',
                      '11indexNDRBI',
                      '12vari',
                      '13green_blue_slope',
                      '14green_red_slope',
                      '15red_blue_slope',
                      '16green_peak_slope',
                      '17red_transformation', 
                      '18green_transformation',
                      '19blue_transformation',
                      '20exg',
                      '21exr',
                      '22exgexr']
        
        if not os.path.exists(os.path.join(self.outpath,meta_name)):
            with open(meta_name, 'w') as metafile:
                metafile.writelines(blurb)
                i = 1
                for line in band_names:
                    i = i+1
                    metafile.writelines(i+','+'"'+line+'"')
                metafile.close()
                    
class RgbMakeItSo(LoadImage,ImageRGB,WriteImage,FilterImage):
    '''This is a class to make everything happen'''
    def __init__(self, indir, outdir):
        self.file_list = os.listdir(indir)
        os.chdir(indir)
        
        #set up outputs
        filtg_dir = os.path.join(outdir,'gauss')
        
        if not os.path.exists(filtg_dir):
            os.mkdir(filtg_dir)
        
        filtm_dir = os.path.join(outdir,'median')
        
        if not os.path.exists(filtm_dir):
            os.mkdir(filtm_dir)
            
        #instantiate stuff          
        for image in self.file_list:
            os.chdir(indir)
            try:
                LoadImage.__init__(self,
                                   image)
                ImageRGB.__init__(self)
                FilterImage.__init__(self)
                WriteImage.__init__(self,
                                    self.stacked_out,
                                    outdir,
                                    image)
                WriteImage.__init__(self,
                                    self.filtered_g,
                                    filtg_dir,
                                    'g'+image)
                WriteImage.__init__(self,
                                    self.filtered_m,
                                    filtm_dir,
                                    'm'+image)
                                    
            except:
                print ('NOT A VALID IMAGE FILE')
            
        
        
        
###################### MAIN BLAH BLAH BLAH ###################################
if __name__ == "__main__": 
    dir_path = '/home/dav/data/rectified_verticals/sprgb'
    #dir_path = os.path.dirname(os.path.abspath('...'))
    indir = os.path.join(dir_path,'Images')  
    outdir = os.path.join(dir_path,'Output')
    
    RgbMakeItSo(indir,outdir)
        
