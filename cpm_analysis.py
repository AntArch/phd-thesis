# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:05:31 2014

@author: david
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:30:37 2014

@author: david
"""

import numpy as np
from osgeo import gdal
from osgeo import osr



import numpy.ma as ma
import os

import matplotlib.pyplot as plt

from sklearn.decomposition import RandomizedPCA
from sklearn.decomposition import PCA

from sklearn import preprocessing


import emd

import Image

import fiona
from shapely.geometry import shape
from shapely.geometry import asPoint
from scipy import stats


#LoadImage: loads an image using gdal
class LoadImage():
    def __init__(self,infile):
        
        # open the dataset
        self.image_name = infile[:-4]
        self.dataset = gdal.Open(infile) #GA_ReadOnly)
        # if there's nothign there print error
        #self.stacked = None
        if self.dataset is None: 
            print 'BORK: Could not load file: %s' %(infile)
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
            print 'BANDS:',bands
            #get the geotransform Returns a list object. This is standard GDAL ordering:
                #spatial[0] = top left x
                #spatial[1] = w-e pixel size
                #spatial[2] = rotation (should be 0)
                #spatial[3] = top left y
                #spatial[4] = rotation (should be 0)
                #spatial[5] = n-s pixel size
            self.spatial = self.dataset.GetGeoTransform()

            #print some stuff to console to show  we're paying attention
            print 'Found raster in %s format. Raster has %s bands' %(self.driver,bands)
            print 'Projected as %s' %(self.proj)
            print 'Dimensions: %s x %s' %(self.xsize,self.ysize)
            
            #instantiate a counter
            count = 1
            
            #OK. This is the bit that catually loads the bands in in a while loop
            # Loop through bands as long as count is equal to or less than total
            while (count<=bands):
                print 'BANDS less than COUNT'
                #show that your computer's fans are whining for a reason
                print 'Loading band: %s of %s' %(count,bands)
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
                
            self.coords_matrix = self.coords()
            print self.coords_matrix.shape
            #print self.coords_matrix
        
    #******************** MEFFODS ******************************************
    def coords(self):
        '''This gets the geographic coordinates of each cell in the raster and 
        returns a list of tuples containing the y and x array references for each pixel
        and the geographic x and y coordinates for each pixel'''
        
        print 'call to coords'
        
        #get the shape of the array
        matrix_dims = self.stacked.shape
        print matrix_dims
        #get the number of rows
        rows = matrix_dims[0]-1
        print rows
        #get the number of columns
        columns = matrix_dims[1]-1
        print columns
        
        x_coords = np.zeros(shape=(matrix_dims[0],matrix_dims[1]))
        y_coords = np.zeros(shape=(matrix_dims[0],matrix_dims[1]))

        
        #instantiate a counter
        row = 0
        
        #fruity loops
        for row in range(matrix_dims[0]):
            #increment counter
            column = 0
            for column in range(matrix_dims[1]):
                #print row, column
                xgeo = self.spatial[0] + (column*self.spatial[1])+(row*self.spatial[2])
                x_coords[row,column] = xgeo
                ygeo = self.spatial[3] + (column*self.spatial[4])+(row*self.spatial[5])
                y_coords[row,column] = ygeo
                column=column+1
        
        print x_coords.shape, y_coords.shape 
        print np.min(x_coords), np.min(y_coords)
        return np.dstack((x_coords,
                          y_coords,
                          np.zeros(shape=(matrix_dims[0],matrix_dims[1]))))
        
    def coord_test(self):
        print 'call to coord test'

        print self.coords_list[0]
        print self.coords_list[-1]

class ImageAnalysis():      
    def __init__(self,image_dir,image,mask_dir,mask,plot_dir, band_names):
        os.chdir(image_dir)
        loadimage = LoadImage(image)
        print loadimage
        self.multiband_image = loadimage.stacked
        self.name = loadimage.image_name
        
        self.plot_dir = plot_dir
        
        os.chdir(mask_dir)
        loadmask = LoadImage(mask)
        print loadmask
        self.classes = loadmask.stacked
        
        self.arc = self.multiband_image[np.where(self.classes==1)]
        self.bac = self.multiband_image[np.where(self.classes==2)]
        
        print 'ARC',self.arc.shape
        print 'BAC',self.bac.shape
        
        self.band_names = band_names
        
        if len(self.multiband_image.shape) == 3:
            self.bands = self.multiband_image.shape[2]
        elif len(self.multiband_image.shape) == 2:
            self.bands = 1
        
        
    def hist_compare(self):
        
        
        print 'BANDS',self.bands
        
        bins = 32
        
        emd_list = []
        
        out_dir = os.path.join(self.plot_dir,self.name)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            
        hist_dir = os.path.join(out_dir,'histograms')
        if not os.path.exists(hist_dir):
            os.mkdir(hist_dir)
                
                

        archaeology = self.arc
        background = self.bac
        
        #print archaeology.shape
        #print background.shape
        
        minima = np.min(archaeology)
        if minima > np.min(background):
            minima = np.min(background)
            
        maxima = np.min(archaeology)
        if maxima < np.max(background):
            maxima = np.max(background)
            
        hist_arch = np.histogram(archaeology,
                                 bins=bins, 
                                 range=(minima,maxima))
        hist_back = np.histogram(background,
                                 bins=bins, 
                                 range=(minima,maxima))
                                 
        #print hist_arch[0]
        #print hist_back[0]
        
        print 'Totals'             
        print 'hist_arch', np.sum(hist_arch[0])
        print archaeology.shape                          
        
        #print hist_arch[0].shape
        #print hist_arch[0].shape
                                 
                                 
        hist_arch_norm = np.true_divide(hist_arch[0],archaeology.shape)
        hist_back_norm = np.true_divide(hist_back[0],background.shape)
        
        #hist_arch_norm = hist_arch
        #hist_back_norm = hist_back           
        
        #print hist_arch_norm
        

        
        #print x_vals.shape
        
        os.chdir(hist_dir)
      
                                             
      
        
        sum_of_difference = np.sum(np.abs(hist_arch_norm-hist_back_norm))
                                          
        print sum_of_difference       

        contrast_emd = emd.emd(range(bins),range(bins),hist_arch_norm, hist_back_norm)          
        print 'EMD',contrast_emd
        emd_list.append(contrast_emd)
            
        emd_comp = np.array(emd_list)
        print emd_comp.shape

        
        os.chdir(self.plot_dir)
        np.savetxt(self.name+'_emd.txt',emd_comp, delimiter=',')
        
    def anova(self):

        f = []        
        
        p = []
        
        
        for band in range(self.bands):
            archaeology = self.arc[:,band]
            background = self.bac[:,band]
            
            anova = stats.f_oneway(archaeology, background)
            f.append(anova[0])
            p.append(anova[1])
            
        f_list = np.array(f)
        p_list = np.array(p)
        
        pos = np.arange(self.bands)
        #plt.xlim([0,6])
        plt.barh(pos,f_list, color='#3A01DF')
        plt.margins(0.01)
        
        plt.yticks(pos+0.5, band_names, style ='italic')
        
        plt.xlabel("Contrast (ANOVA F)",fontweight='bold')
        #plt.ylabel('RGB')
        plt.title(self.name, fontweight='bold')
        
        
        fig=plt.gcf()
        fig.subplots_adjust(left=0.18)
        
        
        plt.savefig(self.name+'_ANOVA_F')
        plt.close()
        
        os.chdir(self.plot_dir)
        np.savetxt(self.name+'_anova_f.txt',f_list, delimiter=',')
        
        
        pos = np.arange(self.bands)
        #plt.xlim([0,6])
        plt.barh(pos,p_list, color='#3A01DF')
        plt.margins(0.01)
        
        plt.yticks(pos+0.5, band_names, style ='italic')
        
        plt.xlabel("Contrast (ANOVA P)",fontweight='bold')
        #plt.ylabel('RGB')
        plt.title(self.name, fontweight='bold')
        
        
        fig=plt.gcf()
        fig.subplots_adjust(left=0.18)
        
        
        plt.savefig(self.name+'_ANOVA_P')
        plt.close()
        
        os.chdir(self.plot_dir)
        np.savetxt(self.name+'_anova_p.txt',p_list, delimiter=',')
            
            
            
    def pca(self):
        print 'PCA'
        
        
        
        
        data = np.vstack((self.arc[:,0:2],self.bac[:,0:2]))
        
        #bad = data[np.where(np.isfinite(data))]
        
        masked = np.ma.array(data, 
                             mask = np.isfinite(data))        
        
        #print data.shape
        #print bad.shape, 'BAD'
        
        fixed = np.ma.fix_invalid(masked)
        
        imp = preprocessing.Imputer(missing_values='NaN', 
                                    strategy='mean',
                                    axis=0)
        
        imp.fit_transform(fixed.data)

        pca = PCA(n_components=3)        
        rgb_pca = pca.fit_transform(fixed.data)
        
        print rgb_pca.shape
        print rgb_pca[0],rgb_pca[1]
        
        
        
        
        plt.scatter(rgb_pca[:self.arc.shape[0],0],rgb_pca[:self.arc.shape[0],1],color='red',alpha=0.5)
        plt.scatter(rgb_pca[self.arc.shape[0]:,0],rgb_pca[self.arc.shape[0]:,1],color='blue',alpha=0.5)
        
        #plt.show()
        
        os.chdir(self.plot_dir)
        
        plt.savefig('1111pca'+self.name)
        
        
        plt.close()
        
    def ttest(self):
        
        t =[]
        p = []
              
        
        archaeology = self.arc
        background = self.bac
        
        print 'ttestloop'
        print archaeology.shape
        print background.shape
        
        t_test = stats.ttest_1samp(archaeology, background)
        #t_test = stats.ttest_ind(archaeology, background)
        t.append(t_test[0])
        p.append(t_test[1])
            
        t_list = np.array(t)
        print 'TSHAPE', t_list.shape
        
        
        p_list = np.array(p)
        print 'PSHAPE', p_list.shape
        
               
        os.chdir(self.plot_dir)
        np.savetxt(self.name+'_t.txt',np.mean(t_list, axis=1), delimiter=',')
        np.savetxt(self.name+'_p.txt',np.mean(p_list, axis=1), delimiter=',')
            
        #np.savetxt(self.name+'_t.txt',t_list, delimiter=',')
        #np.savetxt(self.name+'_p.txt',p_list, delimiter=',')
        
    
        

if __name__ == "__main__": 
    dir_path = os.path.dirname(os.path.abspath('...'))
    
    plot_dir = os.path.join(dir_path,'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
        
    band_names = ('intensity',
                  'peak_start',
                  'peak_end',
                  'peak_location',
                  'peak_width',
                  'max_intensity',
                  'peak_sum',
                  'shoulder_location',
                  'shoulder_intensity')

    dir_list = os.listdir(dir_path)
    for directory in dir_list:
        if not directory == 'plots':
            current = os.path.join(dir_path,directory)
            image_dir = os.path.join(current, 'image')
            class_dir = os.path.join(current, 'grid_output')
              
            img_list = os.listdir(image_dir)
            class_list = os.listdir(class_dir)    
        
            for image in img_list:
                    #img = os.path.join(image_dir,image)
                    #print image
                    img = image
                
                    for classification in class_list:
                        clim = classification
                        #print clim
                        blah = ImageAnalysis(image_dir,
                                             img,class_dir,clim, 
                                             plot_dir,
                                             band_names)
                        blah.hist_compare()
                        
                        #blah.pca()
                      
                        blah.ttest()