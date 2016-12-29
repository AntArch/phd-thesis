# -*- coding: utf-8 -*-
"""
Created on Thu May 29 12:28:08 2014

@author: david
"""

import numpy as np
import sys
import os

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo import gdal_array
from osgeo import gdalconst

from scipy.ndimage import filters

class Empty_Grid():
    def __init__(self, minx, maxx, miny, maxy, pix):
        self.rows = int((maxy-miny)/pix)
        print maxx,minx
        print 'rows', self.rows
        self.cols = int((maxx-minx)/pix)
        print 'cols',self.cols
        self.empty = np.zeros((self.rows,self.cols))
        
        print 'Array dimensions', self.empty.shape
        
        self.x_vals = np.arange(minx,maxx,pix)
        self.y_vals = np.arange(miny,maxy,pix)
        
        print self.y_vals.shape, self.x_vals.shape
        
        
        if not self.empty[1].shape == self.x_vals.shape:
            if self.empty[1].shape < self.x_vals.shape:
                print 'x empty < xvals'
                diff = self.empty[1].shape[0]-self.x_vals.shape[0]
                self.x_vals = self.x_vals[0:-diff]
            if self.empty.shape[1] > self.x_vals.shape:
                print 'x empty > xvals'
                diff = self.empty[1].shape[0]-self.x_vals.shape[0]
                newmax = self.x_vals[-1]+(diff*pix)
                self.x_vals = np.append((self.x_vals,np.arange(self.x_vals[-1],newmax,pix)))
                
        if not self.empty[0].shape[0] == self.y_vals.shape[0]:
            if self.empty[0].shape[0] < self.y_vals.shape[0]:
                print 'y empty < yvals'
                print self.empty[0].shape, self.y_vals.shape
                diff = self.empty.shape[0]-self.y_vals.shape[0]
                self.y_vals = self.y_vals[0:-diff]
            if self.empty[0].shape > self.y_vals.shape:
                print 'y empty > yvals'
                diff = self.empty[0].shape[0] - self.y_vals.shape[0]
                print self.y_vals.shape[0], self.empty[0].shape[0]
                print diff
                newmax = self.y_vals[-1]+(diff*pix)
                print y_vals[-1],newmax
                self.y_vals = np.hstack((self.y_vals,np.arange(self.y_vals[-1],newmax,pix)))  
                                                  
                                                 
class Grid_Data():
    def __init__(self, points, pix, rad):
        minx = np.min(points[:,0])
        maxx = np.max(points[:,0])
        miny = np.min(points[:,1])
        maxy = np.max(points[:,1])
        
        grid = Empty_Grid(minx,maxx,miny,maxy,pix)
        
        print 'shapes:',grid.empty.shape, grid.x_vals.shape, grid.y_vals.shape
        

        #instantiate counters
        direct = 0
        null = 0
        interpolated = 0
        
        row_ref = 0
        
        
        self.gridded = grid.empty
        
        for row in grid.y_vals:
            col_ref = 0
            #define the minimum & maximum coordinates of the row        
            cellymin = row
            cellymax = row+pix
        
            #define the centre of the cell
            cell_centy = cellymin+(pix/2)
            #use this to define search radius for interpolation later
            rad_miny = cell_centy-(pix*rad)
            rad_maxy = cell_centy+(pix*rad)
            
            # use this to find poits along the row within the radius
            # this will constrain the search space for later
            rad_row_idx = np.where((points[:,1]>=rad_miny)&
                                   (points[:,1]<=rad_maxy))
            #slice                       
            rad_row_pts = points[rad_row_idx]
            
            # find points coincident with the cells in the row
            #doing this with the search-radius subset to keep everything efficient
            row_idx = np.where((rad_row_pts[:,1]>=cellymin)&
                               (rad_row_pts[:,1]<=cellymax))
            #slice
            row_pts = rad_row_pts[row_idx]
           
            # iterate through the columns at each y value
            for column in grid.x_vals:
                # define the boundaries of th cell            
                cellxmin = column        
                cellxmax = column+pix
                
                #find points coincident with cell
                col_idx = np.where((row_pts[:,0]>=cellxmin)&
                                   (row_pts[:,0]<=cellxmax))
                
                # create an array of z values in each cell       
                col_pts = row_pts[col_idx] 
                cell_zvals = col_pts[:,2]
                         
                # get the shape of this            
                cell_znum = cell_zvals.shape[0]
                       
                # if there's only one point the becomes the z value of the cell            
                if cell_znum == 1:
                    cell_z = cell_zvals[0]
                    direct = direct+1
                    #method = 1
                    
                # if there's more than one point z = the mean of the points            
                elif cell_znum >1:
                    cell_z = np.mean(cell_zvals)
                    direct = direct+1
                    #method = 2
                    
                # if there's no points...
                else:
                    # find the centre of the cell                
                    cell_centx = cellxmin+(pix/2)
                    # define a search radius                
                    rad_minx = cell_centx-(pix*rad)
                    rad_maxx = cell_centx+(pix*rad)
                    # find the lidar points within the search radius                
                    rad_points = np.where((rad_row_pts[:,0]>=rad_minx)&
                                          (rad_row_pts[:,0]<=rad_maxx))
                    rad_zs = rad_row_pts[rad_points]
                    # work out how many points fall within this
                    rad_num = rad_zs.shape[0]
                    #if no points within the radius cell value = no data                
                    if rad_num == 0:
                        cell_z = -999
                        null = null+1
                        #method = 0
                    #otherwise pass the points to be interpolated using idw                
                    else:
                        cell_z = self.interpolate_idw(rad_points, 
                                                      cellxmin, 
                                                      cellymin, 
                                                      rad_row_pts)
                        interpolated = interpolated+1
                        #method = 3
                #print row_ref,col_ref
                self.gridded[row_ref,col_ref] = cell_z
                
                #sys.stdout.write("\r\Direct %s Interpolated %s Null %s" %(direct, interpolated, null))
                #sys.stdout.flush()
                col_ref = col_ref+1
               
            row_ref = row_ref+1    
            #stack z value into columns along row
            
            
    def distance_matrix(self,rad_points, column, row, flightline):
        #slice the input points using the indices of points in the search radius
        points = flightline[rad_points]    
        #make a matrix
        cell = np.vstack((column,row)).T
        #unfuncs....
        d0 = np.subtract.outer(points[:,0], cell[:,0])
        d1 = np.subtract.outer(points[:,1], cell[:,1])
        #return distance between points and points
        return np.hypot(d0,d1)


    def interpolate_idw(self,rad_points, column, row, flightline):
        #make the distance matrix
        d = self.distance_matrix(rad_points, column, row, flightline)
        #slice the points
        points = flightline[rad_points]
        #define distance weights using the distance matrix
        weights = 1.0/d
        #divide and resassign weights using average
        weights /= weights.sum(axis=0)
        #matrix multiplication
        cell_z = np.dot(weights.T, points[:,2])
        
        return cell_z


class WriteImage():
    def __init__(self,
                 outpath, 
                 outname,
                 xsize,
                 ysize,
                 array,
                 minx,
                 maxy,
                 pix):
        
        print 'call to WriteImage'

        data_out = np.flipud(array)
        print 'ROWS,COLS,BANDS',array.shape
        print 'Call to write image'
        os.chdir(outpath)

        print 'OUTPATH',outpath
        print 'OUTNAME',outname
        
        os.chdir(outpath)
        
        #load the driver for the format of choice
        driver = gdal.GetDriverByName("Gtiff") 
        #create an empty output file
        #get the number of bands we'll need:
        
        if len(array.shape) == 2:
            bands = 1
        elif len(array.shape)==3:
            bands = array.shape[2]
        else:
            print ('ERROR: Input array dimensions are crazy!')
            raise Exception()
        
        print 'BANDS OUT', bands
        
        #file name, x columns, y columns, bands, dtype
        out = driver.Create(outname, array.shape[1], array.shape[0], bands, gdal.GDT_Float32)
        #define the location using coords of top-left corner
        # minimum x, e-w pixel size, rotation, maximum y, n-s pixel size, rotation
        out.SetGeoTransform((minx,pix,0,maxy,0,0-pix))

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
            
            print 'Saving %s/%s' % (band,bands)
            
        else:
            while (band<=bands):
                data = data_out[:,:,band-1]
                #write values to empty array
                out.GetRasterBand(band).WriteArray( data )    
                #set the no data value
                out.GetRasterBand(band).SetNoDataValue(-999)
                #apend the statistics to dataset
                out.GetRasterBand(band).GetStatistics(0,1)  
                print 'Saving %s/%s' % (band,bands)
                band = band+1  
        out = None     
        
        print 'Processing of %s complete' % (outname)       
           


    
if __name__ == "__main__": 
    dir_path = os.path.dirname(os.path.abspath('...'))
    output_path = os.path.join(dir_path, 'output')
    if not os.path.exists(output_path):
        os.mkdir(output_path)
        
    pix = 0.5
    
    rad = 2
        
    for file in os.listdir(dir_path):
        if not os.path.isdir(file):
            points = np.genfromtxt(file, delimiter=',')
            
            print 'DATA DIMENSIONS', points.shape
            
            name = file[:-4]
            print name

            if points.shape[1]>=3:
                x_vals = points[:,0]
                y_vals = points[:,1]
                
                
                for col in range(3, points.shape[1]):
                    data = np.column_stack((x_vals,y_vals,points[:,col]))
                    print data.shape
                    interpolate = Grid_Data(data,pix,rad)
                    filtered = filters.gaussian_filter(interpolate.gridded,3)
                    print filtered.shape
                    if col == 3:
                        array = filtered
                    elif col>3:
                        array = np.dstack((array,filtered))
                    print 'array',array.shape
                    
                image = WriteImage(output_path, 
                                   name, 
                                   array.shape[1],
                                   array.shape[0],
                                   array,
                                   np.min(x_vals),
                                   np.max(y_vals),
                                   pix)
                os.chdir(dir_path)