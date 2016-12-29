# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 20:32:25 2016

@author: dav
"""

import os

import numpy as np

import matplotlib.pyplot as plt


################# GLOBALS #####################################################
index_names = ['SR 700 800',
               'NDVI 760 694',
               'NDVI 805 695',
               'NDVI 800 700',
               'NDVI 750 705',
               'RDVI',
               'SAVI',
               'MSAVI',
               'MSR',
               'MSRVI',
               'MDVI',
               'TVI',
               'MTVI',
               'MTVI2',
               'VOG1',
               'VOG2',
               'PRSI',
               'PRI',
               'SIPI',
               'CARI',
               'MCARI1',
               'MCARI2',
               'NPCI',
               'NPQI',
               'CRI1',
               'CRI2',
               'ARI1',
               'ARI2',
               'WBI',
               'NDWI',
               'MSI',
               'NDII',
               'NDNI',
               'NDLI']
                            
fluo_names = ['$680/630$',
              '$685/630$',
              '$685/655$',
              '$687/630$',
              '$690/630$',
              '$750/800$',
              '$685^{2}/(675-690)$',
              '$(685-690)/683^{2}$',
              "$f'705/f'722$",
              "$f'730/f'706$",
              "$(f'688-f'710)/f'697^{2}$",
              "$max f'/f'720$",
              "$max f'/f'703$",
              "$max f'/ max f'+12$"]
                           
ft470_names = ['Refined start 470nm',
               'Refined end 470nm',
               'Minima $\lambda$  470nm',
               'Minima Reflectance 470nm',
               'Feature Area 470nm',
               'Left TBVI 470nm',
               'Right TBVI 470nm',
               'Gradient 470nm',
               'Area 470nm',
               'Maxima reflectance 470nm',
               'Maxima $\lambda$ 470nm',
               'Area Left 470nm',
               'Area Right 470nm']
ft670_names = ['Refined start 670nm',
               'Refined end 670nm',
               'Minima $\lambda$ 670nm',
               'Minima Reflectance 670nm',
               'Feature Area 670nm',
               'Left TBVI 670nm',
               'Right TBVI 670nm',
               'Gradient 670nm',
               'Area 670nm',
               'Maxima reflectance 670nm',
               'Maxima $\lambda$ 670nm',
               'Area Left 670nm',
               'Area Right 670nm']

ft970_names = ['Refined start 970nm',
               'Refined end 970nm',
               'Minima $\lambda$ 970nm',
               'Minima Reflectance 970nm',
               'Feature Area 970nm',
               'Left TBVI 970nm',
               'Right TBVI 970nm',
               'Gradient 970nm',
               'Area 970nm',
               'Maxima reflectance 970nm',
               'Maxima $\lambda$ 970nm',
               'Area Left 970nm',
               'Area Right 970nm']

ft1200_names = ['Refined start 1200nm',
               'Refined end 1200nm',
               'Minima $\lambda$ 1200nm',
               'Minima Reflectance 1200nm',
               'Feature Area 1200nm',
               'Left TBVI 1200nm',
               'Right TBVI 1200nm',
               'Gradient 1200nm',
               'Area 1200nm',
               'Maxima reflectance 1200nm',
               'Maxima $\lambda$ 1200nm',
               'Area Left 1200nm',
               'Area Right 1200nm']

index_colour = '#8c510a'
redge_colour = '#bf812d'
fluo_colour = '#dfc27d'
mtbvi_colour = '#f6e8c3'
ft470_colour = '#c7eae5'
ft670_colour = '#80cdc1'
ft970_colour = '#35978f'
ft1200_colour = '#01665e'

index_names = np.asarray(index_names)
fluo_names = np.asarray(fluo_names)
ft470_names = np.asarray(ft470_names)
ft670_names = np.asarray(ft670_names)
ft970_names = np.asarray(ft970_names)
ft1200_names = np.asarray(ft1200_names)


def parse_names(s):
    if s == 'ddt1':
        n = 'Diddington Transect 1'
    
    elif s == 'ddcf':
        n = 'Diddington Clay Field'
        
    elif s == 'hhcc':
        n = 'Harnhill Cherry Copse'
        
    elif s == 'hhqf':
        n = 'Harnhill Quarry Field'
        
    return n


class Plotting():
    
    def __init__(self, name, outdir):
        
        name = name.split('_')
        
        self.site = name[0]
        
        self.date= '%s-%s-%s' %(name[1][:4],name[1][4:6],name[1][6:8])
        
        self.outdir = outdir
        
    '''def appendix_plots(self,
                       indices,
                       redge,
                       fluo,
                       ft670):
                           
        
        colours = ['#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#4dac26',
                   '#d01c8b',
                   '#d01c8b',
                   '#d01c8b',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#b8e186',
                   '#f1b6da']
                   
        xes = np.arange(0,indices.shape[0],1)
                    
        fig=plt.figure()
        ax1 = plt.subplot()
        ax1.bar(xes, indices[:,1], color=colours)
        ax1.margins(0.01)
        plt.xticks(xes+0.5, index_names, style ='italic', rotation='vertical', size='smaller')
        plt.xlabel('Vegetation Indices')
        plt.ylabel('Contrast (T value)')
        plt.title(parse_names(self.site)+' '+self.date+':\n Vegetation Indices')
        plt.tight_layout()
        
        
        lgd_veg = plt.Rectangle((0,0),1,1, fc='#4dac26')
        lgd_redge = plt.Rectangle((0,0),1,1, fc= '#d01c8b')
        lgd_pig = plt.Rectangle((0,0),1,1, fc='#b8e186')
        lgd_chem = plt.Rectangle((0,0),1,1, fc='#f1b6da')
        
        box=ax1.get_position()
        ax1.set_position([box.x0,box.y0,box.width*0.8, box.height])
        
        plt.legend([lgd_veg,lgd_redge,lgd_pig,lgd_chem],
                   ['Generic','Red Edge','Pigments','Biochemicals'],
                   loc='center left',
                   bbox_to_anchor=(1,0.5),
                   prop={'size':10})
        plt.savefig(os.path.join(self.outdir,self.site+'_'+self.date+'_indices.png'))
        
        plt.close()
        
        fluo = np.vstack((fluo, redge))
        
        names = np.append(fluo_names,'$REIP$')
        
        #names = np.vstack((np.asarray(fluo_names),'$REIP$'))
        
        xes = np.arange(0,fluo.shape[0],1)
        plt.bar(xes, fluo[:,1], color='#4dac26')
        plt.margins(0.01)
        plt.xticks(xes+0.5, names, style ='italic',rotation='vertical',size='smaller')
        plt.xlabel('Fluorescence Indices & REIP')
        plt.ylabel('Contrast (T value)')
        plt.title(parse_names(self.site)+' '+self.date+':\n Fluorescence Indices & REIP')
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir,self.site+'_'+self.date+'_fluo.png'))
        plt.close()
        
        print (ft670)
        
        xes = np.arange(0, ft670.shape[0], 1)
                    
        fig=plt.figure()
        ax1 = plt.subplot()
        ax1.bar(xes, ft670[:,1], color='#4dac26')
        ax1.margins(0.01)
        plt.xticks(xes+0.5, ft670_names, style ='italic', rotation='vertical', size='smaller')
        plt.xlabel('Metric')
        plt.ylabel('Contrast (T)')
        plt.title(parse_names(self.site)+' '+self.date+':\n Continuum removal metrics: '+'670nm feature')
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir,self.site+'_'+self.date+'_ft670nm.png'))
        
        plt.close()'''
        
        
                           
    def summary_plot_data(self,
                          indices,
                          redge,
                          fluo,
                          mtbvi,
                          ft470,
                          ft670,
                          ft970,
                          ft1200):
        data = None
        
        names = None
                          
 
        
        #best_indices = np.where(np.abs(indices[:,1])>4 and indices[:,2]<0.6)
        print (indices.shape)
        best_indices = np.logical_and(np.abs(indices[:,0])>5,indices[:,1]<0.05)
        
        print(best_indices)
        
        #best_fluo = np.where(np.abs(fluo[:,1])>4 and fluo[:,2]<0.6)
        best_fluo = np.logical_and(np.abs(fluo[:,0])>5,fluo[:,1]<0.05)         
                                     
        #best_ft670 = np.where(np.abs(fluo[:,1])>4 and fluo[:,2]<0.6)
        best_ft470 = np.logical_and(np.abs(ft470[:,0])>5,ft470[:,1]<0.05)
        best_ft670 = np.logical_and(np.abs(ft670[:,0])>5,ft670[:,1]<0.05)
        best_ft970 = np.logical_and(np.abs(ft970[:,0])>5,ft970[:,1]<0.05)
        best_ft1200 = np.logical_and(np.abs(ft1200[:,0])>5,ft1200[:,1]<0.05)
        
        indices = indices[:,0]
        #redge = redge[0]
        fluo = fluo[:,0]
        ft470 = ft470[:,0]
        ft670 = ft670[:,0]
        ft970 = ft970[:,0]
        ft1200 = ft1200[:,0]
        
        print (indices.shape, fluo.shape)
        print ('*&*&*&*&*&*&*&*&*&*')
        print (indices[best_indices].size)
        print (fluo[best_fluo].size)
        

        if not indices[best_indices].size == 0:
            if not data is None:
                data = np.append(data,indices[best_indices])
                names = np.append(names,index_names[best_indices])
            else:                           
                data = indices[best_indices]
                names = np.array((index_names[best_indices]))
                
            print ('I', data.shape, names.shape)
                
        if not fluo[best_fluo].size == 0:
            print ('F', fluo[best_fluo].shape)
            if not data is None:
                print ('D', data.shape)
        
                data = np.append(data,fluo[best_fluo])
                names = np.append(names,fluo_names[best_fluo])
            else:                           
                data = fluo[best_fluo]
                names = np.array((fluo_names[best_fluo]))
                
            print ('F', data.shape, names.shape)
                
        if not ft470[best_ft470].size == 0:
            if not data is None:
                data = np.append(data,ft470[best_ft470])
                names = np.append(names,ft470_names[best_ft470])
            else:                           
                data = ft470[best_ft470]
                names = np.array(ft470_names[best_ft470])
                
            print ('F47', data.shape, names.shape)
                
        if not ft670[best_ft670].size == 0:
            if not data is None:
                data = np.append(data,ft670[best_ft670])
                names = np.append(names,ft670_names[best_ft670])
            else:                           
                data = ft670[best_ft670]
                names = np.array(ft670_names[best_ft670])
                
            print ('F670', data.shape, names.shape)
                
        if not ft970[best_ft970].size == 0:
            if not data is None:
                data = np.append(data,ft970[best_ft970])
                names = np.append(names,ft970_names[best_ft970])
            else:                           
                data = ft970[best_ft970]
                names = np.array((ft970_names[best_ft970]))
                
            print ('f970', data.shape, names.shape)
        
        if not ft1200[best_ft1200].size == 0:
            if not data is None:
                data = np.append(data,ft1200[best_ft1200])
                names = np.append(names,ft1200_names[best_ft1200])
            else:                           
                data = ft1200[best_ft1200]
                names = np.array((ft1200_names[best_ft1200]))
        
            print ('f1200', data.shape, names.shape)        
        
        print (redge)
        if np.abs(redge[0])>5 and redge[1]<0.05:
            if not data is None:
                data = np.append(data,redge[0])
                names = np.append(names,'$REIP$')
            else:                           
                data = redge[0]
                names = np.array(('$REIP$'))
            print ('r', data.shape, names.shape)
        
        if np.abs(mtbvi[2])>5 and mtbvi[3]<0.05:
            
            mtbvi_name = 'MTBVI $%s$ $%s$' %(str(int(mtbvi[0])),str(int(mtbvi[1])))
            
            if not data is None:
                data = np.append(data,mtbvi[2])
                names = np.append(names,mtbvi_name)
            else:                           
                data = mtbvi[2]
                names = np.array((mtbvi_name))
            
            print ('m', data.shape, names.shape)
        
        if data is not None:
            
            print (names)
            
            print (len(data.shape))
            
            if len(data.shape)>0:
            
                sorting = np.argsort(data)
        
                sorted_data = data[sorting]
                
                print ('NNNNNNN',names.shape, sorting.shape)
        
                sorted_names = names[sorting]
        
                return sorted_names, sorted_data
                    
        else:
            return None
        
class MultiPlot():
    
    def __init__(self):
        
        self.data ={}
    
    def append_data(self,
                    site,
                    date,
                    names,
                    data,
                    outdir):
                        
        if not site in self.data:
            self.data[site] = {}
            
        self.data[site][date]=[names,data]
        
    #def plot_multi(self):
                   
        for key in self.data:
            
            #site_name = names(key)
            
            #n_flights = len(self.data[key])
            
            #cfig = plt.subplots(n_flights)
            
            #splots = cfig[1]           
            
            #p = 0
            
            for d in sorted(self.data[key]):
                
                data = self.data[key][d][1]
                
                n = self.data[key][d][0]
                
                print (data)
                
                print (n)
                
                x = data.shape[0]
                x = np.arange(0,data.shape[0])
                
                print (x)
                
                for r in range(data.shape[0]):
                    
                    nn = n[r]  
                    
                    if nn in index_names:
                        c = index_colour
                        l = 'Vegetation indices'
                    elif nn in fluo_names:
                        c = fluo_colour
                        l = 'Fluorescence indices'
                    elif nn in ft470_names:
                        c = ft470_colour
                        l = 'Continuum removal 470nm'
                    elif nn in ft670_names:
                        c = ft670_colour
                        l = 'Continuum removal 670nm'
                    elif nn in ft970_names:
                        c = ft970_colour
                        l = 'Continuum removal 970nm'
                    elif nn in ft1200_names:
                        c = ft1200_colour
                        l = 'Continuum removal 1200nm'
                    elif 'MTBVI' in nn:
                        c = mtbvi_colour
                        l = 'MTBVI'
                    else:
                        c = redge_colour
                        l = 'Red Edge'
                        
                    
                    #splots[p].bar(r,d[r,0],color=d)
                    
                    plt.barh(r,data[r],color=c,label=l)
                
                '''plt.legend(loc=2,
                           bbox_to_anchor=(1.05,1),
                           prop={'size':10})'''
                
                plt.margins(0.01)
                plt.xticks(x+0.5, n, style ='italic',rotation='vertical',size='smaller')
                plt.xlabel('Method')
                plt.ylabel('Contrast $(T)$')
                plt.title('%s %s \n Performant methods ($T$ > 5 $P$ <0.05)' %(parse_names(key), d))
                #plt.tight_layout(rect=[0,0,0.7,0.95])
                
                lgd_ind = plt.Rectangle((0,0),1,1, fc=index_colour)
                lgd_redge = plt.Rectangle((0,0),1,1, fc= redge_colour)
                lgd_fluo = plt.Rectangle((0,0),1,1, fc=fluo_colour)
                lgd_mtbvi = plt.Rectangle((0,0),1,1, fc=mtbvi_colour)
                lgd_ft470 = plt.Rectangle((0,0),1,1, fc=ft470_colour)
                lgd_ft670 = plt.Rectangle((0,0),1,1, fc=ft670_colour)
                lgd_ft970 = plt.Rectangle((0,0),1,1, fc=ft970_colour)
                lgd_ft1200 = plt.Rectangle((0,0),1,1, fc=ft1200_colour)
                #box=ax1.get_position()
                #ax1.set_position([box.x0,box.y0,box.width*0.8, box.height])
                
                
        
                plt.legend([lgd_ind,lgd_redge,lgd_fluo,lgd_mtbvi,lgd_ft470,lgd_ft670,lgd_ft970,lgd_ft1200],
                           ['Indices','Red Edge','Fluorescence Indices','MTBVI','Continuum removal 470nm','Continuum removal 670nm','Continuum removal 970nm','Continuum removal 1200nm'],
                           loc='lower center',
                           bbox_to_anchor=(0.5,-0.9),
                           prop={'size':8},
                           ncol=4)
                
                
                plt.savefig(os.path.join(outdir,'asd_'+key+d+'.png'))
                #plt.close()     
                
                #plt.tight_layout(rect=[0,0,0.98,0.98])
                
                plt.show()
                
                plt.close()               
                
                #splots[p].set_xticks(x, minor=False)
                #splots.set_xticklabels(n)
                
                #p+=1
                
                


if __name__ == '__main__':
    
    hs_dir = '/home/dav/data/Dropbox/data_analysis/input_data/asd/surveys/plots'
    
    #mtdir = '/home/dav/data/Dropbox/hs_reprocessed/mtbvi/subsets'
    
    plotsdir = '/home/dav/data/Dropbox/data_analysis/input_data/asd/surveys/plots/performant_plots'
    
    if not os.path.exists(plotsdir):
        os.mkdir(plotsdir)

    survey_list = []
    
    for file in os.listdir(hs_dir):
        elements = file.split('_')
        extension=None
        
        print (file)
        
        if '.' in file:
        
            extension = file.split('.')[1]
        
        if extension == 'txt' and not file=='mtbvi_list.txt':
            
            s = ''
           
            if elements[6][0]=='T':
                e = 7
                
            else:
                e = 6
                
            for item in elements[:e]:
                if len(s)==0:
                    s = item
                else:
                    s=s+'_%s' %(item)
                    
            if not s in survey_list:
                survey_list.append(s)
                
    for survey in survey_list:
        
        n = survey.split('_')
        
        name = '%s_%s' %(n[2].lower(),n[3])
        
        print (name)
                
        
        indices = np.genfromtxt(os.path.join(hs_dir,survey+'_indicies.txt'), delimiter=',')
        
        redge = np.genfromtxt(os.path.join(hs_dir,survey+'_redge.txt'), delimiter=',')[2,:]
        
        fluo = np.genfromtxt(os.path.join(hs_dir,survey+'_fluo.txt'), delimiter=',')
        
        ft470 = np.genfromtxt(os.path.join(hs_dir,survey+'_470nm_featdef.txt'), delimiter=',')
        ft670 = np.genfromtxt(os.path.join(hs_dir,survey+'_670nm_featdef.txt'), delimiter=',')
        ft970 = np.genfromtxt(os.path.join(hs_dir,survey+'_970nm_featdef.txt'), delimiter=',')
        ft1200 = np.genfromtxt(os.path.join(hs_dir,survey+'_1200nm_featdef.txt'), delimiter=',')
        
        
        mtbvi = None
        
        with open(os.path.join(hs_dir,survey+'_specread_best.txt'),'r') as mt_file:
            for row in mt_file.readlines():
                mt = row.split(' ')
                if mt[0]=='MTBVI':
                    mtbvi = np.array((float(mt[1]),float(mt[2]),float(mt[4]),float(mt[6])))
                    
        plots = Plotting(name,plotsdir)
        
        #plots.appendix_plots(indices,redge,fluo,ft470,ft670,ft970,ft1200)

        p = plots.summary_plot_data(indices,redge,fluo,mtbvi,ft470,ft670,ft970,ft1200)
        
        if not p is None:

            mplots = MultiPlot()
        
            mplots_data = mplots.append_data(plots.site,plots.date,p[0],p[1],plotsdir)