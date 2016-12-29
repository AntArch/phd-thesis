# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 15:20:05 2014

@author: david
exit

jhkj"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import os
import csv

from scipy import stats

print 'STUFF'



def edge_colours(pvals):
    # helper function to take P values for the bar plots and derive edge colour 
    # no colour = P>0.001
    # black = P<0.001
    edge_colours = []
    
    for val in pvals:
        if val < 0.001:
            edge_colours.append('#000000')
        else:
            edge_colours.append(None)
            
    return edge_colours
    
def edge_width(pvals):
    edge_widths = []
    
    for val in pvals:
        if val < 0.001:
            edge_widths.append(2)
        else:
            edge_widths.append(0)
            
    return edge_widths

class crem_plot():
    def __init__(self, 
                 arc_list, 
                 bac_list, 
                 start, 
                 end, 
                 title, 
                 fname, 
                 absfeature, 
                 wavelengths=None):
        self.arc_list = arc_list
        self.bac_list = bac_list
        self.start = start
        self.end = end
        
        self.wavelengths = wavelengths
        
        if wavelengths == None:
            self.wavelengths = np.arange(start,end+1,1)
        else:
            self.wavelengths = wavelengths
            
        print self.wavelengths
            
        arc = self.list_fixer(self.arc_list)
        self.arc_crem = arc[0]
        print 'ARC CREM',self.arc_crem.shape
        self.arc_bna = arc[1]
        
        bac = self.list_fixer(self.bac_list)
        self.bac_crem = bac[0]
        self.bac_bna = bac[1]
        
        print 'ABSFEATURESTRING*********: ', absfeature
        
        for item in arc[0]:
            print item.shape
        
        bna_plot = self.plot(title+'\nBand normalised by area: '+absfeature, 
                             self.arc_bna,
                             self.bac_bna,
                             fname+'_'+absfeature+'_bna.png')
                             
        crm_plot = self.plot(title+'\nContinuum removed spectra: '+absfeature,
                             self.arc_crem,
                             self.bac_crem,
                             fname+'_'+absfeature+'_crem.png')
                              
        bna_plot
        crm_plot
        
    def list_fixer(self, alist):
        crem_outlist = []
        bna_outlist = []
        for item in alist:
            print item.shape
            feat_start = item[0,0]
            crem = np.ones(self.wavelengths.shape[0])
            #print crem.shape
            bna = np.zeros(self.wavelengths.shape[0])
            start_index = np.argmin(np.abs(self.wavelengths-feat_start))
            end_index = start_index+item[:,0].shape[0]
            crem[start_index:end_index] = np.transpose(item[:,3])
            bna[start_index:end_index] = np.transpose(item[:,4])
            crem_outlist.append(crem)
            bna_outlist.append(bna)
        crem_out = np.asarray(crem_outlist)
        bna_out = np.asarray(bna_outlist)
        
        return crem_out, bna_out
        
    def plot(self, title, arc, bac, fname):
        print 'PLOT', arc.shape, bac.shape, arc.dtype
        arc_m = np.mean(arc, axis=0)
        print arc_m
        arc_s = np.std(arc, axis=0)
        bac_m = np.mean(bac, axis=0)
        bac_s = np.std(bac, axis=0)
        
        '''print self.wavelengths.shape, arc_m.shape
        
        plt.plot(self.wavelengths, arc_m, color='Green')
        #plt.errorbar(self.wavelengths, yerr=arc_s, color='Green', fmt='o')
        plt.plot(self.wavelenghts, bac_m)
        #plt.errorbar(self.wavelengths, yerr=bac_s, color='Green', fmt='o')
        plt.xlabel('$\lambda$ (nm)')
        plt.ylabel('Continuum Removed Reflectance')             
        plt.title(title)
        plt.savefig(fname)
        plt.close()'''
        
        fig=plt.figure()
        ax1=plt.subplot()
        
        ax1.fill_between(self.wavelengths, 
                         arc_m-(arc_s/2),
                         arc_m+(arc_s/2),
                         color = 'Green',
                         alpha=0.3,
                         label = 'Background $\sigma$')
                         
        ax1.fill_between(self.wavelengths, 
                         bac_m-(bac_s/2),
                         bac_m+(bac_s/2),
                         color = 'Blue',
                         alpha=0.3,
                         label='Archaeology $\sigma$')                      
        #plt.errorbar(wl,spec_bac_mean, yerr=spec_bac_stdev, color='Blue')
        print 'Plot again', arc_m.shape, arc_s.shape, bac_m.shape, bac_s.shape                       
        ax1.plot(self.wavelengths, arc_m, color='Green', lw=3)
        ax1.plot(self.wavelengths, bac_m, color='Blue', lw=3)
    
        ax1.axis('tight')
        plt.xlabel('$\lambda$ Wavelength (nm)')
        plt.ylabel('Continuum removed reflectance')
        plt.title(title)
        
        arc_lgnd, = ax1.plot([],[], color='Green', lw=3)
        bac_lgnd, = ax1.plot([],[], color='Blue', lw=3)
        arc_std_lgnd, = ax1.plot([],[], color='Green', lw=10, alpha=0.3)
        bac_std_lgnd, = ax1.plot([],[], color='Blue',lw=10, alpha=0.3)
        
        box = ax1.get_position()
        ax1.set_position([box.x0,box.y0,box.width*0.8, box.height])
        
        ax1.legend([arc_lgnd,bac_lgnd,arc_std_lgnd,bac_std_lgnd],
                   ['Archaeology mean','Background mean', 'Archaeology $\sigma$', 'Background $\sigma$'],
                   loc='center left',
                   bbox_to_anchor=(1,0.5),
                   prop={'size':10})
        
        #plt.tight_layout()
        #plt.show()
        
        plt.savefig(fname)
        plt.close()
            

class load_processed_spectra():
    def __init__(self, raw_path, processed_path, tran_path, plot_path):
        print 'call to processed_spectra'
        os.chdir(tran_path)
        
        
        tran = Transects('surveys.csv','features.csv')
        
        raw_list = os.listdir(raw_path)
        for folder in raw_list:
            
            features = tran.feature_pos(folder)
            print features
            
            tname = tran.title(folder)
            
            root_dir = os.path.join(raw_path,folder)
            pos_dir = os.path.join(root_dir, 'reading_info')
            
            
            
            os.chdir(pos_dir)
            
            
            reading_info_file = open('reading_atributes.txt','rb')
            
            reading_info = csv.DictReader(reading_info_file)         
            
            
            readings_list = [row for row in reading_info]
            
            spec =[]
            
            indices = []
            mtbvi = []
            redge = []
            fluo =[]
            
            cont_rem470 = []
            cont_rem670 = []
            cont_rem970 = []
            cont_rem1200 = []
            
            for reading in readings_list[:]:
                
                '''reading_filename = str(reading['reading_id']+'.txt')
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
                                                     reading['white_ref']))'''
                #print reading_info_line
                #print reading['reading_id']
                
                if reading['reading_type']== 'REF':
                    processed_readings_path = os.path.join(processed_path,folder)
                    os.chdir(processed_readings_path)
                    
                    spec_reading = np.genfromtxt(reading['reading_id']+'_spec.txt', delimiter =',')
                    #print spec_reading.shape
                    spec.append((reading['transectPosition'],spec_reading))
                    
                    indices_reading = np.genfromtxt(reading['reading_id']+'_indices.txt', skip_header=1, delimiter=',')
                    indices.append((reading['transectPosition'],indices_reading))
                    
                    fluo_reading = np.genfromtxt(reading['reading_id']+'_flou.txt', skip_header=1, delimiter=',')
                    fluo.append((reading['transectPosition'],fluo_reading))
                    
                    redge_reading = np.genfromtxt(reading['reading_id']+'_redge.txt',delimiter=',')
                    redge.append((reading['transectPosition'],redge_reading))
                    
                    mtbvi_reading = np.genfromtxt(reading['reading_id']+'_mtbvi.txt', delimiter=',')
                    mtbvi.append((reading['transectPosition'],mtbvi_reading))
                    
                    try:
                        abs470_reading = np.genfromtxt(reading['reading_id']+'_abs470_ftdef.txt', skip_header=1, delimiter=',')
                        crem470_reading = np.genfromtxt(reading['reading_id']+'_abs470_crem.txt', skip_header=1, delimiter=',')
                        cont_rem470.append((reading['transectPosition'],abs470_reading, crem470_reading))
                        
                    except:
                        continue
                    
                    
                    try:
                        abs670_reading = np.genfromtxt(reading['reading_id']+'_abs670_ftdef.txt', skip_header=1, delimiter=',')
                        crem670_reading = np.genfromtxt(reading['reading_id']+'_abs670_crem.txt', skip_header=1, delimiter=',')
                        cont_rem670.append((reading['transectPosition'],abs670_reading, crem670_reading))
                        
                    except:
                        continue
                    
                    try:
                        abs970_reading = np.genfromtxt(reading['reading_id']+'_abs970_ftdef.txt', skip_header=1, delimiter=',')
                        crem970_reading = np.genfromtxt(reading['reading_id']+'_abs970_crem.txt', skip_header=1, delimiter=',')
                        cont_rem970.append((reading['transectPosition'],abs970_reading, crem970_reading))
                        
                    except:
                        continue
                    
                    try:
                        abs1200_reading = np.genfromtxt(reading['reading_id']+'_abs1200_ftdef.txt', skip_header=1, delimiter=',')
                        crem1200_reading = np.genfromtxt(reading['reading_id']+'_abs1200_crem.txt', skip_header=1, delimiter=',')
                        cont_rem1200.append((reading['transectPosition'],abs1200_reading, crem1200_reading))
                        
                    except:
                        continue
            
            os.chdir(plot_path)
            
            spec_arc =[]
            spec_bac =[]
            for entry in spec:
                #print entry
                arc_flag = False
                tpos = float(entry[0])
                data = entry[1][:,1]
                wl = np.transpose(entry[1][:,0])
                print 'WL',wl.shape
                for feat in features:
                    print 'FEAT',feat[0], feat[1], '|||', tpos
                    print type(feat[0]), type(tpos)
                    if tpos > feat[0] and tpos < feat[1]:
                        print 'arc'
                        arc_flag = True
                        spec_arc.append(data)
                if arc_flag == False:
                    print 'bac'
                    spec_bac.append(data)
            
            
            print len(spec_arc)
            for item in spec_arc:
                print item.shape
                if item.shape[0]==2152:
                    spec_arc.remove(item)
            
            
            spec_arc = np.transpose(np.asarray(spec_arc))
            print spec_arc.shape
            
            print len(spec_bac)
            for item in spec_bac:
                print item.shape
                if item.shape[0]==2152:
                    spec_bac.remove(item)
            spec_bac = np.transpose(np.asarray(spec_bac))
            
            print 'BAC', spec_bac.shape
            

            spec_arc_mean = np.mean(spec_arc, axis=1)
            print spec_arc_mean.shape
            spec_arc_stdev =np.std(spec_arc, axis=1)
            
            spec_bac_mean = np.mean(spec_bac, axis=1)
            
            spec_bac_stdev =np.std(spec_bac, axis=1)
            
            #print wl.shape, spec_arc_mean.shape
            
            nmask = np.zeros(wl.shape[0])
            
            
            
            nmask[0:990]=1
            nmask[1080:1438]=1
            nmask[1622:2100]=1
            
            print nmask
            
            spec_arc_mean = ma.array(spec_arc_mean, mask=nmask==0)
            spec_bac_mean = ma.array(spec_bac_mean, mask=nmask==0)
            
            spec_arc_stdev = ma.array(spec_arc_stdev, mask=nmask==0)
            spec_bac_stdev = ma.array(spec_bac_stdev, mask=nmask==0)
            
            #plt.errorbar(wl,spec_arc_mean, yerr=spec_arc_stdev, color='Green')
            fig=plt.figure()
            ax1=plt.subplot()
            
            ax1.fill_between(wl, 
                             spec_bac_mean-(spec_bac_stdev/2),
                             spec_bac_mean+(spec_bac_stdev/2),
                             color = 'Blue',
                             alpha=0.3,
                             label = 'Background $\sigma$')
                             
            ax1.fill_between(wl, 
                             spec_arc_mean-(spec_arc_stdev/2),
                             spec_arc_mean+(spec_arc_stdev/2),
                             color = 'Green',
                             alpha=0.3,
                             label='Archaeology $\sigma$')                         
            #plt.errorbar(wl,spec_bac_mean, yerr=spec_bac_stdev, color='Blue')
            ax1.plot(wl, spec_arc_mean, color='Green', lw=3)
            ax1.plot(wl, spec_bac_mean, color='Blue', lw=3)
        
            ax1.axis('tight')
            plt.xlabel('$\lambda$ Wavelength (nm)')
            plt.ylabel('Reflectance (%)')
            plt.title(tname+': Mean reflectance')
            
            box = ax1.get_position()
            ax1.set_position([box.x0,box.y0,box.width*0.8, box.height])
            
            arc_lgnd, = ax1.plot([],[], color='Green', lw=3)
            bac_lgnd, = ax1.plot([],[], color='Blue', lw=3)
            arc_std_lgnd, = ax1.plot([],[], color='Green', lw=10, alpha=0.3)
            bac_std_lgnd, = ax1.plot([],[], color='Blue',lw=10, alpha=0.3)
            
            ax1.legend([arc_lgnd,bac_lgnd,arc_std_lgnd,bac_std_lgnd],
                       ['Archaeology mean','Background mean', 'Archaeology $\sigma$', 'Background $\sigma$'],
                       loc='center left',
                       bbox_to_anchor=(1,0.5),
                       prop={'size':10})
            
            plt.savefig(folder+'_spec.png')
            plt.close()
            
            '''spec_t = stats.ttest_ind(spec_arc, spec_bac, equal_var=False)
            
            plt.plot(wl, spec_t)
            plt.xlabel('$\lambda$ Wavelength (nm)')
            plt.ylabel('Contrast (T score)')
            plt.axis('tight')
            plt.title(tname+'Contrast by wavelength')'''
            
        
            indices_arc =[]
            indices_bac =[]
            for entry in indices:
                arc_flag = False
                tpos = float(entry[0])
                data = entry[1]
                for feat in features:
                    if tpos > feat[0] and tpos < feat[1]:
                        arc_flag = True
                        indices_arc.append(data)
                if arc_flag == False:
                    indices_bac.append(data)
            indices_arc=np.asarray(indices_arc)[0:-1]
            indices_bac=np.asarray(indices_bac)[0:-1]
            
            print 'INDICES SHAPES', indices_arc.shape,indices_bac.shape
            
            
            
            index_names = ('SR 700 800',
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
                           'NDLI')
            
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
                       '#f1b6da',
                       '#f1b6da',
                       '#f1b6da',
                       '#f1b6da',
                       '#f1b6da',
                       '#f1b6da']
            
            indices_t = stats.ttest_ind(indices_arc, indices_bac, equal_var=False)
            print indices_t[0].shape
            xes = np.arange(0,indices_t[0].shape[0],1)
            
            best_index_idx = np.nanargmax(indices_t[0])
            best_index_name = index_names[best_index_idx]
            best_index_t = indices_t[0][best_index_idx]
            best_index_p = indices_t[1][best_index_idx]
            
            np.savetxt(folder+'_indicies.txt', 
                       np.column_stack((indices_t[0],indices_t[1])), 
                       delimiter=',')
            
            best_index_string = 'Best Index: %s & %s & %s \\\\' %(best_index_name, str(np.around(best_index_t, decimals=2)), str(best_index_p))
            
            
            fig=plt.figure()
            ax1 = plt.subplot()
            ax1.bar(xes, np.abs(indices_t[0]), color=colours, linewidth=edge_width(indices_t[1]))
            ax1.margins(0.01)
            plt.xticks(xes+0.5, index_names, style ='italic', rotation='vertical', size='smaller')
            plt.xlabel('Vegetation Indices')
            plt.ylabel('Contrast (T value)')
            plt.title(tname+': Vegetation Indices')
            plt.tight_layout()
            
            '''lgd_veg = plt.plot([],[], color='#4dac26', lw=10)
            lgd_redge = plt.plot([],[], color= '#d01c8b', lw=10)
            lgd_pig = plt.plot([],[], color='#b8e186',lw=10)
            lgd_chem = plt.plot([],[], color='#f1b6da',lw=10)'''
            
            lgd_veg = plt.Rectangle((0,0),1,1, fc='#4dac26', linewidth=0)
            lgd_redge = plt.Rectangle((0,0),1,1, fc= '#d01c8b', linewidth=0)
            lgd_pig = plt.Rectangle((0,0),1,1, fc='#b8e186', linewidth=0)
            lgd_chem = plt.Rectangle((0,0),1,1, fc='#f1b6da', linewidth=0)
            lgd_sig = plt.Rectangle((0,0),1,1, fc='#FFFFFF', linewidth=2)
            
            box=ax1.get_position()
            ax1.set_position([box.x0,box.y0,box.width*0.8, box.height])
            
            plt.legend([lgd_veg,lgd_redge,lgd_pig,lgd_chem,lgd_sig],
                       ['Generic','Red Edge','Pigments','Biochemicals','P<0.001'],
                       loc='center left',
                       bbox_to_anchor=(1,0.5),
                       prop={'size':10})
            
            plt.savefig(folder+'_indices.png')
            
            plt.close()
            
                    
            fluo_arc =[]
            fluo_bac =[]
            for entry in fluo:
                arc_flag = False
                tpos = float(entry[0])
                data = entry[1]
                for feat in features:
                    if tpos > feat[0] and tpos < feat[1]:
                        arc_flag = True
                        fluo_arc.append(data)
                if arc_flag == False:
                    fluo_bac.append(data)
            fluo_arc = np.asarray(fluo_arc)
            fluo_bac = np.asarray(fluo_bac)
            
            fluo_names = ('$680/630$',
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
                          "$max f'/ max f'+12$")
            
            
            fluo_t = stats.ttest_ind(fluo_arc, fluo_bac, equal_var=False)
            xes = np.arange(0,fluo_arc.shape[1],1)
            plt.bar(xes, np.abs(fluo_t[0]), color='#4dac26', linewidth=edge_width(fluo_t[1]))
            plt.margins(0.01)
            plt.xticks(xes+0.5, fluo_names, style ='italic',rotation='vertical',size='smaller')
            plt.xlabel('Fluorescence Indices')
            plt.ylabel('Contrast (T value)')
            plt.title(tname+': Fluorescence Indices')
            plt.tight_layout()
            plt.savefig(folder+'_fluo_t.png')
            plt.close()
            #plt.show() 
            
            best_fluo_t_idx = np.nanargmax(np.abs(fluo_t[0]))
            best_fluo_t = fluo_t[0][best_fluo_t_idx]
            best_fluo_p = fluo_t[1][best_fluo_t_idx]
            best_fluo_name = fluo_names[best_fluo_t_idx]
            
            
            np.savetxt(folder+'_fluo.txt', 
                       np.column_stack((fluo_t[0],fluo_t[1])), 
                       delimiter=',')
            
            best_fluo_string = 'Fluorescence %s & %s & %s \\\\' %(str(best_fluo_name), str(np.around(best_fluo_t, decimals=2)), str(best_fluo_p))
             
            redge_arc = []
            redge_bac = []
            for entry in redge:
                arc_flag = False
                tpos = float(entry[0])
                data = entry[1]
                for feat in features:
                    if tpos > feat[0] and tpos < feat[1]:
                        arc_flag = True
                        redge_arc.append(data)
                if arc_flag == False:
                    redge_bac.append(data)
                    
            redge_tt = stats.ttest_ind(np.asarray(redge_arc),np.asarray(redge_bac), equal_var=False)
            
            print len(redge_arc), redge_arc
            redge_mean = np.mean(np.asarray(redge_arc)[:,2])
            
            np.savetxt(folder+'_redge.txt', 
                       np.column_stack((redge_tt[0],redge_tt[1])), 
                       delimiter=',')
            
            best_redge_string = 'REIP %s & %s & %s \\\\' %(str(np.around(redge_mean, decimals=1)), str(np.around(redge_tt[0][2], decimals=2)),str(redge_tt[1][2]))
                    
            mtbvi_arc= []
            mtbvi_bac = []
            mtbvi_hdr_flag = False
            for entry in mtbvi:
                arc_flag = False
                tpos=float(entry[0])
                data=np.transpose(entry[1])
                if mtbvi_hdr_flag == False:
                    hdr = data[0:2,:]
                    mtbvi_hdr_flag=True
                data = data[2,:]
                for feat in features:
                    if tpos > feat[0] and tpos < feat[1]:
                        arc_flag = True
                        mtbvi_arc.append(data)
                if arc_flag == False:
                    mtbvi_bac.append(data)                
            mtbvi_arc = np.asarray(mtbvi_arc)
            mtbvi_bac = np.asarray(mtbvi_bac)
            mtbvi_tt = stats.ttest_ind(mtbvi_arc, mtbvi_bac, equal_var=False)
            hdr = np.asarray(hdr)   
            x = np.unique(hdr[0,:])
            y = np.unique(hdr[1,:])
            
            xno = x.shape[0]
            yno = y.shape[0] 
       
            print 'MTBVI T', mtbvi_tt[0].shape, mtbvi_tt[1].shape
            
            mtbvi_t = np.abs(mtbvi_tt[0])
            mtbvi_p = np.abs(mtbvi_tt[1])
            
            mtbvi_matrix = np.flipud(np.transpose(np.reshape(mtbvi_t, (xno,yno))))
            mtbvi_pmatrix = np.flipud(np.transpose(np.reshape(mtbvi_p, (xno,yno))))
            
            np.savetxt(folder+'_mtbvi_eval_t.txt', mtbvi_matrix,delimiter=',')
            np.savetxt(folder+'_mtbvi_eval_p.txt', mtbvi_pmatrix,delimiter=',')
            
            mtbvi_eval = np.nonzero(np.logical_and(mtbvi_matrix>5,mtbvi_pmatrix<0.05))
            
            '''np.savetxt(folder+'_mtbvi_eval.txt', 
                       np.column_stack((mtbvi_eval[0],mtbvi_eval[1])), 
                       delimiter=',')'''
        
            #mtbvi_matrix = np.transpose(np.reshape(mtbvi_t, (xno,yno)))
            #mtbvi_matrix = np.flipud(np.transpose(np.reshape(mtbvi_t, (xno,yno)))) 
            print mtbvi_matrix
            print mtbvi_matrix.shape
            #plt.pcolor(mtbvi_matrix, vmin=np.min(mtbvi_t), vmax=np.max(mtbvi_t))
            plt.imshow(mtbvi_matrix, 
                       extent=[x[0],x[-1],y[0],y[-1]],
                       interpolation='none',
                       cmap='gnuplot2')
            plt.colorbar()
            plt.xlabel('RED $\lambda  (nm)$')
            plt.ylabel('NIR  $\lambda (nm)$')
            plt.title(tname+':\n NDVI variable wavelength: Contrast')
            plt.tight_layout()
            
            plt.savefig(folder+'_mtbvi')
            plt.close()
            
                               
            best_ndvi_idx = np.nanargmax(np.abs(mtbvi_t))
            print best_ndvi_idx, mtbvi_t.shape, mtbvi_t[best_ndvi_idx]
            best_mtbvi_tval = mtbvi_t[best_ndvi_idx]
            print mtbvi_p.shape 
            best_mtbvi_pval = mtbvi_p[best_ndvi_idx]
            print hdr.shape
            #?
            best_mtbvi_name = hdr[:,best_ndvi_idx]
            
            print 'MTBVI***************************NAME', best_mtbvi_name, best_mtbvi_name.shape
            #best_mtbvi_name = hdr[best_ndvi_idx,:]
            
            
            best_mtbvi_string = 'MTBVI %s %s & %s & %s \\\\' %(str(best_mtbvi_name[0]),str(best_mtbvi_name[1]),str(np.around(best_mtbvi_tval, decimals=2)),str(best_mtbvi_pval))
            
            

            cont_rem470_abs_arc = []
            cont_rem470_crm_arc = []
            cont_rem470_abs_bac = []
            cont_rem470_crm_bac = []
            for entry in cont_rem470:
                arc_flag = False
                tpos = float(entry[0])
                abs470 = entry[1]
                crm470=entry[2]
                for feat in features:
                    if tpos > feat[0] and tpos < feat[1]:
                        arc_flag = True
                        cont_rem470_abs_arc.append(abs470)
                        cont_rem470_crm_arc.append(crm470)
                if arc_flag == False:
                    cont_rem470_abs_bac.append(abs470)
                    cont_rem470_crm_bac.append(crm470)
            cr470 = crem_plot(cont_rem470_crm_arc,
                    cont_rem470_crm_bac,
                    400,
                    518,
                    tname,
                    folder,
                    '470nm')
                    
            cr470_t = stats.ttest_ind(cr470.arc_bna, cr470.bac_bna, equal_var=False)
            print cr470_t[0].shape
            
            '''cr470ratio = CRRatio(cr470.arc_bna, cr470.bac_bna)
            
            cr470ratio_t = cr470ratio.ratio()''' 
            
            try:
                cr470_best_idx = np.nanargmax(np.abs(cr470_t[0]))
                print 'CR IDX',cr470_best_idx
            except:
                cr470_best_idx = np.argmin(np.abs(cr470.wavelengths-484))
                
            cr470_best_wl = cr470.wavelengths[cr470_best_idx]
            cr470_best_t = cr470_t[0][cr470_best_idx]
            cr470_best_p = cr470_t[1][cr470_best_idx]
 
            cr470_best_string = '470nm CR BNA Best $\\lambda$ %s & %s & %s \\\\' %(str(np.around(cr470_best_wl)), str(np.around(cr470_best_t, decimals=2)), str(cr470_best_p))
            '''cr470_ratio_string = '470nm CR BNA Ratio $\\lambda_{%s}/lambda_{%s}$  & %s & %s \\\\' %(str(cr470.wavelengths[cr470ratio.min_idx]),
                                                                                                    str(cr470.wavelengths[cr470ratio.max_idx]), 
                                                                                                    str(cr470ratio_t[0]), str(cr470ratio_t[1]))'''
            
            ft470 = ft_def_plot(cont_rem470_abs_arc,
                                cont_rem470_abs_bac,
                                tname,
                                folder,
                                '470nm')
            
            cont_rem670_abs_arc = []
            cont_rem670_crm_arc = []
            cont_rem670_abs_bac = []
            cont_rem670_crm_bac = []
            for entry in cont_rem670:
                arc_flag = False
                tpos = float(entry[0])
                abs670 = entry[1]
                crm670=entry[2]
                for feat in features:
                    if tpos > feat[0] and tpos < feat[1]:
                        arc_flag = True
                        cont_rem670_abs_arc.append(abs670)
                        cont_rem670_crm_arc.append(crm670)
                if arc_flag == False:
                    cont_rem670_abs_bac.append(abs670)
                    cont_rem670_crm_bac.append(crm670)
            cr670 = crem_plot(cont_rem670_crm_arc,
                    cont_rem670_crm_bac,
                    548,
                    800,
                    tname,
                    folder,
                    '670nm')
            ft670 = ft_def_plot(cont_rem670_abs_arc,
                    cont_rem670_abs_bac,
                    tname,
                    folder,
                    '670nm')
             
            
            cr670_t = stats.ttest_ind(cr670.arc_bna, cr670.bac_bna, equal_var=False)
            print 'CR670', cr670_t
            
            
            try:
                cr670_best_idx = np.nanargmax(np.abs(cr670_t[0]))
            except:
                cr670_best_idx = np.argmin(np.abs(cr670.wavelengths-670))
            print 'CR670 BEST WL', cr670_best_idx
            cr670_best_wl = cr670.wavelengths[cr670_best_idx]
            print cr670_best_wl
            cr670_best_t = cr670_t[0][cr670_best_idx]
            cr670_best_p = cr670_t[1][cr670_best_idx]
            print 'cr670 t,p', cr670_best_t, cr670_best_p
            
            '''cr670ratio = CRRatio(cr670.arc_bna, cr670.bac_bna)
            
            cr670ratio_t = cr670ratio.ratio()''' 
 
            cr670_best_string = '670nm CR BNA Best $\\lambda$ %s & %s & %s \\\\' %(str(cr670_best_wl), str(np.around(cr670_best_t, decimals=2)), str(cr670_best_p))
            '''cr670_ratio_string = '670nm CR BNA Ratio $\\lambda_{%s}/lambda_{%s}$  & %s & %s \\\\' %(str(cr670.wavelengths[cr670ratio.min_idx]),
                                                                                                    str(cr670.wavelengths[cr670ratio.max_idx]), 
                                                                                                    str(cr670ratio_t[0]), str(cr670ratio_t[1]))'''
            
            
            cont_rem970_abs_arc = []
            cont_rem970_crm_arc = []
            cont_rem970_abs_bac = []
            cont_rem970_crm_bac = []
            for entry in cont_rem970:
                arc_flag = False
                tpos = float(entry[0])
                abs970 = entry[1]
                crm970=entry[2]
                for feat in features:
                    if tpos > feat[0] and tpos < feat[1]:
                        arc_flag = True
                        cont_rem970_abs_arc.append(abs970)
                        cont_rem970_crm_arc.append(crm970)
                if arc_flag == False:
                    cont_rem970_abs_bac.append(abs970)
                    cont_rem970_crm_bac.append(crm970)
            cr970 = crem_plot(cont_rem970_crm_arc,
                    cont_rem970_crm_bac,
                    880,
                    1115,
                    tname,
                    folder,
                    '970nm')
            ft970 = ft_def_plot(cont_rem970_abs_arc,
                    cont_rem970_abs_bac,
                    tname,
                    folder,
                    '970nm')
                        
                        
            cr970_t = stats.ttest_ind(cr970.arc_bna, cr970.bac_bna, equal_var=False)
            
            try:
                cr970_best_idx = np.nanargmax(np.abs(cr970_t[0]))
            except:
                cr970_best_idx = np.argmin(np.abs(cr970.wavelengths-970))
            cr970_best_wl = cr970.wavelengths[cr970_best_idx]
            cr970_best_t = cr970_t[0][cr970_best_idx]
            cr970_best_p = cr970_t[1][cr970_best_idx]
 
            '''cr970ratio = CRRatio(cr970.arc_bna, cr970.bac_bna)
            
            cr970ratio_t = cr970ratio.ratio()''' 
            
            cr970_best_string = '970nm CR BNA Best $\\lambda$ %s & %s & %s \\\\' %(str(cr970_best_wl), str(np.around(cr970_best_t, decimals=2)), str(cr970_best_p)) 
            '''cr970_ratio_string = '970nm CR BNA Ratio $\\lambda_{%s}/lambda_{%s}$  & %s & %s \\\\' %(str(cr970.wavelengths[cr970ratio.min_idx]),
                                                                                                    str(cr970.wavelengths[cr970ratio.max_idx]), 
                                                                                                    str(cr970ratio_t[0]), str(cr970ratio_t[1]))'''

            cont_rem1200_abs_arc =[]
            cont_rem1200_crm_arc = []
            cont_rem1200_abs_bac =[]
            cont_rem1200_crm_bac = []
            print '1200',tname, len(cont_rem1200)
            for entry in cont_rem1200:
                arc_flag = False
                tpos = float(entry[0])
                abs1200 = entry[1]
                crm1200=entry[2]
                print tpos
                for feat in features:
                    if tpos > feat[0] and tpos < feat[1]:
                        arc_flag = True
                        cont_rem1200_abs_arc.append(abs1200)
                        cont_rem1200_crm_arc.append(crm1200)
                if arc_flag == False:
                    cont_rem1200_abs_bac.append(abs1200)
                    cont_rem1200_crm_bac.append(crm1200)
                    
            print 'CONSTRUCTOR 1220 shapes', len(cont_rem1200_crm_arc), len(cont_rem1200_crm_bac)
            cr1200 = crem_plot(cont_rem1200_crm_arc,
                     cont_rem1200_crm_bac,
                     1080,
                     1300,
                     tname,
                     folder,
                     '1200nm')
            ft1200 = ft_def_plot(cont_rem1200_abs_arc,
                     cont_rem1200_abs_bac,
                     tname,
                     folder,
                     '1200nm')
            
            cr1200_t = stats.ttest_ind(cr1200.arc_bna, cr1200.bac_bna, equal_var=False)
            try:
                cr1200_best_idx = np.nanargmax(np.abs(cr1200_t[0]))
            except:
                cr1200_best_idx = np.argmin(np.abs(cr1200.wavelengths-1200))
            cr1200_best_wl = cr1200.wavelengths[cr1200_best_idx]
            cr1200_best_t = cr1200_t[0][cr1200_best_idx]
            cr1200_best_p = cr1200_t[1][cr1200_best_idx]
 
            #cr1200ratio = CRRatio(cr1200.arc_bna, cr1200.bac_bna)
            
            #cr1200ratio_t = cr1200ratio.ratio()  
             
            cr1200_best_string = '1200nm CR BNA Best $\\lambda$ %s & %s & %s \\\\' %(str(cr1200_best_wl), str(np.around(cr1200_best_t, decimals=2)), str(cr1200_best_p))
            #cr1200_ratio_string = '1200nm CR BNA Ratio $\\lambda_{%s}/lambda_{%s}$  & %s & %s \\\\' %(str(cr1200.wavelengths[cr1200ratio.min_idx]),
            #                                                                                       str(cr1200.wavelengths[cr1200ratio.max_idx]), 
            #                                                                                       str(cr1200ratio_t[0]), str(cr1200ratio_t[1]))
            
            
            '''results_table = str(best_index_string+'\n'+
                                best_mtbvi_string+'\n'+ 
                                best_redge_string+'\n'+
                                best_fluo_string+'\n'+
                                cr470_best_string+'\n'+
                                cr470_ratio_string+'\n'+
                                cr670_best_string+'\n'+
                                cr670_ratio_string+'\n'+
                                cr970_best_string+'\n'+
                                cr970_ratio_string+'\n'+
                                cr1200_best_string+'\n'+
                                cr1200_ratio_string+'\n')'''
                                
            results_table = str(best_index_string+'\n'+
                                best_mtbvi_string+'\n'+ 
                                best_redge_string+'\n'+
                                best_fluo_string+'\n'+
                                cr470_best_string+'\n'+
                                cr670_best_string+'\n'+
                                cr970_best_string+'\n'+
                                cr1200_best_string)
                                
                                
            
            write_table = open(folder+'_specread_best.txt', 'w')
            write_table.writelines(results_table)
            write_table.close()
            
            '''multi_graph = np.column_stack((indices_t[0][0],
                                           best_mtbvi_tval,
                                           indices_t[0][4],
                                           indices_t[0][9],
                                           indices_t[0][17],
                                           indices_t[0][18],
                                           indices_t[0][28],
                                           redge_tt[0][2],
                                           best_fluo_t,
                                           ft470.ftdef_t[0][8],
                                           cr470_best_t,
                                           cr470ratio_t[0],
                                           ft670.ftdef_t[0][8],
                                           cr670_best_t,
                                           cr670ratio_t[0],
                                           ft970.ftdef_t[0][8],
                                           cr970_best_t,
                                           cr970ratio_t[0],
                                           ft1200.ftdef_t[0][8],
                                           cr1200_best_t,
                                           cr1200ratio_t[0]))'''
                                           
            multi_graph = np.column_stack((indices_t[0][0],
                                           best_mtbvi_tval,
                                           indices_t[0][4],
                                           indices_t[0][9],
                                           indices_t[0][17],
                                           indices_t[0][18],
                                           indices_t[0][28],
                                           redge_tt[0][2],
                                           best_fluo_t,
                                           ft470.ftdef_t[0][8],
                                           cr470_best_t,
                                           ft670.ftdef_t[0][8],
                                           cr670_best_t,
                                           ft970.ftdef_t[0][8],
                                           cr970_best_t,
                                           ft1200.ftdef_t[0][8],
                                           cr1200_best_t))                                           
                                           
            np.savetxt(folder+'_multi_graph.txt', 
                       multi_graph, delimiter=',')
            
            
            
class ft_def_plot():
    def __init__(self,
                 ftdef_arc,
                 ftdef_bac,
                 tname,
                 folder,
                 feat_name):
        arc_array = np.asarray(ftdef_arc)
        bac_array   = np.asarray(ftdef_bac)
        
        self.ftdef_t = stats.ttest_ind(arc_array, bac_array, equal_var=False)
        
        names = ["Refined start",
                 "Refined end",
                 "Minima Wavelength",
                 "Minima Reflectance",
                 "Feature Area",
                 "Left TBVI",
                 "Right TBVI",
                 "Continuum Gradient",
                 "CR Area",
                 "CR Maxima",
                 "CR Maxima WL",
                 "CR Area Left",
                 "CR Area Right"]
        
        xes = np.arange(0, self.ftdef_t[0].shape[0], 1)
        
        fig=plt.figure()
        ax1 = plt.subplot()
        ax1.bar(xes, np.abs(self.ftdef_t[0]), color='#4dac26', linewidth=edge_width(self.ftdef_t[1]))
        ax1.margins(0.01)
        plt.xticks(xes+0.5, names, style ='italic', rotation='vertical', size='smaller')
        plt.xlabel('Metric')
        plt.ylabel('Contrast (T)')
        plt.title(tname+'\n'+'Continuum removal metrics: '+feat_name)
        plt.tight_layout()
        plt.savefig(folder+'_'+feat_name+'_featdef.png')
        plt.close()
        
        np.savetxt(folder+'_'+feat_name+'_featdef.txt', 
                   np.column_stack((self.ftdef_t[0],self.ftdef_t[1])), 
                   delimiter=',')
        
        
class CRRatio():
    def __init__(self, bna_arc, bna_bac):   
        self.bna_arc = bna_arc
        self.bna_bac = bna_bac
        ttest = stats.ttest_ind(bna_arc, bna_bac, equal_var=False)
        print ttest
        self.max_idx = np.nanargmax(ttest[0])
        self.min_idx = np.nanargmin(ttest[0])
        
    def ratio(self):
        print self.min_idx, self.max_idx
        arc_ratio = self.bna_arc[:,self.min_idx]/self.bna_arc[:,self.max_idx]
        bac_ratio = self.bna_bac[:,self.min_idx]/self.bna_bac[:,self.max_idx]
        
        print 'ARC CRRR HERE THIS NOW', arc_ratio
        
        ratio_t = stats.ttest_ind(arc_ratio,bac_ratio, equal_var=False)
        
        return ratio_t
        
        

class Transects():
    def __init__(self, tran_file, feat_file):
        print 'call to trnasects'
        self.surveys = []
        self.features = []
        with open(tran_file, 'r') as file:
            reader = csv.reader(file, delimiter=',')
            for row  in reader:
                self.surveys.append(row)
                #print row
                
        with open(feat_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                self.features.append(row)
                
    def feature_pos(self, filename):
        filename = filename[9:]
        print filename
        field = filename[0:4]
        print 'Call to feature pos'
        print 'Field', field
        date = filename[5:13]
        print 'DATE:', date
        
        print filename[-2:]
        if filename[-2:]=='T1':
            t = 'Transect 1'
        
        elif filename[-2:]=='T2':
            t = 'Transect 2'
        
        else:
            t = 'Transect 1'
            
        print 'T',t
            
        for row in self.surveys:
            #print 'ROW', row[0], row[4], row[2]
            #print row[0], row[5], row[2]
            if row[0] == field and row[4] == date and row[2]==t:
                #print row
                if t == 'Transect 1':
                    transect = 'T1'
                else:
                    transect = 'T2'
                    
                version = row[3]
                #print 'version', version
                break
            
        feats = []
         
        for row in self.features:
            #print row
            if row[0]==field and row[1]==transect and row[2]==version:
                #print 'Found row'
                feature = (float(row[3]),float(row[4]))
                feats.append(feature)
        print '||||||FEAT method|||||||||', feats   
        return feats
        
    def title(self, filename):
        filename = filename[9:]
        field = filename[0:4]
        print 'Call to feature pos'
        print 'Field', field
        date = filename[5:13]
        print 'DATE:', date
        
        if filename[-2:]=='T1':
            t = 'Transect 1'
        
        elif filename[-2:]=='T2':
            t = 'Transect 2'
        
        else:
            t = 'Transect 1'
            
        print 'T',t
         
        for row in self.surveys:
            string=None
            
            if row[0] == field and row[4] == date and row[2]==t:
                print 'Titile'
                
                fieldname = row[1]
                datestring = str(date[0:4]+'/'+date[4:6]+'/'+date[6:])
                print 'date^^^^^',date,datestring
                #transect = row[3]
                if row[0]=='HHQF':
                    string = str(fieldname+' '+datestring+' '+t)
                else:
                    string = str(fieldname+' '+datestring)
                
                
                break
            print string
            
        return string

                    

if __name__ == "__main__":
    #dir_path = os.path.dirname(os.path.abspath('...'))

    
    dir_path = '/home/david/Dropbox/data_analysis/input_data/asd/surveys'
    
    print 'dir_path', dir_path
    raw_path = os.path.join(dir_path, 'processed')
    print 'raw', raw_path
    processed_path = os.path.join(dir_path,'output_ascii')
    
    mask_dir = os.path.join(dir_path, 'noise')
    os.chdir(mask_dir)
    mask = np.transpose(np.genfromtxt('Mask.csv',delimiter=','))
    
    trandir = os.path.join(dir_path, 'ground_plots')
    
    plot_path = os.path.join(dir_path, 'plots')
    if not os.path.exists(plot_path):
        os.mkdir(plot_path)
        
    blah = load_processed_spectra(raw_path, processed_path, trandir, plot_path)
    
    blah 
    
