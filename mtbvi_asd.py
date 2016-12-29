# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 14:39:25 2016

@author: dav
"""

from matplotlib.colors import ListedColormap

import os

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm

import matplotlib.ticker as ticker


def mtbvi_plot(files,
               name,
               indir,
               outfile,):

    
    l = len(files)
    
    d= None
    
    
               
    for file in files:
        t = np.fliplr(np.genfromtxt(os.path.join(indir,file),delimiter=','))
        
        t_sig = np.abs(t > 5)
        
        
        if d is None:
            d = np.zeros(t.shape)
            
        # d[t_sig] = d[t_sig]+t[t_sig]
        
        d[t_sig] = d[t_sig]+1
        
    ld = np.max(d)+1
    
    print (ld)
    
    nt = range(int(ld))
    
    cmap = None
    
    if len(range(int(ld)))==2:
        #cmap = ListedColormap(['#636363','#bdbdbd'])
        cmap = ListedColormap(['Black','White'])
    elif len(range(int(ld)))==4:
        cmap = ListedColormap(['#525252','#969696','#cccccc','#f7f7f7'])
    elif len(range(int(ld)))==6:
        cmap = ListedColormap(['#252525','#636363','#969696','#bdbdbd','#d9d9d9','#f7f7f7'])
    
    plt.imshow(np.rot90(d,k=1), 
                        extent=[650,750,700,850],
                        interpolation='none',
                        cmap=cm.get_cmap('gnuplot2',int(ld)))
    cbar = plt.colorbar()
    
    cbar.set_label('Number of significant surveys', rotation=270, labelpad=20)
    
    locator = ticker.MaxNLocator(integer=True, nbins=int(ld))
    
    cbar.locator = locator
    
    #cbar.ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator(‌​))    
    cbar.update_ticks()

    
    
    
    '''cbar.ax.get_yaxis().set_ticks([])
    for j, lab in enumerate(nt):
        cbar.ax.text(.5, (2 * j + 1) / 8.0, lab, ha='center', va='center')
    cbar.ax.get_yaxis().labelpad = 15
    cbar.set_label('Number of significant surveys', rotation=270)'''

    
    plt.xlabel('RED $\lambda  (nm)$')
    plt.ylabel('NIR  $\lambda (nm)$')
    
    plt.title('MTBVI Performance: \n%s' %(name))
    
    #plt.show()
    
    plt.tight_layout()
    
    plt.savefig(outfile)
    
    plt.close()
    
    
            
        
        
                
    
if __name__ == '__main__':
    
    
    dir_path = '/home/dav/data/Dropbox/data_analysis/input_data/asd/surveys/plots'
    
    #dir_path = '/home/david/Dropbox/data_analysis/input_data/asd/surveys/plots'
    out_dir = os.path.join(dir_path,'mtbvi_summary')
    
    wheat = ['DART_ASD_DDT1_20110608_20110608_RAW_mtbvi_eval_t.txt',
             'DART_ASD_DDT1_20110614_20110614_RAW_mtbvi_eval_t.txt',
             'DART_ASD_DDT1_20110629_20110629_RAW_mtbvi_eval_t.txt',
             'DART_ASD_DDT1_20110715_20110715_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHQF_20120524_20120524_RAW_T1_mtbvi_eval_t.txt',
             'DART_ASD_HHQF_20120620_20120620_RAW_T1_mtbvi_eval_t.txt',
             'DART_ASD_HHQF_20120726_20120726_RAW_T1_mtbvi_eval_t.txt',
             'DART_ASD_HHQF_20120823_20120823_RAW_T1_mtbvi_eval_t.txt',
             'DART_ASD_DDCF_20120123_20120123_RAW_mtbvi_eval_t.txt',
             'DART_ASD_DDCF_20120216_20120216_RAW_mtbvi_eval_t.txt',
             'DART_ASD_DDCF_20120320_20120320_RAW_mtbvi_eval_t.txt',
             'DART_ASD_DDCF_20120516_20120516_RAW_mtbvi_eval_t.txt',
             'DART_ASD_DDCF_20120626_20120626_RAW_mtbvi_eval_t.txt',
             'DART_ASD_DDCF_20120724_20120724_RAW_mtbvi_eval_t.txt']
             
    rape = ['DART_ASD_DDT1_20111025_20111025_RAW_mtbvi_eval_t.txt',
            'DART_ASD_DDT1_20120216_20110216_RAW_mtbvi_eval_t.txt',
            'DART_ASD_DDT1_20120320_20110320_RAW_mtbvi_eval_t.txt',
            'DART_ASD_DDT1_20120522_20120522_RAW_mtbvi_eval_t.txt']
            
    '''other = ['DART_ASD_HHCC_20110929_20110929_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20111028_20111028_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20120214_20120214_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20120322_20120322_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20120524_20120524_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20120619_20120619_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20120725_20120725_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHQF_20110929_20110929_RAW_T1_mtbvi_eval_t.txt',
             'DART_ASD_HHQF_20110929_20110929_RAW_T2_mtbvi_eval_t.txt',
             'DART_ASD_HHQF_20120125_20120125_RAW_T2_mtbvi_eval_t.txt']'''
             
    other = ['DART_ASD_HHCC_20120322_20120322_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20120524_20120524_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20120619_20120619_RAW_mtbvi_eval_t.txt',
             'DART_ASD_HHCC_20120725_20120725_RAW_mtbvi_eval_t.txt']

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    mtbvi_plot(wheat,'Wheat',dir_path,os.path.join(out_dir,'wheat.png'))
    mtbvi_plot(rape,'Rape',dir_path,os.path.join(out_dir,'rape.png'))
    mtbvi_plot(other,'Other',dir_path,os.path.join(out_dir,'Pasture.png'))
        
    os.chdir(out_dir)
        
   
    
    
    