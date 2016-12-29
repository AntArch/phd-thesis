1# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 12:41:31 2014

@author: david
"""

import numpy as np
import os
import csv

from scipy import stats

import matplotlib.pyplot as plt

class Transects():
    def __init__(self, tran_file, feat_file):
        self.surveys = []
        self.features = []
        with open(tran_file, 'r') as file:
            reader = csv.reader(file, delimiter=',')
            for row  in reader:
                self.surveys.append(row)
                
        with open(feat_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                self.features.append(row)
                
    def feature_pos(self, filename):
        field = filename[0:4]
        print 'Call to feature pos'
        print 'Field', field
        date = filename[4:12]
        print 'DATE:', date
        
        
        if len(filename)>16:
            t = filename[13:15]
            if t == 'T1':
                t = 'Transect 1'
            elif t == 'T2':
                t = 'Transect 2'
       
        else:
            t = 'Transect 1'
            
        print 'T',t
            
        for row in self.surveys:
            #print row[0], row[5], row[2]
            if row[0] == field and row[4] == date and row[2]==t:
                
                version = row[3]
                print 'version', version
                break
            
        feats = []
         
        for row in self.features:
            if row[0]==field and row[1]==t and row[2]==version:
                print 'Found row'
                feature = (float(row[3]),float(row[4]))
                feats.append(feature)
        print '||||||FEAT method|||||||||', feats   
        return feats
        
    def title(self, filename):
        field = filename[0:4]
        date = filename[4:12]
        
        for row in self.surveys:
            if len(filename)>16:
                t = filename[13:15]
                print t,'!!!!'
                if t == 'T1':
                    t = 'Transect 1'
                elif t == 'T2':
                    t = 'Transect 2'
           
            else:
                t = 'Transect 1'
            
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
        
class plotting():
    def __init__(self,
                 title, 
                 data1, 
                 data1_name, 
                 out_name, 
                 out_dir,
                 features,
                 data2=None, 
                 data2_name=None):
        
        os.chdir(out_dir)    
                     
        fig,ax1 = plt.subplots()
    
        ax1.set_title(title,
                      fontsize=12,
                      fontweight='bold')
        #if len(features)>1
        for feature in features:
            ax1.axvspan(feature[0]+0.5,
                        feature[1]+0.5, 
                        facecolor='Gainsboro', 
                        alpha=0.3, label='Archaelogical feature')
            
            text_anchorx = feature[0]+0.1
            
            text_anchory = data1+((np.max(data1[:,1])-np.min(data1[:,1]))/25)
            
            '''ax1.text(text_anchorx,
                     text_anchory,
                     0.83, 
                     'Archaeological \n Feature', 
                     horizontalalignment='left',
                     verticalalignment='center',
                     rotation='vertical',
                     transform=ax1.transData,
                     color='DimGray')'''
        
        #ax1.plot(ht[:,0],ht[:,1], color='Blue')
        print data1.shape
        data1_axis = [0,np.max(data1[:,0]),np.min(data1[:,1]),np.max(data1[:,1])]
        ax1.axis(data1_axis)
        
        ax1.set_xlabel('Transect Distance (m)')
        ax1.set_ylabel(data1_name, color='Blue')
        
        for t1 in ax1.get_yticklabels():
            t1.set_color('Blue')
         
        data1_binned = self.binning(data1) 
        #ax1.plot(np.arange(0.5,np.max(data1[:,0])+0.5),data1_binned[:,0], color='b', lw=2)
        #ax1.errorbar(np.arange(0.5,data1_binned.shape[0]+0.5),data1_binned[:,0], yerr=data1_binned[:,1], color='b', fmt='o') 
        
        ax1.plot(np.arange(0.5,data1_binned.shape[0]+0.5),data1_binned[:,0], color='b', lw=2)
        ax1.errorbar(np.arange(0.5,data1_binned.shape[0]+0.5),data1_binned[:,0], yerr=data1_binned[:,1], color='b', fmt='o') 
        
        if data2 != None:
            ax2 = ax1.twinx()
            #ax2.plot(lai[:,0],lai[:,1], color='Green')
            data2_axis = [0,np.max(data2[:,0]),np.min(data2[:,1]),np.max(data2[:,1])]
            ax2.axis(data2_axis)
            
            ax2.set_ylabel(data2_name, color='Green')
            
            for t1 in ax2.get_yticklabels():
                t1.set_color('Green')
                
                
            data2_binned = self.binning(data2) 
            print data2_binned.shape
            ax2.plot(np.arange(0.5,data2_binned.shape[0]+0.5),data2_binned[:,0], color='g',lw=2)
            ax2.errorbar(np.arange(0.5,data2_binned.shape[0]+0.5),data2_binned[:,0], yerr=data2_binned[:,1], color='g', fmt='o')
       
        plt.show()
        plt.savefig(out_name)
        plt.close()
            
            
    def binning(self, data):
    
        print 'binning'
        maxi = int(np.max(data[:,0]))
        print 'MAX',maxi
        
        for i in range(0,int(np.max(data[:,0])+1),1):
            print i
            if i < maxi:
                summary = data[np.where(np.logical_and(data[:,0]>i,data[:,0]<i+1)),1]
                print summary
                mean = np.mean(summary)
                std = np.std(summary)
                lower = np.min(summary)
                upper = np.max(summary)
                
                m = np.column_stack((mean,std,lower,upper))
                
                if i == 0:
                    binned = m
                else:
                    binned = np.vstack((binned, m))
        print 'binned shape',binned.shape
            
        return binned
        
if __name__ == "__main__": 
    dir_path = os.path.dirname(os.path.abspath('...'))
    
    trandir = os.path.join(dir_path,'transect')
    #print trandir
    
    datadir = os.path.join(dir_path, 'data')
    
    ht_dir = os.path.join(datadir,'ht')
    ht_list = os.listdir(ht_dir)
    #print 'HT LIST', ht_list
    
    lai_dir = os.path.join(datadir,'lai')
    lai_list = os.listdir(lai_dir)
    #print 'LAI List', lai_list
    
    den_dir = os.path.join(datadir,'den')
    den_list = os.listdir(den_dir)
    #print 'DEN List', den_list
    
    outdir = os.path.join(dir_path,'output')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    
    os.chdir(trandir)
    transect = Transects('surveys.csv','features.csv')
    
    ht_den = set(den_list).intersection(ht_list)
    
    #print '****HT DEN*****', ht_den
    
    #print transect
    
    ht_lai = set(ht_list).intersection(lai_list)
    
    #print ht_lai
    
    for survey in ht_den:
        print '*|*|*|*|*|* SURVEY *|*|*|*|*:', survey
        features = transect.feature_pos(survey)
        print 'FEAT', features
        os.chdir(ht_dir)
        ht_data = np.genfromtxt(survey, delimiter=',')
        os.chdir(den_dir)
        den_data = np.genfromtxt(survey, delimiter=',')
        plotting(transect.title(survey)+'\n Vegetation Height and Coverage',
                 ht_data,
                 'Vegetation Height (m)',
                 survey[0:-4]+'_ht_cover.png',
                 outdir,
                 features,
                 data2=den_data,
                 data2_name='Vegetation Coverage (\%)')
        
    for survey in ht_lai:
        features = transect.feature_pos(survey)
        os.chdir(ht_dir)
        ht_data = np.genfromtxt(survey, delimiter=',')
        os.chdir(lai_dir)
        lai_data = np.genfromtxt(survey, delimiter=',')
        plotting(transect.title(survey)+'\n Vegetation Height and LAI',
                 ht_data,
                 'Vegetation Height (m)',
                 survey[0:-4]+'ht_lai.png',
                 outdir,
                 features,
                 data2=lai_data,
                 data2_name='Leaf Area Index (LAI $m^2$ $m^{-2}$)')
                 
        
    ht_stats = []
    for survey in ht_list:
        
        features = transect.feature_pos(survey)
        os.chdir(ht_dir)
        ht_data = np.genfromtxt(survey, delimiter=',')
        plotting(transect.title(survey)+'\n Vegetation Height',
                 ht_data,
                 'Vegetation Height (m)',
                 survey[0:-4]+'_ht.png',
                 outdir,
                 features) 
                 
        ht_arc = []
        ht_bac = []
        for r in ht_data:
            print r.shape
            for feature in features:
                print feature
                if r[0] >= feature[0] and r[0] <=feature[1]:
                    ht_arc.append(r[1])
                else:
                    ht_bac.append(r[1])

        ht_t = stats.ttest_ind(np.asarray(ht_arc),
                               np.asarray(ht_bac),
                               equal_var=False)
              
        ht_stats.append((survey[0:-4],ht_t[0],ht_t[1]))
        
    with open('ht_stats.txt', 'w') as outfile:
        writer = csv.writer(outfile)
        for row in ht_stats:
            writer.writerow(row)
        outfile.close()       
    
    cover_stats = []
    for survey in den_list:
        features = transect.feature_pos(survey)
        os.chdir(den_dir)
        den_data = np.genfromtxt(survey, delimiter=',')
        plotting(transect.title(survey)+'\n Vegetation Coverage',
                 ht_data,
                 'Vegetation Coverage (%)',
                 survey[0:-4]+'_coverage.png',
                 outdir,
                 features) 
        cover_arc = []
        cover_bac = []
        for r in den_data:
            for feature in features:
                if r[0] >= feature[0] and r[0] <=feature[1]:
                    cover_arc.append(r[1])
                else:
                    cover_bac.append(r[1])
                    
        cover_t = stats.ttest_ind(np.asarray(cover_arc),
                                  np.asarray(cover_bac),
                                  equal_var=False)
                                  
        cover_stats.append((survey[0:-4],cover_t[0],cover_t[1]))
    
    with open('coverage_stats.txt', 'w') as outfile:
        writer = csv.writer(outfile)
        for row in cover_stats:
            writer.writerow(row)
        outfile.close()                      
        
    lai_stats = []   
    for survey in lai_list:
       
       features = transect.feature_pos(survey)
       os.chdir(lai_dir)
       lai_data = np.genfromtxt(survey, delimiter=',')
       plotting(transect.title(survey)+'\n LAI',
                ht_data,
                'Leaf Area Index (LAI $m^2$ $m^{-2}$)',
                survey[0:-4]+'_lai.png',
                outdir,
                features)   
    
       lai_arc = []
       lai_bac = []
       for r in lai_data:
           for feature in features:
               if r[0] >= feature[0] and r[0] <=feature[1]:
                   lai_arc.append(r[1])
           else:
                lai_bac.append(r[1])
            
       lai_t = stats.ttest_ind(np.asarray(lai_arc),
                            np.asarray(lai_bac),
                            equal_var=False)
                              
       lai_stats.append((survey[0:-4],lai_t[0],lai_t[1]))

    with open('lai_stats.txt', 'w') as outfile:
        writer = csv.writer(outfile)
        for row in lai_stats:
            writer.writerow(row)
        outfile.close()                      
        
        
    

    
             
    
        
                
                
