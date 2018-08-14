import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from scipy.stats import linregress
from matplotlib import rc_file

rcFile = 'rc_files/paper.rc'
exp_data = 'Exp_data/all_rate_fits.txt' 

rc_file(rcFile) 

figRows, figCols = 3,6

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

dataDict = {} 
with open(exp_data) as f : 
    for line in f.readlines() : 
        if line.startswith('#') : 
            label = line.strip()[1:] 
            if label == "WT (Ral)" : 
                break  
            data = []
        else : 
            dataPoint = [float(x) for x in line.split()]
            try : 
                data = np.vstack((data,[dataPoint])) 
            except ValueError : 
                data = dataPoint
            if line.startswith('0') : 
                data = np.flip(data,axis=0) 
                dataDict[label] = data 

plotKeys = {}
with open('Plotting_files/ColorMarkerKeys.txt') as f : 
    lines = f.readlines()
    for line in lines :
        if line.startswith('#') :
            continue
        else :
            key = line.split()[0]
            value = line.split()[1:]
            plotKeys[key] = value

print len(dataDict) 
fig, axarr = plt.subplots(figRows,figCols,figsize=[ 6.5,3 ],sharex='all')#,sharey='all') 
fig.subplots_adjust(left=0.08,right=0.95, bottom=0.15,top=0.98,hspace=0,wspace=.35) 
fig.text(0.55,0.04, r"Time (min)",ha='center', va='center') 
fig.text(0.03,0.585, r"[P$_{\rm{i}}$] ($\muup$M)", ha='center', va='center',rotation='vertical') 

for index, mut in enumerate("ADEFGHIKLMNRSTVWYQ") :
    color = plotKeys[mut][0]
    ax = axarr[index/figCols,index%figCols]
    data = dataDict[mut] 
    x     = data[:,0] 
    y1    = data[:,1] 
    yerr1 = data[:,2] 

    ax.scatter(x, y1, marker='.', color=color) 
    ax.errorbar(x, y1, yerr=yerr1 , marker='.', color=color, ls = 'None') 

    endPoint = 5 
    if mut == "S" or mut == "T" or mut == "W" : 
        endPoint = 4  ##These residues are missing the 30 min timepoint

    slope, intercept, r_value, p_value, std_error = linregress(x[:endPoint], y1[:endPoint]) 
    slope2, intercept2, r_value2, p_value2, std_error2 = linregress(x[1:endPoint], y1[1:endPoint]) 
    if std_error2 < std_error : 
        ax.text(0.15,0.83,r"*",transform=ax.transAxes,ha='center',va='top',fontsize='small',color=color) 
        slope, intercept, r_value, p_value, std_error = slope2, intercept2, r_value2, p_value2, std_error2

    ax.text(0.95,0.05,r"%.2f $\pm$ %.2f"%(slope*100, std_error * 100),transform=ax.transAxes,ha='right',va='bottom',fontsize='x-small',color=color) 
    ax.text(0.15,0.90,r"%s"%mut,transform=ax.transAxes,ha='center',va='top',fontsize='small',color=color) 

    ax.plot(x[:endPoint], x[:endPoint]*slope + intercept, ls='-', c=color,label= r"r = %0.3f"%r_value) 


fig.savefig('figures/all_fits.png',format='png') 
plt.close() 

