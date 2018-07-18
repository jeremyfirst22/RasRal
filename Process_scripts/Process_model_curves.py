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
exp_data = 'Exp_data/modelCurvesFigure.txt' 

rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

headlines = 0 
with open(exp_data) as f : 
    for line in f.readlines() : 
        if line.startswith('#') : 
            headlines +=1 
        else : 
            break 
data = np.genfromtxt(exp_data, skip_header=headlines)

x     = data[:,0] 
y1    = data[:,1] 
yerr1 = data[:,2] 
y2    = data[:,3] 
yerr2 = data[:,4] 
y3    = data[:,5] 
yerr3 = data[:,6] 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.12,right=0.95, bottom=0.15) 
fig.text(0.5,0.04, r"Time (min)", ha='center', va='center') 
fig.text(0.03,0.5, r"[ $\rm{P_i}$ ] ($\rm{\mu M})$", ha='center', va='center',rotation='vertical') 

ax.scatter(x, y1, marker='o', color='k',label=r"WT Ras") 
ax.scatter(x, y2, marker='o', color='g',label=r"Ras Q61G") 
ax.scatter(x, y3, marker='o', color='b',label=r"Ras Q61H") 

ax.errorbar(x, y1, yerr=yerr1 , marker='o', color='k' ,capsize=3, ls = 'None') 
ax.errorbar(x, y2, yerr=yerr2 , marker='o', color='g' ,capsize=3, ls = 'None') 
ax.errorbar(x, y3, yerr=yerr3 , marker='o', color='b' ,capsize=3, ls = 'None') 

slope, intercept, r_value, p_value, std_error = linregress(x[-5:], y1[-5:],) 
ax.plot(x[-5:], x[-5:]*slope + intercept, ls='--', c='k') 

slope, intercept, r_value, p_value, std_error = linregress(x[-5:], y2[-5:],) 
ax.plot(x[-5:], x[-5:]*slope + intercept, ls='--', c='g') 

slope, intercept, r_value, p_value, std_error = linregress(x[-5:], y3[-5:],) 
ax.plot(x[-5:], x[-5:]*slope + intercept, ls='--', c='b') 

ax.legend(loc=2) #,fontsize='xx-small') 
#ax.set_xscale('log') 

fig.savefig('figures/modelCurvesFigure.png',format='png') 
plt.close() 

