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
exp_data = 'Exp_data/optical_response.txt' 

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

fig, ax = plt.subplots(1,1,figsize=[3,1.5]) 
fig.subplots_adjust(left=0.15,right=0.95, bottom=0.22,top=0.95) 
fig.text(0.55,0.04, r"[P$_{\rm{i}}$] ($\muup$M)",ha='center', va='center') 
fig.text(0.03,0.585, r"Absorbance (a.u.)", ha='center', va='center',rotation='vertical') 

ax.scatter(x, y1, marker='o', color='k') 
ax.errorbar(x, y1, yerr=yerr1 , marker='o', color='k', ls = 'None') 

slope, intercept, r_value, p_value, std_error = linregress(x, y1) 
ax.plot(x, x*slope + intercept, ls='-', c='k',label= r"r = %0.3f"%r_value) 

ax.legend(loc=2)#,fontsize=9) #,fontsize='xx-small') 
#ax.set_xscale('log') 

#ax.set_ylim([-1.5,19.5]) 

fig.savefig('figures/optical_response.png',format='png') 
plt.close() 

