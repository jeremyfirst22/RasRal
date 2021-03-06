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
exp_data = 'Exp_data/dockedFigure.txt'

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
y2    = data[:,2] 
yerr1 = data[:,3] 
yerr2 = data[:,4] 

fig, ax = plt.subplots(1,1)
fig.subplots_adjust(left=0.15,right=0.95, bottom=0.12,top=0.95)
fig.text(0.55,0.04, r"Time (min)", ha='center', va='center')
fig.text(0.03,0.535, r"[ $\rm{P_i}$ ] ($\rm{\muup M})$", ha='center', va='center',rotation='vertical') 

ax.scatter(x, y1, marker='o', color='b',label=r"WT Ras") 
ax.scatter(x, y2, marker='o', color='g',label=r"WT Ras+Ral$\betaup$I18C$_{\rm{SCN}}$") 

ax.errorbar(x, y1, yerr=yerr1 , marker='o', color='b' ,capsize=3, ls = 'None') 
ax.errorbar(x, y2, yerr=yerr2 , marker='o', color='g' ,capsize=3, ls = 'None') 

#ax1 = fig.add_axes([0.35, 0.20, 0.25, 0.25], frameon='False') 
#ax1.spines['right'].set_visible(False)
#ax1.spines['top'].set_visible(False)
#ax1.errorbar(x, y1, yerr=yerr1 , marker='o', markersize=1.5, color='b' ,capsize=1.5, ls = 'None', elinewidth=0.5, capthick = 0.5 ) 
slope, intercept, r_value, p_value, std_error = linregress(x[-5:], y1[-5:]) 
#ax.plot(x[-5:], x[-5:]*slope+intercept, ls='--', c='b', lw = 0.5, label=r"%.1f $\pm$ %.1f"%(slope*100, std_error * 100) ) 
ax.plot(x[-5:], x[-5:]*slope+intercept, ls='--', c='b') 
##ax1.legend(loc='upper center',fontsize=5,frameon='None',fancybox='None', numpoints='None') 
ax.text(0.95,0.30,r"%.1f $\pm$ %.1f x 10$^{-2} \muup$M min$^{-1}$"%(slope*100, std_error * 100),transform=ax.transAxes,ha='right', va='bottom',color='b') 
#ax1.set_xticks([])
#ax1.set_yticks([])

#ax2 = fig.add_axes([0.65, 0.20, 0.25, 0.25], frame_on='False') 
#ax2.spines['right'].set_visible(False)
#ax2.spines['top'].set_visible(False)
#ax2.errorbar(x, y2, yerr=yerr2 , marker='o', markersize=1.5, color='g' ,capsize=1.5, ls = 'None', elinewidth=0.5, capthick = 0.5 ) 
slope, intercept, r_value, p_value, std_error = linregress(x[-5:], y2[-5:]) 
#ax.plot(x[-5:], x[-5:]*slope+intercept, ls='--', c='g', lw = 0.5, label=r"%.1f $\pm$ %.1f"%(slope*100, std_error * 100) ) 
ax.plot(x[-5:], x[-5:]*slope+intercept, ls='--', c='g') 
##ax1.legend(loc='upper center',fontsize=5,frameon='None',fancybox='None', numpoints='None') 
ax.text(0.95,0.30,r"%.1f $\pm$ %.1f x 10$^{-2}\muup$M min$^{-1}$"%(slope*100, std_error * 100),transform=ax.transAxes,ha='right', va='top',color='g') 
#ax2.set_xticks([])
#ax2.set_yticks([])
#
ax.legend(loc=4) #,fontsize='x-small') 
##ax.set_xscale('log') 

fig.savefig('figures/dockedFigure.png',format='png') 
plt.close() 

