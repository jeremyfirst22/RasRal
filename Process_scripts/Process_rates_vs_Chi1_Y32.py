import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from scipy.stats import linregress
from matplotlib import rc_file

exp_data = 'Exp_data/rates2.txt'
rcFile = 'rc_files/paper.rc'

molecList = [
        "D",
        "E",
        "F",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y"
        ]

rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

##Load in experimental data into dictionary to match with computed values
nameToExpPeak = {} 
with open(exp_data) as f : 
    lines = f.readlines() 
    for line in lines : 
        if line.startswith('#') : 
            continue 
        else : 
            key = line.split()[0]
            value = [float(line.split()[1]), float(line.split()[2])]
            nameToExpPeak[key] = value 

##Load in color and marker keys 
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


fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Avg. Chi1 angle of Y32 (deg)",ha='center', va='center',rotation='vertical') 

avgAccum, rateAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color = plotKeys[mol][0]
    marker = plotKeys[mol][1]

    datafile = '%s/Analysis/boltzmann/angaver_Y32.weighted.out'%molec
    #print datafile
    try : 
        with open(datafile) as f : 
            lines = f.readlines() 
            avg = float(lines[1].split()[-2]) 
            std = float(lines[1].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile)  
        continue 
    except : 
        print "Error importing data from file %s"%(datafile)
        continue 

    #avg *= 1000 ##nm^3 -> A^3
    #std *= 1000 ##nm^3 -> A^3

    rate, error = nameToExpPeak[name] 

    #if mol == "S" or mol == "T" or mol == "N" or mol == "Q" or mol == "Y" : 
    if True : 
        avgAccum = np.append(avgAccum,avg) 
        rateAccum = np.append(rateAccum,rate) 
        marker = '^'
    else : 
        #avgAccum = np.append(avgAccum,avg) 
        #rateAccum = np.append(rateAccum,rate) 
        marker = 'o'

    #ax.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    ax.errorbar(rate, avg,xerr=error, marker=marker,color=color,capsize=3) 
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, avg) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(rateAccum, avgAccum) 
x = np.linspace(np.min(rateAccum), np.max(rateAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_Chi1_Y32.png',format='png') 
plt.close() 



