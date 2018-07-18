import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from scipy.stats import linregress
from matplotlib import rc_file
from scipy.optimize import curve_fit

exp_data = 'Exp_data/rates.txt'
exp_data2 = 'Exp_data/transforming_efficiency.txt' 
rcFile = 'rc_files/paper.rc'

nameToColorKeys = {
        }

molecList = [
        "A", 
#        "D",
        "E",
        "F",
        "G", 
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "Q",
        "R",
#        "S",
        "T",
        "V",
        "W",
        "Y"
        ]

rc_file(rcFile)

if not os.path.isdir('figures') : 
    os.mkdir('figures') 

##Load in experimental data into dictionary to match with computed values
nameToExpRate = {} 
with open(exp_data) as f : 
    lines = f.readlines() 
    for line in lines : 
        if line.startswith('#') : 
            continue 
        else : 
            key = line.split()[0]
            value = [float(line.split()[1]), float(line.split()[2])]
            nameToExpRate[key] = value 

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



nameToExpPeak = {} 
with open(exp_data2) as f : 
    lines = f.readlines() 
    for line in lines : 
        if line.startswith('#') : 
            continue 
        else : 
            key = line.split()[0]
            value = float(line.split()[1])
            nameToExpPeak[key] = value 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Intial rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Transorm efficiency", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol]

    peak = nameToExpPeak[name]

    rate, error = nameToExpRate[name] 

    if color == 'b' : 
        avgbAccum = np.append(avgbAccum,peak) 
        ratebAccum = np.append(ratebAccum,rate) 
    elif color == 'r' : 
        avgrAccum = np.append(avgrAccum,peak) 
        raterAccum = np.append(raterAccum,rate) 
    elif color == 'k' : 
        avgkAccum = np.append(avgkAccum,peak) 
        ratekAccum = np.append(ratekAccum,rate) 

    ax.errorbar(rate,peak, xerr=error,marker=marker,color=color, capsize=3) 
    ax.annotate(mol, (rate,peak) ) 

#slope, intercept, r_value, p_value, std_error = linregress(ratebAccum,avgbAccum) 
#x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
#y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

#slope, intercept, r_value, p_value, std_error = linregress(avgrAccum, raterAccum) 
#x = np.linspace(np.min(avgrAccum), np.max(avgrAccum),100) 
#y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='r')
#
#slope, intercept, r_value, p_value, std_error = linregress(avgkAccum, ratekAccum) 
#x = np.linspace(np.min(avgkAccum), np.max(avgkAccum),100) 
#y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=2) 

ax.set_xscale('log')

fig.savefig('figures/rate_v_transform.png',format='png') 
plt.close() 
