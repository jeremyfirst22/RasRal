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
exp_data = 'Exp_data/ftir_peaks.txt'
#exp_data = 'Exp_data/rates2.txt'

nameToColorKeys = {
        }

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
#fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.05,0.5, r"Ext. field at nitrile ($\frac{k_B T}{e^- \AA}$)", ha='center', va='center',rotation='vertical') 

avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/external_field_nitrile.weighted.out'%molec
    color,marker = plotKeys[mol] 
    #print datafile
    try : 
        with open(datafile) as f : 
            lines = f.readlines() 
            avg = float(lines[0].split()[-2]) 
            std = float(lines[0].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile)  
        continue  
    except : 
        print "Error importing data from file %s"%(datafile)
        continue 

#    avg *= 100 ##nm^2 -> A^2
#    std *= 100 ##nm^2 -> A^2

    peak, error = nameToExpPeak[name] 

    if mol == "M" or mol == "F" or mol == "L" or mol == "V" or mol == "I" : 
        #avg = 0 
        #std = 0 
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        #marker = '^'
    else : 
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        #marker = 'o'

    ax.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    print peak, avg
    #ax.scatter(peak, avg) 
    ax.annotate(mol, (peak, avg) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
#ax.set_xscale('log') 

fig.savefig('figures/peak_v_external_field_nitrile.png',format='png') 
plt.close() 

