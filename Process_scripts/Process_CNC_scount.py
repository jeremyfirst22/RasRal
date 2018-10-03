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
exp_data2 = 'Exp_data/rates.txt'

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

nameToExpRate = {}
with open(exp_data2) as f :
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

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Residue", ha='center', va='center') 
fig.text(0.03,0.5, r"Avg Num waters near CNC", ha='center', va='center',rotation='vertical') 

for index,mol in enumerate(molecList) : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/cnc_scount.weighted.out'%molec
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

    avg /= 3 ##atoms->molecules
    std /= 3 ##atoms->molecules

    #if mol == 'W' or mol == 'D' : 
    #    marker = '^'
    #    color = 'c' 
    #if color == 'b' : 
    #    avgbAccum = np.append(avgbAccum,avg) 
    #    ratebAccum = np.append(ratebAccum,rate) 
    #elif color == 'r' : 
    #    avgrAccum = np.append(avgrAccum,avg) 
    #    raterAccum = np.append(raterAccum,rate) 
    #elif color == 'k' : 
    #    avgkAccum = np.append(avgkAccum,avg) 
    #    ratekAccum = np.append(ratekAccum,rate) 

    ax.bar(index, avg, yerr=std,color=color,capsize=3) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1
fig.savefig('figures/bar_scount_cnc.png',format='png') 
plt.close() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Avg Num waters near CNC", ha='center', va='center',rotation='vertical') 

avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color,marker = plotKeys[mol]
    datafile = '%s/Analysis/boltzmann/cnc_scount.weighted.out'%molec
    print datafile
    try : 
        with open(datafile) as f : 
            lines = f.readlines() 
            avg = float(lines[1].split()[-2]) 
            std = float(lines[1].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    avg /= 3 ##atoms->molecules
    std /= 3 ##atoms->molecules

    peak, error = nameToExpPeak[name] 

    if mol == "M" or mol == "F" or mol == "L" or mol == "V" or mol == "I" : 
        #avg = 0 
        #std = 0 
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        marker = '^'
    else : 
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        marker = 'o'

    ax.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    #ax.scatter(peak, avg) 
    ax.annotate(mol, (peak+0.02, avg+3) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 

fig.savefig('figures/peak_v_cnc_scount.png',format='png') 


