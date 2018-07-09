import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from scipy.stats import linregress
from matplotlib import rc_file

exp_data = 'Exp_data/rates.txt'
rcFile = ('rc_files/paper.rc')

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
fig.text(0.03,0.5, r"Num of hbonds to phosphate oxygens from water", ha='center', va='center',rotation='vertical') 

avgAccum, rateAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/hbnum.weighted.out'%molec
    color,marker = plotKeys[mol] 
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

    #avg *= 10 ##nm^2 -> A^2
    #std *= 10 ##nm^2 -> A^2

    rate, error = nameToExpPeak[name] 

    avgAccum = np.append(avgAccum,avg) 
    rateAccum = np.append(rateAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    ax.annotate(mol, (rate, avg) ) 

slope, intercept, r_value, p_value, std_error = linregress(rateAccum, avgAccum) 
x = np.linspace(np.min(rateAccum), np.max(rateAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=4) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_hbnum_gtp.png',format='png') 
plt.close() 



fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Num of hbonds to Q61X s.c. from water", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/hbnum_Q61.weighted.out'%molec
    color,marker = plotKeys[mol] 
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

    #avg *= 10 ##nm^2 -> A^2
    #std *= 10 ##nm^2 -> A^2

    rate, error = nameToExpPeak[name] 

    if color == 'b' : 
        avgbAccum = np.append(avgbAccum,avg) 
        ratebAccum = np.append(ratebAccum,rate) 
    elif color == 'r' : 
        avgrAccum = np.append(avgrAccum,avg) 
        raterAccum = np.append(raterAccum,rate) 
    elif color == 'k' : 
        avgkAccum = np.append(avgkAccum,avg) 
        ratekAccum = np.append(ratekAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    ax.annotate(mol, (rate, avg) ) 

slope, intercept, r_value, p_value, std_error = linregress(ratebAccum, avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(raterAccum, avgrAccum) 
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_hbnum_Q61.png',format='png') 
plt.close() 



fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Num of hbonds to Y32 s.c. from water", ha='center', va='center',rotation='vertical') 

avgAccum, rateAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/hbnum_Y32.weighted.out'%molec
    color,marker = plotKeys[mol] 
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

    #avg *= 10 ##nm^2 -> A^2
    #std *= 10 ##nm^2 -> A^2

    rate, error = nameToExpPeak[name] 

    avgAccum = np.append(avgAccum,avg) 
    rateAccum = np.append(rateAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    ax.annotate(mol, (rate, avg) ) 

slope, intercept, r_value, p_value, std_error = linregress(rateAccum, avgAccum) 
x = np.linspace(np.min(rateAccum), np.max(rateAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=4) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_hbnum_Y32.png',format='png') 
plt.close() 
