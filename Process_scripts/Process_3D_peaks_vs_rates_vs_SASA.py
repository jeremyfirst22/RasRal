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
from mpl_toolkits.mplot3d import Axes3D

exp_data = 'Exp_data/rates.txt'
exp_data2 = 'Exp_data/ftir_peaks.txt'
rcFile = 'rc_files/paper.rc'

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

#rc_file(rcFile)

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
            value = [float(line.split()[1]), float(line.split()[2])]
            nameToExpPeak[key] = value 

fig = plt.figure() 
ax = plt.axes(projection='3d') 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
#fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
#fig.text(0.03,0.5, r"Intial rate (min-1)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = "%s/Analysis/boltzmann/sidechain.weighted.out"%molec
    color, marker = plotKeys[mol]

    peak, error2 = nameToExpPeak[name]

    rate, error = nameToExpRate[name] 

    try : 
        with open(datafile) as f : 
            lines = f.readlines()
            avg = float(lines[2].split()[-2])
            std = float(lines[2].split()[-1])
    except IOError :
        print "No file found for %s"%(datafile)
        continue
    except :
        print "Error importing data from file %s"%(datafile)
        continue

    avg *= 100 ##nm^2 -> A^2
    std *= 100 ##nm^2 -> A^2

    if mol == 'I' or mol == 'N' : 
        marker = '^' 
        color = 'c'
    elif color == 'b' : 
        avgbAccum = np.append(avgbAccum,peak) 
        ratebAccum = np.append(ratebAccum,rate) 
    elif color == 'r' : 
        avgrAccum = np.append(avgrAccum,peak) 
        raterAccum = np.append(raterAccum,rate) 
    elif color == 'k' : 
        avgkAccum = np.append(avgkAccum,peak) 
        ratekAccum = np.append(ratekAccum,rate) 

    #ax.errorbar(peak,rate,xerr=error2, yerr=error,marker=marker,color=color, capsize=3) 
    ax.scatter(peak, avg, np.log(rate), c=color, marker=marker) 
    print molec, peak, avg, rate
    #ax.annotate(mol, (peak+0.03, rate) ) 

slope, intercept, r_value, p_value, std_error = linregress(avgbAccum, ratebAccum) 
x = np.linspace(np.min(avgbAccum), np.max(avgbAccum),100) 
y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(avgrAccum, raterAccum) 
x = np.linspace(np.min(avgrAccum), np.max(avgrAccum),100) 
y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

slope, intercept, r_value, p_value, std_error = linregress(avgkAccum, ratekAccum) 
x = np.linspace(np.min(avgkAccum), np.max(avgkAccum),100) 
y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

#ax.legend(loc=1) 

#ax.set_xlim([2161,2165]) 
#ax.set_ylim([0.0025, 0.12]) 

#ax.set_yscale('log')

fig.savefig('figures/3D_rate_v_peak.png',format='png') 
plt.show() 
plt.close() 
