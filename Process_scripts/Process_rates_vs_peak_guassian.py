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
exp_data2 = 'Exp_data/ftir_peaks.txt'
rcFile = 'rc_files/paper.rc'

nameToColorKeys = {
        }

molecList = [
        "A", 
        "D",
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

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Initial rate ($\mu$M min$^{-1}$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol]

    peak, error2 = nameToExpPeak[name]

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

    ax.errorbar(peak,rate,xerr=error2, yerr=error,marker=marker,color=color, capsize=3) 
    #ax.annotate(mol, (peak+0.03, rate+0.001) ) 

def normDist(x, a, b, c, d) : 
    return a*np.exp(-(x - b)**2/(2 * c**2)) + d

try : 
    popt, pcov = curve_fit(normDist, avgbAccum, np.log(ratebAccum), p0=(np.max(np.log(ratebAccum)), np.average(avgbAccum), np.std(avgbAccum), np.min(np.log(ratebAccum)))) 
    p1, p2, p3, p4, = popt[:]
    print p1, p2, p3, p4
    x = np.linspace(np.min(avgbAccum), np.max(avgbAccum), 100) 
    y = np.exp(normDist(x, p1, p2, p3, p4) ) 
    #y = normDist(x, 0.1, 2163.5, 1)  
    ax.plot(x, y, color='b') 
except RuntimeError : 
    print "No fit found!" 
    x = np.linspace(np.min(avgbAccum), np.max(avgbAccum), 100) 
    y = (normDist(x, np.max(ratebAccum), np.average(avgbAccum), 0.3, np.min(ratebAccum) ))  
    ax.plot(x, y, color='b') 
    

ax.legend(loc=1) 

#ax.set_xlim([2161,2165]) 
#ax.set_ylim([0.0025, 0.12]) 

ax.set_yscale('log')

fig.savefig('figures/rate_v_peak_guassian.png',format='png') 
plt.close() 

