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

rcFile = 'rc_files/paper.rc'
exp_data = 'Exp_data/ftir_peaks.txt'

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
#        "X",
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

def logFit(x, a, b) :
    return a*np.log(x) + b



fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(left=0.15,right=0.95, bottom=0.12,top=0.95) 
fig.text(0.55,0.04, r"Initial rate ($\muup$M min$^{-1}$)", ha='center', va='center') 
#fig.text(0.03,0.5, r"Num waters within %s A of Q61 s.c. and %s A of O1G"%(dist, dist) , ha='center', va='center',rotation='vertical') 
fig.text(0.03,0.535, r"Average number of waters", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/size_%s.weighted.out'%(molec,dist) 
    #print datafile

    rate, error = nameToExpRate[name] 
    
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

    avg /= 3 ##atoms -> water molecules
    std /= 3 ##atoms -> water molecules 

    if color == 'b' : 
        avgbAccum = np.append(avgbAccum,avg) 
        ratebAccum = np.append(ratebAccum,rate) 
    elif color == 'r' : 
        marker = '^' 
        avgrAccum = np.append(avgrAccum,avg) 
        raterAccum = np.append(raterAccum,rate) 
    elif color == 'k' : 
        marker = '^' 
        avgkAccum = np.append(avgkAccum,avg) 
        ratekAccum = np.append(ratekAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
    #ax.annotate(mol, (rate, avg) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

try : 
    slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
    x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
    y = slope * np.log(x) + intercept 
    ax.plot(x, y, label = "$r = %.2f$"%r_value,color='b')
except ValueError : 
    print "No data found for blue" 

 try : 
     slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
     x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
     y = slope * np.log(x) + intercept 
     ax.plot(x, y, label = "r = %.2f"%r_value,color='r')
 except ValueError : 
     print "No data found for red" 
 
 try : 
     slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
     x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
     y = slope * np.log(x) + intercept 
     ax.plot(x, y, label = "r = %.2f"%r_value,color='k')
 except ValueError : 
     print "No data found for black" 

ax.legend(loc=4) 
ax.set_xscale('log') 

fig.savefig('figures/peak_vs_cnc_scount.png'%dist,format='png') 

