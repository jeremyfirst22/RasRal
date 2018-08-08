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
fig.text(0.5,0.04, r"SASAs (A22)", ha='center', va='center') 
fig.text(0.03,0.5, r"Volume of Q61X s.c. ($\rm{\AA}^3$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/isolated_sc.weighted.out'%molec
    datafile2 = '%s/Analysis/boltzmann/sidechain.weighted.out'%molec
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
    
    try : 
        with open(datafile2) as f : 
            lines = f.readlines() 
            avg2= float(lines[2].split()[-2]) 
            std2= float(lines[2].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile)  
        continue 
    except : 
        print "Error importing data from file %s"%(datafile)
        continue 

    avg *= 10**3 ##nm^3 -> A^3
    std *= 10**3 ##nm^3 -> A^3
    avg2*= 10**2 ##nm^2 -> A^2
    std2*= 10**2 ##nm^2 -> A^2

    if True :# color == 'b' : 
        avgbAccum = np.append(avgbAccum,avg) 
        ratebAccum = np.append(ratebAccum,avg2) 

    ax.errorbar(avg2, avg,xerr=std2, yerr=std,marker=marker,color=color,capsize=3) 
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (avg2, avg+1) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(ratebAccum, avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

ax.legend(loc=1) 
#ax.set_xscale('log') 

fig.savefig('figures/volume_vs_sasa.png',format='png') 
plt.close() 
