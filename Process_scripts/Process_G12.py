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
#fig.text(0.5,0.04, r"Residue", ha='center', va='center') 
fig.text(0.03,0.5, r"G12 - OS3 Distance ($\rm{\AA}$)", ha='center', va='center',rotation='vertical') 

for index,mol in enumerate(molecList) : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/g12_distave.weighted.out'%molec
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

    avg *= 10 #nm -> A
    std *= 10 #nm -> A

    ax.bar(index, avg, yerr=std,color=color,capsize=3) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

ax.set_xticks(np.arange(len(molecList))) 
ax.set_xticklabels(molecList) 

fig.savefig('figures/bar_g12_dist.png') 
plt.close() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
#fig.text(0.5,0.04, r"Residue", ha='center', va='center') 
fig.text(0.03,0.5, r"Avg. Num. of Hbond per frame", ha='center', va='center',rotation='vertical') 

for index,mol in enumerate(molecList) : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/g12_hbnum.weighted.out'%molec
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

    ax.bar(index, avg, yerr=std,color=color,capsize=3) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

ax.set_xticks(np.arange(len(molecList))) 
ax.set_xticklabels(molecList) 

fig.savefig('figures/bar_g12_hbnum.png') 
plt.close() 
