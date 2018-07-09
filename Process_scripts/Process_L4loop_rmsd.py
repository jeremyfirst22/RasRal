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
fig.text(0.5,0.04, r"Residue", ha='center', va='center') 
fig.text(0.03,0.5, r"RMSD of l4loop ($\rm{\AA}$)", ha='center', va='center',rotation='vertical') 

for index,mol in enumerate(molecList) : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/l4loop.weighted.out'%molec
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

    avg *= 10 ##nm -> A
    std *= 10 ##nm -> A

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
#    ax.axhline(avg, color=color,linestyle='--') 
    ax.annotate(mol, (index, 0.8) )

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1
ax.set_xticks(np.arange(len(molecList)), (molecList) ) 
fig.savefig('figures/bar_rmsd_l4loop.png',format='png') 
plt.close() 
