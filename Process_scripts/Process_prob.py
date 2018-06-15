import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from matplotlib import rc_file

figCols=4
figRows=4

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


if not os.path.isdir('figures') : 
    os.mkdir('figures') 

fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Angle (deg)", ha='center', va='center') 
fig.text(0.03,0.5, r"Probability", ha='center', va='center',rotation='vertical') 

index=0
for mol in molecList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    molec = 'RasRalC18CNC_Q61%s'%mol
    datafile = '%s/Analysis/wham/%s.output.prob'%(molec, molec) 
    print datafile
    try : 
        data = np.genfromtxt(datafile) 
        angles = np.linspace(-180,180,len(data) ) 
        ax.plot(angles,data) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    ax.set_title("Q61%s"%mol) 
    ax.set_xlim([-180,180])
    ax.set_ylim([0,0.05]) 

    index +=1

fig.savefig('figures/wham_probs.pdf',format='pdf') 
plt.close() 



##PMF plots 
fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Angle (deg)", ha='center', va='center') 
fig.text(0.03,0.5, r"PMF (kJ/mol)", ha='center', va='center',rotation='vertical') 

index=0
for mol in molecList : 
    print index, index/figCols, index%figRows
    ax = axarr[index/figCols,index%figCols]

    molec = 'RasRalC18CNC_Q61%s'%mol
    datafile = '%s/Analysis/wham/%s.output.mean'%(molec, molec) 
    print datafile
    try : 
        data = np.genfromtxt(datafile) 
        angles = np.linspace(-180,180,len(data) ) 
        ax.plot(angles,data) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    ax.set_title("Q61%s"%mol) 
    ax.set_xlim([-180,180])
    ax.set_ylim([0,50]) 

    index +=1

fig.savefig('figures/wham_pmf.pdf',format='pdf') 
plt.close() 
