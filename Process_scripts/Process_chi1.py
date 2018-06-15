import matplotlib.pyplot as plt 
import numpy as np 
import os as os
import glob as glob 
from sys import exit 


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
        "Y"]

saveDir = 'figures'

if not os.path.isdir(saveDir) : 
    os.mkdir(saveDir) 

for molec in molecList : 
    fig, ax = plt.subplots(1,1) 
    ax.set_xlabel('Angle (deg)') 
    ax.set_ylabel('Occurances') 

    for angle in np.arange(0,359,30) : 
        fileName = 'RasRalC18CNC_Q61%s/Analysis/chi1/angdist.%i.xvg'%(molec, angle) 
        if os.path.isfile(fileName) : 
            headlines = 0 
            with open(fileName) as f : 
                lines = f.readlines() 
            for line in lines : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        else : 
            print "Warning: No file %s found" %fileName 
            continue 

        try : 
            data = np.genfromtxt(fileName, skip_header=headlines) 
        except : 
            print "ERROR: Data import failed for %s" %fileName 
            sys.exit() 

        ax.plot(data[:,0],data[:,1],label=angle) 

    fig.savefig("%s/chi1_hist_%s.pdf"%(saveDir, molec),format='pdf' ) 
    plt.close() 

    fig, ax = plt.subplots(1,1) 
    ax.set_xlabel('Time (ns)') 
    ax.set_ylabel('Angle (deg)') 

    angAccum = [] 
    angAccum = np.array(angAccum) 

    for angle in np.arange(0,359,30) : 
        fileName = 'RasRalC18CNC_Q61%s/Analysis/chi1/angaver.%i.xvg'%(molec, angle) 
        if os.path.isfile(fileName) : 
            headlines = 0 
            with open(fileName) as f : 
                lines = f.readlines() 
            for line in lines : 
                if line.startswith('#') or line.startswith('@') : 
                    headlines += 1 
                else : 
                    break 
        else : 
            print "Warning: No file %s found" %fileName 
            continue 

        try : 
            data = np.genfromtxt(fileName, skip_header=headlines) 
        except : 
            print "ERROR: Data import failed for %s" %fileName 
            sys.exit() 

        ax.scatter(data[:,0],data[:,1],label=angle, s = 0.5) 
        angAccum = np.append(angAccum, data[:,1]) 

    fig.savefig("%s/chi1_time_%s.png"%(saveDir, molec) , format='png') 
    plt.close() 

    fig, ax = plt.subplots(1,1) 
    ax.set_xlabel('Angle (deg)') 
    ax.set_ylabel('Total occurances') 

    hist, bins = np.histogram(angAccum, bins = np.linspace(-180,180,360) ) 
    ax.plot(hist) 
    fig.savefig("%s/chi1_total_%s.pdf"%(saveDir,molec),format='pdf') 




