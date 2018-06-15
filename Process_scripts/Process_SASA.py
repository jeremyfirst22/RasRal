import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from scipy.stats import linregress
from matplotlib import rc_file

exp_data = 'Exp_data/ftir_peaks.txt'

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


fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (whole residue) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/area.weighted.out'%molec
    #print datafile
    try : 
        with open(datafile) as f : 
            lines = f.readlines() 
            avg = float(lines[2].split()[-2]) 
            std = float(lines[2].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    avg *= 100 ##nm^2 -> A^2
    std *= 100 ##nm^2 -> A^2

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

    ax.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color='k',capsize=3) 
    #ax.scatter(peak, avg) 
    ax.annotate(mol, (peak+0.02, avg+3) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 

fig.savefig('figures/peak_v_sasa_whole_residue.pdf',format='pdf') 
plt.close() 



fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (side chain) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/sidechain.weighted.out'%molec
    print datafile
    try : 
        with open(datafile) as f : 
            lines = f.readlines() 
            avg = float(lines[2].split()[-2]) 
            std = float(lines[2].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    avg *= 100 ##nm^2 -> A^2
    std *= 100 ##nm^2 -> A^2

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

    ax.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color='k',capsize=3) 
    #ax.scatter(peak, avg) 
    ax.annotate(mol, (peak+0.02, avg+3) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 

fig.savefig('figures/peak_v_sasa_sidechain.pdf',format='pdf') 



fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (polar atoms only) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/polar.weighted.out'%molec
    print datafile

    peak, error = nameToExpPeak[name] 
    
    if False : #mol == "M" or mol == "F" or mol == "L" or mol == "V" or mol == "I" : 
        avg, std = 0, 0 
        marker = '^'
    else : 
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
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        marker = 'o'

    ax.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color='k',capsize=3) 
    #ax.scatter(peak, avg) 
    ax.annotate(mol, (peak+0.02, avg+3) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 

fig.savefig('figures/peak_v_sasa_polar.pdf',format='pdf') 
plt.close() 
plt.close() 



fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (side chain polar atoms only) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/sc_polar.weighted.out'%molec
    #print datafile

    peak, error = nameToExpPeak[name] 
    
    if mol == "M" or mol == "F" or mol == "L" or mol == "V" or mol == "I" : 
        avg, std = 0, 0 
        marker = '^'
    else : 
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
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        marker = 'o'

    ax.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color='k',capsize=3) 
    #ax.scatter(peak, avg) 
    ax.annotate(mol, (peak+0.02, avg+3) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 

fig.savefig('figures/peak_v_sasa_sc_polar.pdf',format='pdf') 
plt.close() 
plt.close() 



fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (David's sidechain polar atoms only) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/boltzmann/davids.weighted.out'%molec
    #print datafile

    peak, error = nameToExpPeak[name] 
    
    if mol == "M" or mol == "F" or mol == "L" or mol == "V" or mol == "I" : 
        avg, std = 0, 0 
        marker = '^'
    else : 
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

        avg *= 100 ##nm^2 -> A^2
        std *= 100 ##nm^2 -> A^2
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        marker = 'o'

    ax.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color='k',capsize=3) 
    #ax.scatter(peak, avg) 
    ax.annotate(mol, (peak+0.02, avg+3) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 

fig.savefig('figures/peak_v_sasa_davids.pdf',format='pdf') 
plt.close() 
plt.close() 
