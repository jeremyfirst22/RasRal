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
exp_data = 'Exp_data/rates.txt'

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

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (whole residue) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/area.weighted.out'%molec
    #print datafile

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

    if mol == 'W' or mol == 'D' : 
        marker = '^'
        color = 'c' 
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
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, avg+1) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(ratebAccum, avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(raterAccum, avgrAccum) 
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

slope, intercept, r_value, p_value, std_error = linregress(ratekAccum, avgkAccum) 
x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_sasa_whole_residue.png',format='png') 
plt.close() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (side chain) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/sidechain.weighted.out'%molec
    #print datafile

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

    if mol == 'W' or mol == 'D' : 
        marker = '^'
        color = 'c' 
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
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, avg+1) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_sasa_sidechain.png',format='png') 
plt.close() 


fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (side chain) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/isolated_sc.weighted.out'%molec
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

    avg *= 100 ##nm^2 -> A^2
    std *= 100 ##nm^2 -> A^2

    if mol == 'W' or mol == 'D' : 
        marker = '^'
        color = 'c' 
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
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, avg+1) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_sasa_isolated.png',format='png') 
plt.close() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (side chain) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/sidechain.weighted.out'%molec
    datafile2 = '%s/Analysis/boltzmann/isolated_sc.weighted.out'%molec
    #print datafile

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

    try : 
        with open(datafile2) as f : 
            lines = f.readlines() 
            avgIso = float(lines[1].split()[-2]) 
            stdIso = float(lines[1].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile2)  
        continue 
    except : 
        print "Error importing data from file %s"%(datafile2)
        continue 

    avg /= avgIso ##nm^2 -> A^2
    #std /= stdIso ##nm^2 -> A^2

    if mol == 'W' or mol == 'D' : 
        marker = '^'
        color = 'c' 
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
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, avg+1) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
y = slope * np.log(x) + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_sasa_relative_sidechain.png',format='png') 
plt.close() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (polar atoms only) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/polar.weighted.out'%molec
    #print datafile

    rate, error = nameToExpRate[name] 
    
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

        #if mol == 'W' : 
        #    marker = '^'
        #    color = 'c' 
        if color == 'b' : 
            avgbAccum = np.append(avgbAccum,avg) 
            ratebAccum = np.append(ratebAccum,rate) 
        elif color == 'r' : 
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        #elif color == 'k' : 
        #    avgkAccum = np.append(avgkAccum,avg) 
        #    ratekAccum = np.append(ratekAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, avg+1) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(ratebAccum, avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(raterAccum, avgrAccum) 
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

#slope, intercept, r_value, p_value, std_error = linregress(ratekAccum, avgkAccum) 
#x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
#y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_sasa_polar.png',format='png') 
plt.close() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (side chain polar atoms only) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/sc_polar.weighted.out'%molec
    #print datafile

    rate, error = nameToExpRate[name] 
    
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

        #if mol == 'W' : 
        #    marker = '^'
        #    color = 'c' 
        if color == 'b' : 
            avgbAccum = np.append(avgbAccum,avg) 
            ratebAccum = np.append(ratebAccum,rate) 
        elif color == 'r' : 
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        #elif color == 'k' : 
        #    rate, std = 0,0
        #    avgkAccum = np.append(avgkAccum,avg) 
        #    ratekAccum = np.append(ratekAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, avg+1) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(ratebAccum, avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(raterAccum, avgrAccum) 
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

#slope, intercept, r_value, p_value, std_error = linregress(ratekAccum, avgkAccum) 
#x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
#y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_sasa_sc_polar.png',format='png') 
plt.close() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"SASA (David's side chain polar atoms only) ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/davids.weighted.out'%molec
    #print datafile

    rate, error = nameToExpRate[name] 
    
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

        #if mol == 'W' : 
        #    marker = '^'
        #    color = 'c' 
        if color == 'b' : 
            avgbAccum = np.append(avgbAccum,avg) 
            ratebAccum = np.append(ratebAccum,rate) 
        elif color == 'r' : 
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        #elif color == 'k' : 
        #    rate, std = 0,0
        #    avgkAccum = np.append(avgkAccum,avg) 
        #    ratekAccum = np.append(ratekAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color,capsize=3) 
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, avg+1) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(ratebAccum, avgbAccum) 
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='b')

slope, intercept, r_value, p_value, std_error = linregress(raterAccum, avgrAccum) 
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
y = slope * x + intercept 
ax.plot(x, y, label = "r = %.3f"%r_value,color='r')

#slope, intercept, r_value, p_value, std_error = linregress(ratekAccum, avgkAccum) 
#x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
#y = slope * x + intercept 
#ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_sasa_davids.png',format='png') 
plt.close() 


