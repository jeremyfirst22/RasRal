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
exp_data = 'Exp_data/ftir_peaks.txt'
exp_data2= 'Exp_data/rates.txt'

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
#        "X", 
        "Y"
        ]

rc_file(rcFile) 

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

nameToExpRate = {}
with open(exp_data2) as f :
    lines = f.readlines()
    for line in lines :
        if line.startswith('#') :
            continue
        else :
            key = line.split()[0]
            value = [float(line.split()[1]), float(line.split()[2])]
            nameToExpRate[key] = value

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


fig, axarr = plt.subplots(1,2,figsize=[6,2.5]) 
fig.subplots_adjust(wspace=0.20,hspace=0.00,top=0.95,bottom=0.15,left=0.08,right=0.98) 
fig.text(0.280,0.050, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig.text(0.780,0.050, r"Initial rate ($\muup$M min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.550, r"Polar atom SASA ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 
fig.text(0.515,0.550, r"Side chain SASA ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 
#fig.text(0.03,0.27, r"Sidechain", ha='center', va='center',rotation='vertical') 
#fig.text(0.03,0.75, r"Polar atom", ha='center', va='center',rotation='vertical') 

fig2, axarr2 = plt.subplots(1,2,figsize=[6,2.5]) 
fig2.subplots_adjust(wspace=0.20,hspace=0.00,top=0.95,bottom=0.15,left=0.08,right=0.98) 
fig2.text(0.780,0.050, r"Absorption frequency (cm$^{-1}$)", ha='center', va='center') 
fig2.text(0.280,0.050, r"Initial rate ($\muup$M min$^{-1}$)", ha='center', va='center') 
fig2.text(0.030,0.550, r"Polar atom SASA ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 
fig2.text(0.515,0.550, r"Side chain SASA ($\rm{\AA}^2$)", ha='center', va='center',rotation='vertical') 
#fig.text(0.03,0.27, r"Sidechain", ha='center', va='center',rotation='vertical') 
#fig.text(0.03,0.75, r"Polar atom", ha='center', va='center',rotation='vertical') 

ax1 = axarr[0] 
avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color,marker = plotKeys[mol]
    datafile = '%s/Analysis/boltzmann/sc_polar.weighted.out'%molec
    #print datafile
    try : 
        with open(datafile) as f : 
            lines = f.readlines() 
            avg = float(lines[2].split()[-2]) 
            std = float(lines[2].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile)  
        avg, std = 0,0
    except : 
        print "Error importing data from file %s"%(datafile)

    avg *= 100 ##nm^2 -> A^2
    std *= 100 ##nm^2 -> A^2

    peak, error = nameToExpPeak[name] 

    if mol == "M" or mol == "F" or mol == "L" or mol == "V" or mol == "I" : 
        #avg = 0 
        #std = 0 
        marker = '^'
    else : 
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        marker = 'o'

    ax1.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color=color) 
    #ax.scatter(peak, avg) 
    #ax1.annotate(mol, (peak+0.02, avg+3) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 
ax1.plot(x, y, label = "$r = %.2f$"%r_value,color='k')
ax1.legend(loc=1) 
ax1.set_xlim([2161.8,2165.2]) 
#ax1.text(0.9,0.9,'A') 
#ax1.text(0.05,0.95,"A",transform=ax1.transAxes,ha='left',va='top') 

ax2 = axarr2[1]
avgAccum, peakAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color,marker = plotKeys[mol]
    datafile = '%s/Analysis/boltzmann/sidechain.weighted.out'%molec
    #print datafile
    try : 
        with open(datafile) as f : 
            lines = f.readlines() 
            avg = float(lines[2].split()[-2]) 
            std = float(lines[2].split()[-1])
    except IOError : 
        print "No file found for %s"%(datafile)  
        avg, std = 0,0
    except : 
        print "Error importing data from file %s"%(datafile)

    avg *= 100 ##nm^2 -> A^2
    std *= 100 ##nm^2 -> A^2

    peak, error = nameToExpPeak[name] 

    if False : #mol == "M" or mol == "F" or mol == "L" or mol == "V" or mol == "I" : 
        #avg = 0 
        #std = 0 
        marker = '^'
    else : 
        avgAccum = np.append(avgAccum,avg) 
        peakAccum = np.append(peakAccum,peak) 
        marker = 'o'

    ax2.errorbar(peak, avg,xerr=error, yerr=std,marker=marker,color=color) 
    #ax.scatter(peak, avg) 
    #ax2.annotate(mol, (peak+0.02, avg+3) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(peakAccum, avgAccum) 
x = np.linspace(np.min(peakAccum), np.max(peakAccum),100) 
y = slope * x + intercept 
ax2.plot(x, y, label = "$r = %.2f$"%r_value,color='k')
ax2.legend(loc=3) 
#ax2.text(0.05,0.95,"B",transform=ax2.transAxes,ha='left',va='top') 
ax2.set_xlim([2161.8,2165.2]) 

ax3 = axarr2[0]
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

    try :
        with open(datafile) as f :
            lines = f.readlines()
            avg = float(lines[2].split()[-2])
            std = float(lines[2].split()[-1])
    except IOError :
        print "No file found for %s"%(datafile)
        avg, std = 0, 0 
    except :
        print "Error importing data from file %s"%(datafile)
        continue

    avg *= 100 ##nm^2 -> A^2
    std *= 100 ##nm^2 -> A^2

    #if mol == 'W' : #or mol == 'D' :
    #    marker = '^'
    #    color = 'c'
    if color == 'b' :
        avgbAccum = np.append(avgbAccum,avg)
        ratebAccum = np.append(ratebAccum,rate)
    elif color == 'r' :
        marker = '^' 
    #    avgrAccum = np.append(avgrAccum,avg)
    #    raterAccum = np.append(raterAccum,rate)
    elif color == 'k' :
        marker = '^'
        #avgkAccum = np.append(avgkAccum,avg)
        #ratekAccum = np.append(ratekAccum,rate)

    ax3.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color) 
    #ax.scatter(rate, avg)
    #ax3.annotate(mol, (rate, avg+1) )

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std)

    #index +=1

try :
    slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum)
    x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100)
    y = slope * np.log(x) + intercept
    ax3.plot(x, y, label = "$r = %.2f$"%r_value,color='b')
except ValueError :
    print "No blue"

try :
    slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum)
    x = np.linspace(np.min(raterAccum), np.max(raterAccum),100)
    y = slope * np.log(x) + intercept
    ax3.plot(x, y, label = "$r = %.2f$"%r_value,color='r')
except ValueError :
    print "No red"

try :
    slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum)
    x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100)
    y = slope * np.log(x) + intercept
    ax3.plot(x, y, label = "$r = %.2f$"%r_value,color='k')
except ValueError :
    print "No black"

ax3.legend(loc=1)
#ax3.axes.get_yaxis().set_ticklabels([]) 
#ax3.text(0.05,0.95,"C",transform=ax3.transAxes,ha='left',va='top') 
ax3.set_xscale('log')

ax4 = axarr[1] 
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

    if mol == 'W' : #or mol == 'D' :
        marker = '^'
        color = 'c'
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

    ax4.errorbar(rate, avg,xerr=error, yerr=std,marker=marker,color=color) 
    #ax.scatter(rate, avg)
    #ax4.annotate(mol, (rate, avg+1) )

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std)

    #index +=1

try :
    slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum)
    x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100)
    y = slope * np.log(x) + intercept
    ax4.plot(x, y, label = "$r = %.2f$"%r_value,color='b')
except ValueError :
    print "No blue"

#try :
#    slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum)
#    x = np.linspace(np.min(raterAccum), np.max(raterAccum),100)
#    y = slope * np.log(x) + intercept
#    ax4.plot(x, y, label = "r = %.3f"%r_value,color='r')
#except ValueError :
#    print "No red"

#try :
#    slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum)
#    x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100)
#    y = slope * np.log(x) + intercept
#    ax4.plot(x, y, label = "r = %.3f"%r_value,color='k')
#except ValueError :
#    print "No black"

ax4.legend(loc=4)
#ax4.axes.get_yaxis().set_ticklabels([]) 
#ax4.set_xlim([0.9*10**-3,1.3*10**-1]) 
#ax4.text(0.05,0.95,"D",transform=ax4.transAxes,ha='left',va='top') 
ax4.set_xscale('log')


fig.savefig('figures/paper_sasa_figure.png',format='png') 
fig2.savefig('figures/paper_sasa_si_figure.png',format='png') 
plt.close() 




