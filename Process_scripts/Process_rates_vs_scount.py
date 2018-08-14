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


for dist in ['245', '3', '4', '5', '6'] : 
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
        ax.plot(x, y, label = "r = %.3f"%r_value,color='b')
    except ValueError : 
        print "No data found for blue" 
    
#    try : 
#        slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
#        x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
#        y = slope * np.log(x) + intercept 
#        ax.plot(x, y, label = "r = %.3f"%r_value,color='r')
#    except ValueError : 
#        print "No data found for red" 
#    
#    try : 
#        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
#        x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
#        y = slope * np.log(x) + intercept 
#        ax.plot(x, y, label = "r = %.3f"%r_value,color='k')
#    except ValueError : 
#        print "No data found for black" 
    
    ax.legend(loc=4) 
    ax.set_xscale('log') 
    
    fig.savefig('figures/rate_v_scount_%s.png'%dist,format='png') 
sys.exit() 

for dist in ['245', '3', '4', '5', '6'] : 
    fig, ax = plt.subplots(1,1) 
    fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
    fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
    fig.text(0.03,0.5, r"Num waters within %s A of Q61 s.c. and %s A of O2G"%(dist, dist) , ha='center', va='center',rotation='vertical') 
    
    avgbAccum, ratebAccum = [], [] 
    avgrAccum, raterAccum = [], [] 
    avgkAccum, ratekAccum = [], [] 
    for mol in molecList : 
        molec = 'RasRalC18CNC_Q61%s'%mol
        name = "Q61%s"%mol 
        color, marker = plotKeys[mol] 
        datafile = '%s/Analysis/boltzmann/O2G_size_%s.weighted.out'%(molec,dist) 
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
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        elif color == 'k' : 
            avgkAccum = np.append(avgkAccum,avg) 
            ratekAccum = np.append(ratekAccum,rate) 
    
        ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
        ax.annotate(mol, (rate, avg) ) 
    
        print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
        x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='b')
    except ValueError : 
        print "No data found for blue" 
    
#    try : 
#        slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
#        x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
#        y = slope * np.log(x) + intercept 
#        ax.plot(x, y, label = "r = %.3f"%r_value,color='r')
#    except ValueError : 
#        print "No data found for red" 
#    
#    try : 
#        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
#        x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
#        y = slope * np.log(x) + intercept 
#        ax.plot(x, y, label = "r = %.3f"%r_value,color='k')
#    except ValueError : 
#        print "No data found for black" 
    
    ax.legend(loc=4) 
    ax.set_xscale('log') 
    
    fig.savefig('figures/rate_v_O2G_scount_%s.png'%dist,format='png') 
    plt.close() 
sys.exit() 

for dist in ['245', '3', '4', '5', '6'] : 
    fig, ax = plt.subplots(1,1) 
    fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
    fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
    fig.text(0.03,0.5, r"Num waters within 5 A of Q61 s.c. and %s A of O1G"%(dist) , ha='center', va='center',rotation='vertical') 
    
    avgbAccum, ratebAccum = [], [] 
    avgrAccum, raterAccum = [], [] 
    avgkAccum, ratekAccum = [], [] 
    for mol in molecList : 
        molec = 'RasRalC18CNC_Q61%s'%mol
        name = "Q61%s"%mol 
        color, marker = plotKeys[mol] 
        datafile = '%s/Analysis/boltzmann/O1G_size_%s.weighted.out'%(molec,dist) 
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
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        elif color == 'k' : 
            avgkAccum = np.append(avgkAccum,avg) 
            ratekAccum = np.append(ratekAccum,rate) 
    
        ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
        ax.annotate(mol, (rate, avg) ) 
    
        print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
        x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='b')
    except ValueError : 
        print "No data found for blue" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
        x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='r')
    except ValueError : 
        print "No data found for red" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
        x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='k')
    except ValueError : 
        print "No data found for black" 
    
    ax.legend(loc=4) 
    ax.set_xscale('log') 
    
    fig.savefig('figures/rate_v_O1G_scount_%s.png'%dist,format='png') 
    plt.close() 

for dist in ['245', '3', '4', '5', '6'] : 
    fig, ax = plt.subplots(1,1) 
    fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
    fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
    fig.text(0.03,0.5, r"Num waters within %s A of Q61 s.c. and 5 A of O1G"%(dist) , ha='center', va='center',rotation='vertical') 
    
    avgbAccum, ratebAccum = [], [] 
    avgrAccum, raterAccum = [], [] 
    avgkAccum, ratekAccum = [], [] 
    for mol in molecList : 
        molec = 'RasRalC18CNC_Q61%s'%mol
        name = "Q61%s"%mol 
        color, marker = plotKeys[mol] 
        datafile = '%s/Analysis/boltzmann/SCD_size_%s.weighted.out'%(molec,dist) 
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
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        elif color == 'k' : 
            avgkAccum = np.append(avgkAccum,avg) 
            ratekAccum = np.append(ratekAccum,rate) 
    
        ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
        ax.annotate(mol, (rate, avg) ) 
    
        print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
        x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='b')
    except ValueError : 
        print "No data found for blue" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
        x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='r')
    except ValueError : 
        print "No data found for red" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
        x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='k')
    except ValueError : 
        print "No data found for black" 
    
    ax.legend(loc=4) 
    ax.set_xscale('log') 
    
    fig.savefig('figures/rate_v_SCD_scount_%s.png'%dist,format='png') 
    plt.close() 

sys.exit() 

for dist in ['245', '3', '4', '5'] : 
    fig, ax = plt.subplots(1,1) 
    fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
    fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
    fig.text(0.03,0.5, r"Num waters within %s A of Q61 s.c. and %s A of O1G or O2G"%(dist, dist) , ha='center', va='center',rotation='vertical') 
    
    avgbAccum, ratebAccum = [], [] 
    avgrAccum, raterAccum = [], [] 
    avgkAccum, ratekAccum = [], [] 
    for mol in molecList : 
        molec = 'RasRalC18CNC_Q61%s'%mol
        name = "Q61%s"%mol 
        color, marker = plotKeys[mol] 
        datafile = '%s/Analysis/boltzmann/size_polar_%s.weighted.out'%(molec,dist) 
        #print datafile
    
        rate, error = nameToExpRate[name] 
        
        try : 
            with open(datafile) as f : 
                lines = f.readlines() 
                avg = float(lines[1].split()[-2]) 
                std = float(lines[1].split()[-1])
        except IOError : 
            print "No file found for %s"%(datafile)  
            avg, std = 0, 0 
            #continue     
        except : 
            print "Error importing data from file %s"%(datafile)
            continue 
    
        avg /= 3 ##atoms -> water molecules
        std /= 3 ##atoms -> water molecules 
    
        if color == 'b' : 
            avgbAccum = np.append(avgbAccum,avg) 
            ratebAccum = np.append(ratebAccum,rate) 
        elif color == 'r' : 
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        elif color == 'k' : 
            avgkAccum = np.append(avgkAccum,avg) 
            ratekAccum = np.append(ratekAccum,rate) 
    
        ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
        ax.annotate(mol, (rate, avg) ) 
    
        print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
        x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='b')
    except ValueError : 
        print "No data found for blue" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
        x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='r')
    except ValueError : 
        print "No data found for red" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
        x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='k')
    except ValueError : 
        print "No data found for black" 
    
    ax.legend(loc=4) 
    ax.set_xscale('log') 
    
    fig.savefig('figures/rate_v_polar_scount_%s.png'%dist,format='png') 
    plt.close() 

for dist in ['245', '3', '4', '5'] : 
    fig, ax = plt.subplots(1,1) 
    fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
    fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
    fig.text(0.03,0.5, r"Num waters within %s A of Q61 s.c. and %s A of O1G or O2G"%(dist, dist) , ha='center', va='center',rotation='vertical') 
    
    avgbAccum, ratebAccum = [], [] 
    avgrAccum, raterAccum = [], [] 
    avgkAccum, ratekAccum = [], [] 
    for mol in molecList : 
        molec = 'RasRalC18CNC_Q61%s'%mol
        name = "Q61%s"%mol 
        color, marker = plotKeys[mol] 
        datafile = '%s/Analysis/boltzmann/both_size_%s.weighted.out'%(molec,dist) 
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
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        elif color == 'k' : 
            avgkAccum = np.append(avgkAccum,avg) 
            ratekAccum = np.append(ratekAccum,rate) 
    
        ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
        ax.annotate(mol, (rate, avg) ) 
    
        print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
        x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='b')
    except ValueError : 
        print "No data found for blue" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
        x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='r')
    except ValueError : 
        print "No data found for red" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
        x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='k')
    except ValueError : 
        print "No data found for black" 
    
    ax.legend(loc=4) 
    ax.set_xscale('log') 
    
    fig.savefig('figures/rate_v_both_scount_%s.png'%dist,format='png') 
    plt.close() 

for dist in ['245', '3', '4', '5'] : 
    fig, ax = plt.subplots(1,1) 
    fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
    fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
    fig.text(0.03,0.5, r"Num waters within %s A of Q61 s.c. and %s A of O1G or O2G"%(dist, dist) , ha='center', va='center',rotation='vertical') 
    
    avgbAccum, ratebAccum = [], [] 
    avgrAccum, raterAccum = [], [] 
    avgkAccum, ratekAccum = [], [] 
    for mol in molecList : 
        molec = 'RasRalC18CNC_Q61%s'%mol
        name = "Q61%s"%mol 
        color, marker = plotKeys[mol] 
        datafile = '%s/Analysis/boltzmann/size_polar_%s.weighted.out'%(molec,dist) 
        #print datafile
    
        rate, error = nameToExpRate[name] 
        
        try : 
            with open(datafile) as f : 
                lines = f.readlines() 
                avg = float(lines[1].split()[-2]) 
                std = float(lines[1].split()[-1])
        except IOError : 
            print "No file found for %s"%(datafile)  
            avg, std = 0, 0 
            #continue     
        except : 
            print "Error importing data from file %s"%(datafile)
            continue 
    
        avg /= 3 ##atoms -> water molecules
        std /= 3 ##atoms -> water molecules 
    
        if color == 'b' : 
            avgbAccum = np.append(avgbAccum,avg) 
            ratebAccum = np.append(ratebAccum,rate) 
        elif color == 'r' : 
            avgrAccum = np.append(avgrAccum,avg) 
            raterAccum = np.append(raterAccum,rate) 
        elif color == 'k' : 
            avgkAccum = np.append(avgkAccum,avg) 
            ratekAccum = np.append(ratekAccum,rate) 
    
        ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
        ax.annotate(mol, (rate, avg) ) 
    
        print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratebAccum), avgbAccum) 
        x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='b')
    except ValueError : 
        print "No data found for blue" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(raterAccum), avgrAccum) 
        x = np.linspace(np.min(raterAccum), np.max(raterAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='r')
    except ValueError : 
        print "No data found for red" 
    
    try : 
        slope, intercept, r_value, p_value, std_error = linregress(np.log(ratekAccum), avgkAccum) 
        x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100) 
        y = slope * np.log(x) + intercept 
        ax.plot(x, y, label = "r = %.3f"%r_value,color='k')
    except ValueError : 
        print "No data found for black" 
    
    ax.legend(loc=4) 
    ax.set_xscale('log') 
    
    fig.savefig('figures/rate_v_polar_scount_%s.png'%dist,format='png') 
    plt.close() 


sys.exit() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Num waters within 4A of G60 backbone and 3A of O1G", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/size2.weighted.out'%molec
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

    #avg *= 10**3 ##nm^3 -> A^3
    #std *= 10**3 ##nm^3 -> A^3

    if color == 'b' : 
        avgbAccum = np.append(avgbAccum,avg) 
        ratebAccum = np.append(ratebAccum,rate) 
    elif color == 'r' : 
        avgrAccum = np.append(avgrAccum,avg) 
        raterAccum = np.append(raterAccum,rate) 
    elif color == 'k' : 
        avgkAccum = np.append(avgkAccum,avg) 
        ratekAccum = np.append(ratekAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
    ax.annotate(mol, (rate, avg) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

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

ax.legend(loc=4) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_scount_G60_linear.png',format='png') 
plt.close() 

fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (min$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Num waters within 4A of G60 backbone and 3A of O1G", ha='center', va='center',rotation='vertical') 

avgbAccum, ratebAccum = [], [] 
avgrAccum, raterAccum = [], [] 
avgkAccum, ratekAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    color, marker = plotKeys[mol] 
    datafile = '%s/Analysis/boltzmann/size2.weighted.out'%molec
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

    #avg *= 10**3 ##nm^3 -> A^3
    #std *= 10**3 ##nm^3 -> A^3

    if color == 'b' : 
        avgbAccum = np.append(avgbAccum,avg) 
        ratebAccum = np.append(ratebAccum,rate) 
    elif color == 'r' : 
        avgrAccum = np.append(avgrAccum,avg) 
        raterAccum = np.append(raterAccum,rate) 
    elif color == 'k' : 
        avgkAccum = np.append(avgkAccum,avg) 
        ratekAccum = np.append(ratekAccum,rate) 

    ax.errorbar(rate, avg,xerr=error, marker=marker,color=color) 
    ax.annotate(mol, (rate, avg) ) 

    print "Q61%s\t%3.2f\t%.2f"%(mol,avg,std) 

popt, pcov = curve_fit(logFit, ratebAccum, avgbAccum, p0=(1, -1, ))
p1, p2, = popt[:]
x = np.linspace(np.min(ratebAccum), np.max(ratebAccum),100)
y = logFit(x, p1, p2)
ax.plot(x,y,color='b', label = "E = %.3f"%np.sum(np.sqrt(np.diag(pcov) ))  )

popt, pcov = curve_fit(logFit, raterAccum, avgrAccum, p0=(1, -1, ))
p1, p2, = popt[:]
x = np.linspace(np.min(raterAccum), np.max(raterAccum),100)
y = logFit(x, p1, p2)
ax.plot(x,y,color='r', label = "E = %.3f"%np.sum(np.sqrt(np.diag(pcov) ))  )

popt, pcov = curve_fit(logFit, ratekAccum, avgkAccum, p0=(1, -1, ))
p1, p2, = popt[:]
x = np.linspace(np.min(ratekAccum), np.max(ratekAccum),100)
y = logFit(x, p1, p2)
ax.plot(x,y,color='k', label = "E = %.3f"%np.sum(np.sqrt(np.diag(pcov) ))  )

ax.legend(loc=4) 
ax.set_xscale('log') 

fig.savefig('figures/rate_v_scount_G60_nonlinear.png',format='png') 
plt.close() 
