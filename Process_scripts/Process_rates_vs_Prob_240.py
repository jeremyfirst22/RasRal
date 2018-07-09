import glob 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from os import sys
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines 
from scipy.stats import linregress
from matplotlib import rc_file

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
colorKeys = {
        'Q':'b',
        'D':'r',
        'E':'r',
        'H':'r',
        'K':'r',
        'R':'r',
        'N':'b',
        'S':'b',
        'T':'b',
        'W':'b',
        'Y':'b',
        'A':'k',
        'F':'k',
        'G':'k',
        'I':'k',
        'L':'k',
        'M':'k',
        'V':'k'}



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


fig, ax = plt.subplots(1,1) 
fig.subplots_adjust(wspace=0.1,hspace=0.35,left=0.15,right=0.95) 
fig.text(0.5,0.04, r"Rate (s$^{-1}$)", ha='center', va='center') 
fig.text(0.03,0.5, r"Probability of angle -120", ha='center', va='center',rotation='vertical') 

probAccum, rateAccum = [], [] 
for mol in molecList : 
    molec = 'RasRalC18CNC_Q61%s'%mol
    name = "Q61%s"%mol 
    datafile = '%s/Analysis/wham/%s.output.prob'%(molec,molec) 
    #print datafile
    try : 
        probs = np.genfromtxt(datafile) 
        prob = np.max(probs[100:140]) 
    except IOError : 
        print "No file found for %s"%(datafile)  
    except : 
        print "Error importing data from file %s"%(datafile)

    rate, error = nameToExpRate[name] 

    if mol == "M" or mol == "F" or mol == "L" or mol == "V" or mol == "I" : 
        probAccum = np.append(probAccum,prob) 
        rateAccum = np.append(rateAccum,rate) 
        marker = '^'
    else : 
        probAccum = np.append(probAccum,prob) 
        rateAccum = np.append(rateAccum,rate) 
        marker = 'o'

    ax.errorbar(rate, prob,xerr=error, marker=marker,color=colorKeys[mol],capsize=3) 
    #ax.scatter(rate, avg) 
    ax.annotate(mol, (rate, prob+0.001) ) 

    #index +=1

slope, intercept, r_value, p_value, std_error = linregress(rateAccum, probAccum) 
x = np.linspace(np.min(rateAccum), np.max(rateAccum),100) 
y = slope * x + intercept 

ax.plot(x, y, label = "r = %.3f"%r_value,color='k')

ax.legend(loc=1) 
#ax.set_xscale('log') 

fig.savefig('figures/rate_v_prob_chi1_240.png',format='png') 
plt.close() 

