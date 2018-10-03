import glob 
import numpy as np 
import math
import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm
import os
import matplotlib as mpl 
from scipy.stats import linregress
from matplotlib import rc_file
from sys import exit

rcFile = 'rc_files/paper.rc'  

nameToColorKeys = {
        "F145X":'#A93226',
        "M218X":'r',
        "F165X":'#E67E22',
        "F114X":'y',
        "N198X":'#2ECC71',
        "Y143X":'g',
        "N212X":'c',
        "D117X":'b',
        "Y92X":'m'
        }

molecList = [
        "RasRalC18CNC_Q61D", 
        "RasRalC18CNC_Q61E", 
        "RasRalC18CNC_Q61F", 
        "RasRalC18CNC_Q61H", 
        "RasRalC18CNC_Q61I", 
        "RasRalC18CNC_Q61K", 
        "RasRalC18CNC_Q61L", 
        "RasRalC18CNC_Q61M", 
        "RasRalC18CNC_Q61N", 
        "RasRalC18CNC_Q61Q", 
        "RasRalC18CNC_Q61R", 
        "RasRalC18CNC_Q61S", 
        "RasRalC18CNC_Q61T", 
        "RasRalC18CNC_Q61V", 
        "RasRalC18CNC_Q61W",
        "RasRalC18CNC_Q61Y" 
        ]

figCols=4
figRows=4

rc_file(rcFile) 

if not os.path.isdir('figures') : 
    os.mkdir('figures') 


#for angle in np.arange(0,359,30) : 
#    print angle 
#    fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
#    fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
#    fig.text(0.5,0.02, r"$d_{\rm{NH}}$ ($\rm{\AA}$)", ha='center', va='center') 
#    fig.text(0.01,0.5, r"$\theta_1$ (deg)", ha='left', va='center',rotation='vertical') 
#    
#    xmin,xmax = 1.45, 2.45
#    ymin,ymax = 99,180 
#    
#    for index,molec in enumerate(molecList) : 
#        print "\t%s"%molec
#        file="%s/Analysis/hbond/geometry.%s.xvg"%(molec,angle) 
#        try : 
#            data = np.genfromtxt(file,skip_header=23) 
#        except :
#            print "\t\tNo file found" 
#            continue 
#        try :
#            data[:,0] = data[:,0] / 1000 * 4
#        except IndexError : 
#            print "%s is empty. Skipping"%file      
#            data = np.empty((0,4),int)  
#            data = np.append(data,[[0,0,0,0]],axis=0) 
#    
#        for i in range(len(data[:,3]))  : 
#            if data[i,3] > 180 : 
#                data[i,3] -= 180 
#        ax = axarr[index/figCols,index%figCols]
#    
#        xbins, ybins = np.arange(xmin, xmax, 0.01), np.arange(ymin,ymax,1)
#        z, x, y = np.histogram2d(data[:,2],data[:,3],[xbins,ybins])
#    
#        im = ax.pcolor(x,y,z.T, cmap='plasma', vmin = 0, vmax = 2) 
#    
#        ## Overlay contour lines 
#        x = data[:,2] ; y = data[:,3]
#        counts,  xbins,  ybins  = np.histogram2d(x, y,range=((xmin,xmax),(ymin,ymax))) #,bins=(64,64)) 
#        #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
#    
#    #    ax.imshow(heatmap) 
#    #    plt.close() 
#    
#        ax.contour(counts.transpose(),extent=[xmin,xmax,ymin,ymax],
#           #ybins.min(),ybins.max()],linewidths=1,colors='black',
#           linewidths=1,colors='black',
#               linestyles='solid',levels=np.arange(-1,1500,125))
#    
#        ax.set_title(molec.split('_')[1],color='k') 
#        ax.set_ylim([ymin,179])
#        ax.set_xlim([xmin,xmax]) 
#    
#    fig.subplots_adjust(right=0.88) 
#    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7]) 
#    fig.text(0.98,0.5, r"Counts per bin", ha='center', va='center',rotation='vertical') 
#    fig.colorbar(im, cax=cbar_ax) 
#    
#    
#    fig.savefig('figures/Geometries_heat_%s.png'%angle,format='png') 
#    plt.close() 



def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))


###################################
fig, axarr = plt.subplots(figRows,figCols,sharex='col',sharey='row') 
fig.subplots_adjust(wspace=0.10,hspace=0.35,left=0.12,right=0.98,top=0.93,bottom=0.1) 
fig.text(0.5,0.02, r"Angle (deg)",ha='center', va='center') 
fig.text(0.01,0.5, r"Probability", ha='left', va='center',rotation='vertical') 

for index, molec in enumerate(molecList) : 
    #print molec 
    
    xmin,xmax = 1.45, 2.45
    ymin,ymax = 99,180 

    ax = axarr[index/figCols,index%figCols]
    

    probFile="%s/Analysis/wham/%s.output.prob"%(molec,molec) 
    probs = np.genfromtxt(probFile ) 

#    totHBonds, totAngle = 0,0
    totProb, totAngle = 0,0
    accumAngle,accumDist,accumAngle2, accumProbs = [], [], [], []
    for index, angle in enumerate(np.arange(0,359,30)) : 
        binFile="%s/Analysis/wham/%s.output.%s.bin"%(molec,molec,index)
        bins = np.genfromtxt(binFile ) 

        dataFile="%s/Analysis/hbond/geometry.%s.xvg"%(molec,angle) 
        try : 
            data = np.genfromtxt(dataFile,skip_header=23) 
        except :
            print "\t\tNo file found" 
            continue 

        for i in range(len(data)) : 
            frame, dist, angle,angle2 = int(data[i,0]),data[i,2],data[i,3],data[i,4]
            whamBin = int(bins[frame]) 
            prob = probs[whamBin]

            accumAngle.append(angle) 
            accumAngle2.append(angle2) 
            accumDist.append(dist) 
            accumProbs.append(prob) 

            totProb += prob
            totAngle += angle*prob

        #file="%s/Analysis/hbond/theta1.%s.poly"%(molec,angle) 
        #try : 
        #    data = np.genfromtxt(file) 
        #except :
        #    #print "\t\tNo file found" 
        #    continue 
    
        #file="%s/Analysis/hbond/theta1.%s.his"%(molec,angle) 
        #try : 
        #    data = np.genfromtxt(file) 
        #except :
        #    #print "\t\tNo file found" 
        #    continue 
        ##ax.scatter(data[:,0], data[:,1],s=0.10) 
    
    #print "%s\t%0.1f\t%0.1f\t%0.1f"%(molec, np.average(totData), np.std(totData),totAngle/totHBonds ) 
    ang, std = weighted_avg_and_std(accumAngle, accumProbs) 
    dist, stdD = weighted_avg_and_std(accumDist, accumProbs) 
    ang2, std2 = weighted_avg_and_std(accumAngle2, accumProbs) 
    print "%s\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f"%(molec, ang,std,dist,stdD, ang2, std2) 
    #size, bins = np.histogram(totData[:,3]) 
    #ax.scatter(size, bins[:-1]) 
    #ax.scatter(totData[:,0],totData[:,3],s=0.1) 



    #del totData

    ax.set_title(molec.split('_')[1],color='k') 
    #ax.set_ylim([ymin,179])
    #ax.set_xlim([xmin,xmax]) 

    #print "%s\t%0.1f" %(molec, totAngle / totHBonds) 
    
fig.savefig('figures/Geometries_angles.png',format='png') 
plt.close() 

