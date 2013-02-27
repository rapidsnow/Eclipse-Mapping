'''
Created on Jun 23, 2012

@author: rapidsnow
'''

# -------------------
# Third party imports
# -------------------

import matplotlib.pyplot as plt
import numpy as np
import csv
from mpl_toolkits.basemap import Basemap

# --------------------------
# Define Classes and Methods
# --------------------------
def plot_lc_and_brights(path, bright_file, model_curve, data_curve, bb, Rp, nstripes, outfile, number):
    #Read in the values for model and data light curves
    model_time = []
    model_flux = []
    data_time_on = []
    data_flux_on = []
    data_time_off = []
    data_flux_off = []
    
    #Read in file
    temp = []
    for line in open(bright_file):                           
        temp.append(float(line))
    brights = np.array(temp)
   
    #Define latitudes of boxes
    lat1 = bb - Rp
    lat2 = bb + Rp

    for line in open(model_curve):                            
        new_array = line.split(' ')
        model_time.append(float(new_array[0]))
        model_flux.append(float(new_array[1]))
        
    for line in open(data_curve):                            
        new_array = line.split(' ')
        if float(new_array[0]) >= model_time[0] and float(new_array[0]) <= model_time[len(model_time) - 1]:
            data_time_on.append(float(new_array[0]))
            data_flux_on.append(float(new_array[1]))
        else:
            data_time_off.append(float(new_array[0]))
            data_flux_off.append(float(new_array[1]))
    
    #Create the plot and label it
    plt.plot(model_time, model_flux, c='black', label="Model")
    plt.scatter(data_time_on, data_flux_on, c='red', marker='+', label="Brightness Map Region")
    plt.scatter(data_time_off, data_flux_off, c='green', marker = '+', label="Data")
    minimum = min(data_flux_on + data_flux_off)
    maximum = max(data_flux_on + data_flux_off)
    diff = maximum - minimum
    plt.ylim(ymax=maximum + (diff * .85))
    plt.ylim(ymin=minimum - (diff * .1))
    plt.xlabel("Orbital Phase")
    plt.ylabel("Relative Flux")
    plt.title("Model Light Curve Vs. Data - Window %d" % number)

       #Set up arrays for brightnesses and positions
    nx = nstripes * 50
    ny = nstripes * 50
    bright_stripe = np.array(brights)
    img = np.zeros((nx,ny))
    idx1_stripe = np.arange(0, nstripes, 1.)/nstripes * nx
    idx2_stripe = (np.arange(0, nstripes, 1.) + 1)/nstripes * nx - 1

    #Populate the pixel values
    for jj in range(nstripes):
        idx1 = idx1_stripe[jj]
        idx2 = idx2_stripe[jj]
        idy1 = 0
        idy2 = ny-1
        img[idy1:idy2, idx1:idx2] = bright_stripe[jj]

    ax = plt.axes([.15, .6, .32, .32])
    plt.imsave(path + "/temp.png", img, cmap='hot', vmin=0.95, vmax=1.05)
    plt.imshow(img, cmap='hot')    
    plt.clim(0.95,1.05)
    plt.colorbar(shrink=0.5)    

    #Create the plot
    bmap = Basemap(projection='moll', lon_0 = 0)
    bmap.warpimage(path + "/temp.png")

    #Create the plots of the transits
    fig = plt.gcf()
    fig.set_size_inches(18.5,11.5)
    plt.savefig(path + "/" + outfile + ".png", dpi=120)
    plt.close()
if __name__ == '__main__':
    for ii in range(78):
        plot_lc_and_brights(".", "input_files/brights_%d.out" % ii, "input_files/model_%d.out" % ii, "input_files/binned.out", .268, .13031, 12, "images/brights_lc_test_%02d" % ii, ii)
