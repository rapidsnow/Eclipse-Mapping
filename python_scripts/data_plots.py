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

def plot_brights_stripes_only(path, bright_file, nstripes, output):
    '''
    Take in a set of brightness values and output a brightness map.
    This snippet of code is adapted from Leslie's IDL code
    '''
    
    #Read in file
    temp = []
    for line in open(bright_file):                           
        temp.append(float(line))
    brights = np.array(temp)
    
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
    
    plt.imsave(path + "/" + output, img, cmap='hot')
    plt.imshow(img, cmap='hot')
    plt.clim(.9,1.25)
    plt.title("Surface Brightness Map")
    #Create the plot
    bmap = Basemap(projection='moll', lon_0 = 0)
    #plt.savefig(path + "/" + output)
    bmap.warpimage(path + "/" + output)
    #Place the plot in the appropriate place
    plt.savefig(path + "/" + output)
    plt.close()

def plot_lc_and_brights(path, bright_file, chi_file, model_curve, data_curve, bb, Rp, nstripes, nboxes, outfile, number, trans_info):
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

    try:
        information = csv.reader(open(trans_info))
    except:
        information = []
    transit_count = 1
    
    for row in information:
        trans_time = []
        trans_flux = []
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
        nx = nboxes * 50
        ny = nboxes * 50
        bright_stripe = np.array(brights[nboxes:nboxes+nstripes])
        bright_box = np.array(brights[0:nboxes])
        img = np.zeros((nx,ny))
        idx1_stripe = np.arange(0, nstripes, 1.)/nstripes * nx
        idx2_stripe = (np.arange(0, nstripes, 1.) + 1)/nstripes * nx - 1
        idx1_box = np.arange(0, nboxes, 1.)/nboxes * nx
        idx2_box = (np.arange(0, nboxes, 1.) + 1)/nboxes * nx - 1
        idy1_lat = int(lat1 * ny/2 + ny/2 + 0.5)
        idy2_lat = int(lat2 * ny/2 + ny/2 + 0.5) - 1
    
        #Populate the pixel values
        for jj in range(nstripes):
            idx1 = idx1_stripe[jj]
            idx2 = idx2_stripe[jj]
            idy1 = 0
            idy2 = idy1_lat-1
            img[idy1:idy2, idx1:idx2] = bright_stripe[jj]
            idy1 = idy2_lat+1
            idy2 = ny-1
            img[idy1:idy2, idx1:idx2] = bright_stripe[jj]
        for jj in range(nboxes):
            idx1 = idx1_box[jj]
            idx2 = idx2_box[jj]
            idy1 = idy1_lat
            idy2 = idy2_lat
            img[idy1:idy2, idx1:idx2] = bright_box[jj]
    
        ax = plt.axes([.15, .6, .32, .32])
        plt.imsave(path + "/temp.png", img, cmap='hot', vmin=0.85, vmax=1.15)
        plt.imshow(img, cmap='hot')    
        plt.clim(0.85,1.15)
        plt.colorbar(shrink=0.5)    

        #Create the plot
        bmap = Basemap(projection='moll', lon_0 = 0)
        bmap.warpimage(path + "/temp.png")

        #Create the plots of the transits
        for line in open(path + "/" + "input_files/binned_%d.out" % number):
            new_array = line.split(' ')
            trans_time.append(float(new_array[0]))
            trans_flux.append(float(new_array[1]))
        minimum = min(trans_flux[int(row[0]):int(row[1])])
        maximum = max(trans_flux[int(row[0]):int(row[1])])
        diff = maximum - minimum
        trans_ax = plt.axes([.58,.6,.27,.27])
        plt.plot(model_time[int(row[0]):int(row[1])], model_flux[int(row[0]):int(row[1])], c='black', label="Model")
        plt.scatter(trans_time[int(row[0]):int(row[1])], trans_flux[int(row[0]):int(row[1])], c='red', marker='+', label="Data")
        plt.xlabel("Orbital Phase")
        plt.ylabel("Relative Flux")
        plt.title("Transit %d" % transit_count)
        plt.ylim(ymax=maximum + (diff * .08))
        plt.ylim(ymin=minimum - (diff * .08))
        #plt.setp(trans_ax)
        fig = plt.gcf()
        fig.set_size_inches(18.5,11.5)
        plt.savefig(path + "/" + outfile + "_t%02d" % transit_count + ".png", dpi=120)
        print "Window:  %d Trans: %d" % (number, transit_count - 1)
        transit_count += 1
        del trans_time[:]
        del trans_flux[:]
        plt.close()
        
def regions_over_time(path, outfile, nfiles, nboxes, nstripes):
    nsb = nboxes + nstripes
    box_per_stripe = nboxes/nstripes
    #brights = [[] for ii in range(nsb)]
    brights_time = np.arange(float(nsb * nfiles)).reshape(nsb, nfiles)
    window_number = np.arange(nfiles) #Array to match number of files so that we can keep track of windows
    fig = plt.gcf()
    fig.set_size_inches(18.5,11.5)

    for ii in xrange(nfiles):
        jj = 0
        for line in open(path + "brights_%d.out" % ii):
            line = line.rstrip('\r\n')
            brights_time[jj][ii] = float(line)
            jj += 1

    #Okay, I now have all of the brightnesses for all files in a matrix called brights.
    #The indices of brights are brights[region_number][file_number] where region number goes 1 -> (nboxes - 1), nboxes -> (nsb - 1)

    #Plotting utility arrays (to cycle through)
    color_array = ('magenta', 'blue', 'red', 'green', 'cyan', 'orange')
    lstyles = ['-.', '--', ':']
    y_ticks = np.arange(1.0, 1 + (nstripes+1) * .3, .3)
    y_labels = np.ones(nstripes)

    #Plot the stripe and box groups together, scaling up each time so that they don't overlap
    scale = 0
    for ii in range(0, nboxes, box_per_stripe):
        plt.plot(window_number, np.ones(nfiles, dtype='float') + scale, c='black')
        for jj in range(box_per_stripe):
            plt.plot(window_number, brights_time[ii + jj] + scale, c=color_array[((ii+1)/box_per_stripe) % len(color_array)], linestyle=lstyles[jj], linewidth=2.0, marker='o')
        plt.plot(window_number, brights_time[ii/nboxes + nboxes] + scale, c=color_array[((ii+1)/box_per_stripe) % len(color_array)], linestyle='-', marker='o')    
        scale += .3 #Addition factor so that subsequent graphs don't overlap
        

    plt.xlabel("Window Number")
    plt.ylabel("Relative Brightness")
    plt.ylim(ymax=1.15 + (nstripes - 1) * .3) #Min and Max are equal to the scale factor on the brightness maps
    plt.ylim(ymin=.85)
    plt.xticks(window_number) #So that we don't get half value tickmarks (make no sense)
    plt.yticks(y_ticks, y_labels) #No misleading values on scaled y axis
    plt.title(path)

    plt.savefig(path + "/" + outfile + ".png", dpi=120)
    plt.close()
    
    #Plot all of the stripes together for all datasets
    plt.plot(window_number, np.ones(nfiles, dtype='float'), c='black')
    for ii in range(nboxes, nsb):
        plt.plot(window_number, brights_time[ii], c=color_array[ii % len(color_array)], marker='o')
    plt.xlabel("Window Number")
    plt.ylabel("Relative Brightness")
    plt.xticks(window_number)
    plt.title(path)
    plt.ylim(ymax= max([max(bright) for bright in brights_time]))
    plt.ylim(ymin= min([min(bright) for bright in brights_time]))
    
    plt.savefig(path + "/stripes_over_time.png", dpi=120)
    plt.close()
    
    #Plot the boxes together for each group of boxes
    lstyles[0] = '-'
    for ii in range(0, nboxes, box_per_stripe):
        plt.plot(window_number, np.ones(nfiles, dtype='float'), c='black')
        for jj in range(box_per_stripe):
            plt.plot(window_number, brights_time[ii + jj], marker='o', linestyle=lstyles[jj], color=color_array[((ii % (box_per_stripe+1)) + nboxes) % len(color_array)], label="Box %d" % jj)
        plt.xlabel("Window Number")
        plt.ylabel("Relative Brightness")
        plt.legend()
        plt.xticks(window_number)
        plt.title(path)
        plt.ylim(ymax=1.15)
        plt.ylim(ymin=.85)
        
        plt.savefig(path + "/box_%d.png" % ii, dpi=120)
        plt.close()

if __name__ == '__main__':
    for ii in range(0,64):
        plot_lc_and_brights(".", "input_files/brights_%d.out" % ii, "input_files/chi_%d.out" % ii, "input_files/model_%d.out" % ii, "input_files/binned.out", .268, .13031, 7, 14, "images/brights_lc_test_%02d" % ii, ii, "input_files/trans_%d.out" % ii)
