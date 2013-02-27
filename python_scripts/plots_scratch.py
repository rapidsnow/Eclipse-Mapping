import matplotlib.pyplot as plt
import numpy as np
import math
import csv
from mpl_toolkits.basemap import Basemap
 
class Star(object):
    def __init__(self, b, Rp, nStripes, nBoxes):
        self.b = b
        self.Rp = Rp
        self.nStripes = nStripes
        self.nBoxes = nBoxes
        self.nsb = nBoxes + nStripes
        self.boxPerStripe = nBoxes/nStripes
        self.lat1 = b - Rp
        self.lat2 = b + Rp

class Region(object):
    def __init__(self, nFiles, goalVal=1):
        self.brights = [0] * nFiles
        self.goalVal = goalVal
    def add_bright_value(self, val, window):
        self.brights[window] = val
    def get_bright_at_time(self, window):
        try:
            return self.brights[window]
        except IndexError:
            raise IndexError('Window not yet defined')
    def get_brights(self):
        return np.array(self.brights[:])
    def rms(self):
        sumVals = sum([(b - self.goalVal) ** 2 for b in self.brights])
        rms = math.sqrt((1.0/len(self.brights)) * sum([(b - self.goalVal) ** 2 for b in self.brights]))
        return rms
    def get_mean(self):
        return sum(self.brights)/len(self.brights)
    def get_max(self):
        return max(self.brights)
    def get_min(self):
        return min(self.brights)
    def get_window_number(self):
        return len(self.brights)
        
class Light_Curve(object):
    '''
    Takes a filename and creates a lightcurve object
    '''
    def __init__(self, filename):
        temp = []
        self.time = []
        self.flux = []
        self.error = []
        
        for line in open(filename):
            temp = line.split(' ')
            self.time.append(float(temp[0]))
            self.flux.append(float(temp[1]))
            self.error.append(float(temp[2]))
        self.time = np.array(self.time)
        self.flux = np.array(self.flux)
        self.error = np.array(self.error)
    def get_flux(self):
        return np.array(self.flux[:])
    def get_time(self):
        return np.array(self.time[:])
    def get_error(self):
        return np.array(self.error[:])
    def flux_max(self):
        return max(self.flux)
    def flux_min(self):
        return min(self.flux)
    
#class Data_Container(object):
#    def __init__ (self, total_binned_file):#, model_file, data_file, trans_file, bright_file
#        self.total_binned_file = Light_Curve(total_binned_file)
#        self.models = []
#        self.dataCurves = []
        
def make_brightness_structures(goalFile, fileList):
    '''
    Allocates a set of brightness values for each region and returns a list of Region objects
    '''
    #Creates a list of regions with goals defined by an external file
    regionList = [Region(len(fileList), float(goal)) for goal in open(goalFile)]
    
    #Adds brightness values to the regionList
    for i in range(len(fileList)):
        valList = [float(line) for line in open(fileList[i])]
        for j in range(len(regionList)):
            regionList[j].add_bright_value(valList[j], i)
    return regionList

def regions_over_time(path, outfile, star, regionList):
    '''
    This method is old and gross :(
    '''
    nWindows = regionList[0].get_window_number()
    window_number = np.arange(nWindows) #Array to match number of files so that we can keep track of windows
    fig = plt.gcf()
    fig.set_size_inches(18.5,11.5)
    scaleFactor = .3
    yMin = .85
    yMax = 1.15
    

    #Plotting utility arrays (to cycle through)
    color_array = ['magenta', 'blue', 'red', 'green', 'cyan', 'orange']
    lstyles = ['-.', '--', ':']
    y_ticks = np.arange(1.0, 1 + (star.nStripes+1) * scaleFactor, scaleFactor)
    y_labels = np.ones(star.nStripes)

#Plot the stripe and box groups together, scaling up each time so that they don't overlap    
    scale = 0
    for ii in range(0, star.nBoxes, star.boxPerStripe):
        color = color_array[((ii+1)/star.boxPerStripe) % len(color_array)]
        plt.plot(window_number, np.ones(nWindows, dtype='float') + scale, c='black')
        for jj in range(star.boxPerStripe):
            box = regionList[ii + jj].get_brights()
            plt.plot(window_number, box + scale, c=color, linestyle=lstyles[jj], linewidth=2.0, marker='o')
        plt.plot(window_number, regionList[(ii+1)/star.nBoxes + star.nBoxes].get_brights() + scale, c=color, linestyle='-', marker='o')    
        scale += scaleFactor #Addition factor so that subsequent graphs don't overlap
        

    plt.xlabel("Window Number")
    plt.ylabel("Relative Brightness")
    plt.ylim(ymax=yMax + (star.nStripes - 1) * scaleFactor) #Min and Max are equal to the scale factor on the brightness maps
    plt.ylim(ymin=yMin)
    plt.xticks(window_number) #So that we don't get half value tickmarks (make no sense)
    plt.yticks(y_ticks, y_labels) #No misleading values on scaled y axis
    plt.title(path)

    plt.savefig(path + "/" + outfile + ".png", dpi=120)
    plt.close()
    
#Plot all of the stripes together for all datasets
    plt.plot(window_number, np.ones(nWindows, dtype='float'), c='black')
    for ii in range(star.nBoxes, star.nsb):
        plt.plot(window_number, regionList[ii].get_brights(), c=color_array[ii % len(color_array)], marker='o')
    plt.xlabel("Window Number")
    plt.ylabel("Relative Brightness")
    plt.xticks(window_number)
    plt.title(path)
    plt.ylim(ymax=1.1)
    plt.ylim(ymin=0.9)
    #plt.ylim(ymax= max([bright.get_max() for bright in regionList]))
    #plt.ylim(ymin= min([bright.get_min() for bright in regionList]))
    
    plt.savefig(path + "/stripes_over_time.png", dpi=120)
    plt.close()
    
#Plot the boxes together for each group of boxes
    lstyles[0] = '-'
    for ii in range(0, star.nBoxes, star.boxPerStripe):
        plt.plot(window_number, np.ones(nWindows, dtype='float'), c='black', linewidth=2.0)
        color = color_array[((ii+1)/star.boxPerStripe) % len(color_array)]
        for jj in range(star.boxPerStripe):
            average = regionList[ii + jj].get_mean()
            box = regionList[ii + jj].get_brights()
            
            plt.plot(window_number, box, marker='o', linestyle=lstyles[jj], c=color, label="Box %d" % (jj + ii))
            plt.plot(window_number, np.ones(nWindows) * average, c='black', linestyle=lstyles[jj])
        plt.xlabel("Window Number")
        plt.ylabel("Relative Brightness")
        plt.legend()
        plt.xticks(window_number)
        plt.title(path)
        plt.ylim(ymax=yMax)
        plt.ylim(ymin=yMin)
        
        plt.savefig(path + "/box_%d.png" % ii, dpi=120)
        plt.close()
        
def rms_vs_region(path, outfile, star, regionList):
    plt.scatter(range(star.nsb), [region.rms() for region in regionList])
    
    plt.xlabel("Region Number")
    plt.ylabel("RMS")
    plt.title(path)
        
    plt.savefig(path + outfile + ".png", dpi=120)
    plt.close()
    
def rms_over_time(path, outfile, star, regionMatrix):
    simCount = 0
    simTracker = []
    rmsTracker = []
    
    for regionList in regionMatrix:
        simCount += 1
        for region in regionList:
            simTracker.append(simCount)
            rmsTracker.append(region.rms())
            
    plt.scatter(simTracker, rmsTracker)
    plt.xlabel('Simulation ID')
    plt.ylabel('RMS')
    plt.title('RMS vs Simulation ID')
    
    plt.savefig(path + outfile + ".png", dpi=120)
    plt.close()
    
def twoD_brights_over_time(path, outfile, star, regionList):
    DIMENSION_SCALE = 50
    
    #Set up arrays for brightnesses and positions
    nx = len(regionList[0].get_brights()) * DIMENSION_SCALE #Number of files * 50
    ny = len(regionList[:star.nBoxes]) * DIMENSION_SCALE #Number of boxes and stripes * 50
    img = np.zeros((ny,nx))

    #Populate the pixel values
    idy1 = 0
    idy2 = DIMENSION_SCALE
    idx1 = 0
    idx2 = DIMENSION_SCALE
    
    
    for region in regionList[:star.nBoxes]:
        idx1 = 0
        idx2 = DIMENSION_SCALE
        for brightness in region.get_brights():
            img[idy1:idy2, idx1:idx2] = brightness
            idx1 += DIMENSION_SCALE
            idx2 += DIMENSION_SCALE
        idy1 += DIMENSION_SCALE
        idy2 += DIMENSION_SCALE

    plt.imshow(img, cmap='hot', vmin=0.88, vmax=1.05)
    plt.clim(0.88,1.05)
    plt.colorbar(shrink=0.55)
    plt.xlabel('Window Number')
    plt.xticks(np.arange(0, len(regionList[0].get_brights()), 5) * 50, range(0,len(regionList[0].get_brights()), 5))
    plt.yticks(np.arange(0, star.nBoxes * 50, 100), np.arange(0, 360, (360/star.nBoxes) * 2))
    plt.ylabel('Longitude')
    plt.title('Brightness values over time: ' + path)
    
    plt.savefig(path + outfile + ".png", dpi=180)
    plt.close()
    
    
def make_bright_image(star, regionList, currentWindow):
    
    brights = np.array([region.get_bright_at_time(currentWindow) for region in regionList])
    
    #Set up arrays for brightnesses and positions
    nx = star.nBoxes * 50
    ny = star.nBoxes * 50
    
    bright_stripe = np.array(brights[star.nBoxes:star.nBoxes+star.nStripes])
    bright_box = np.array(brights[0:star.nBoxes])
    
    img = np.zeros((nx,ny))
    
    idx1_stripe = np.arange(0, star.nStripes, 1.)/star.nStripes * nx
    idx2_stripe = (np.arange(0, star.nStripes, 1.) + 1)/star.nStripes * nx - 1
    
    idx1_box = np.arange(0, star.nBoxes, 1.)/star.nBoxes * nx
    idx2_box = (np.arange(0, star.nBoxes, 1.) + 1)/star.nBoxes * nx - 1
    
    idy1_lat = int(star.lat1 * ny/2 + ny/2 + 0.5)
    idy2_lat = int(star.lat2 * ny/2 + ny/2 + 0.5) - 1

    #Populate the pixel values
    for jj in range(star.nStripes):
        idx1 = idx1_stripe[jj]
        idx2 = idx2_stripe[jj]
        idy1 = 0
        idy2 = idy1_lat-1
        img[idy1:idy2, idx1:idx2] = bright_stripe[jj]
        idy1 = idy2_lat+1
        idy2 = ny-1
        img[idy1:idy2, idx1:idx2] = bright_stripe[jj]
    for jj in range(star.nBoxes):
        idx1 = idx1_box[jj]
        idx2 = idx2_box[jj]
        idy1 = idy1_lat
        idy2 = idy2_lat
        img[idy1:idy2, idx1:idx2] = bright_box[jj]
    return img
    
def plot_lc_and_brights(path, outfile, modelLCList, binnedLC, star, regionList, transFileBaseName):
    '''
    Components of this routine:
        Projected brightness map

        Overall binned data in one color
        Windowed binned data in another color
        Model lightcurve in another color

        Transit inset:
            Model portion
            Binned portion        
    '''
    currentWindow = 0
    binMax = max(binnedLC.flux)
    binMin = min(binnedLC.flux)
    
    binDiff = binMax - binMin

    for modelLC in modelLCList:
        #Determine the current window start and stop indices
        start = 0
        while binnedLC.time[start] < modelLC.time[0]:
            start += 1
        stop = start
        while binnedLC.time[stop] < modelLC.time[-1]:
            stop += 1    
        
        #TODO: Change this to actually be nice and pretty rather than dirty and hackish
        #   However, in the meantime... read in the transit files the old way(ew)
        trans_info = path + transFileBaseName + "_%d.out" % currentWindow
        try:
            information = csv.reader(open(trans_info))
        except:
            information = []
        transitCount = 0
        
        for row in information:
            ###############################
            # Make the overall lightcurve #
            ###############################
            plt.plot(modelLC.time, modelLC.flux, c='black', label="Model")
            
            #Order matters here. MPL arranges things such that the most recent is on the top
            plt.scatter(binnedLC.time, binnedLC.flux, c='green', marker = '+', label="Data")
            plt.scatter(binnedLC.time[start:stop], binnedLC.flux[start:stop], c='red', marker='+', label="Brightness Map Region")
            
            plt.ylim(ymax=binMax + (binDiff * .85))
            plt.ylim(ymin=binMin - (binDiff * .1))
            
            plt.xlabel("Orbital Phase")
            plt.ylabel("Relative Flux")
            plt.title("Model Light Curve Vs. Data - Window %d" % currentWindow)
            
            
            ###########################
            # Make the brightness map #
            ###########################
            img = make_bright_image(star, regionList, currentWindow)
            
            ax = plt.axes([.15, .6, .32, .32])
            plt.imsave(path + "/temp.png", img, cmap='hot', vmin=0.85, vmax=1.15)
            plt.imshow(img, cmap='hot')    
            plt.clim(0.85,1.15)
            plt.colorbar(shrink=0.5)    
    
            #Create the plot
            bmap = Basemap(projection='moll', lon_0 = 0)
            bmap.warpimage(path + "/temp.png")
        
            ##################
            # Transit Insets #
            ##################
            transStart = int(row[0]) - start
            transStop = int(row[1]) - start
            
            #TODO: Fix this in the C code.
            #There is a leftover transit relic that is incorrect
            #when starting in transit or ending in transit
            if transStop < transStart:
                plt.close()
                continue  
                          
            transMax = max(modelLC.flux[transStart:transStop])
            transMin = min(modelLC.flux[transStart:transStop])
            
            transDiff = transMax - transMin
            
            transAx = plt.axes([.58,.6,.27,.27])
            plt.plot(modelLC.time[transStart:transStop], modelLC.flux[transStart:transStop], c='black', label="Model")
            plt.scatter(binnedLC.time[transStart + start:transStop + start], binnedLC.flux[transStart + start:transStop + start], c='red', marker='+', label="Data")
            plt.xlabel("Orbital Phase")
            plt.ylabel("Relative Flux")
            plt.title("Transit %d" % transitCount)
            
            plt.ylim(ymax=transMax + (transDiff * .08))
            plt.ylim(ymin=transMin - (transDiff * .08))
            
            #################
            # Save the Plot #
            #################
            fig = plt.gcf()
            fig.set_size_inches(18.5,11.5)
            plt.savefig(path + outfile + "_w%02d_t%02d" % (currentWindow, transitCount) + ".png", dpi=120)
            print "Window:  %d Trans: %d" % (currentWindow, transitCount)
            transitCount += 1
            plt.close()
        
        #Increment the window count for the brightness map maker
        currentWindow += 1
    
