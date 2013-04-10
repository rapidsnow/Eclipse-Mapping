import matplotlib.pyplot as plt
import numpy as np
import math
import csv
from mpl_toolkits.basemap import Basemap
 
class Star(object):
    '''
    Container for useful information about the star that can be passed through methods easier than the parameters themselves. 
    Also, calculates lat1 and lat2 for us.
    '''
    def __init__(self, b, Rp, nStripes, nBoxes):
        self.b = b
        self.Rp = Rp
        self.nStripes = nStripes
        self.nBoxes = nBoxes
        self.nsb = nBoxes + nStripes
        self.boxPerStripe = nBoxes/nStripes
        self.lat1 = math.acos(b + Rp)
        self.lat2 = math.acos(b - Rp)

class Region(object):
    '''
    Container that represents a region over time. 
    Each value in the brights array represents the brightness of this region for a different (usually sequential) window.
    It is best to create this via the make_brightness_structures method rather than doing it by hand.
    goalVal only applies to synthetic lightcurves. It corresponds to the region's value that was used to create a given synthetic curve
    
    Has a few useful methods to operate on the data. The most useful one is rms which returns the RMS of the brightness of the region over all windows versus the stated goalVal
    '''
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
    Takes a filename and creates a lightcurve object with flux, time, and error arrays
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
    
class Transit_Container(object):
    def __init__(self, start, stop):
        self.timeList = []
        self.fluxList = []
        self.errorList = []
        self.startOffsetList = []
        self.start = start
        self.stop = stop
    def add_transit(self, lightCurve, indexDiff):
        transStart = max([self.start - indexDiff, 0])
        transStop = min([self.stop - indexDiff, len(lightCurve.get_flux()) - 1])
    
        self.timeList.append(lightCurve.get_time()[transStart:transStop])
        self.fluxList.append(lightCurve.get_flux()[transStart:transStop])
        self.errorList.append(lightCurve.get_error()[transStart:transStop])
        self.startOffsetList.append(transStart - self.start + indexDiff) #I think this should give the initial offset
    def get_mean_time(self):
        '''
        Returns the longest list of times
        (so that you don't get a partial transit as the representative timeList)
        '''
        lenList = [len(time) for time in self.timeList]
        index = lenList.index(max(lenList))
        return self.timeList[index]
    def get_mean_flux(self):
        return self.get_mean(self.fluxList)
    def get_mean_error(self):
        return self.get_mean(self.errorList)
    def get_mean(self, twoDArray):
        '''
        Returns a mean array that averages every value in the 2d array
        The 2d array does not need to have the same number of elements in every list
        '''
        meanList = []
        for j in range(len(twoDArray[0])):
            total = 0
            count = 0
            for i in range(len(twoDArray)):
                try:
                    total += twoDArray[i][j + self.startOffsetList[i]]
                    count += 1
                except(IndexError):
                    pass
            meanList.append(total/float(count))
        return meanList
#################################################################################################################
    
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

#################################################################################################################

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
        
#################################################################################################################

def rms_vs_region(ax, path, outfile, star, regionList):
    regionNum = 0
    
    regionTracker = []
    regionAbbTracker = []
    rmsTracker = []
    abberantTracker = []
    
    for region in regionList:
        regionNum += 1
        if abs(1.0 - region.goalVal) < .000000001:
            rmsTracker.append(region.rms())
            regionTracker.append(regionNum)
        else:
            abberantTracker.append(region.rms())
            regionAbbTracker.append(regionNum)
                
    ax.scatter(regionTracker, rmsTracker, marker='o', color='b')    
    ax.scatter(regionAbbTracker, abberantTracker, marker='D', color='R')
    ax.plot([star.nBoxes + .5, star.nBoxes + .5], [0, 1.1 * max(rmsTracker + abberantTracker)], color='black')
    
    ax.set_xlim(xmin = -1)
    ax.set_xlim(xmax = star.nsb + 1)
    ax.set_ylim(ymin = 0)
    ax.set_ylim(ymax = 1.1 * max(rmsTracker + abberantTracker))
    
    ax.set_xlabel("Region")
    ax.set_ylabel("RMS")
    ax.set_title("RMS")
        
#    plt.savefig(path + outfile + ".png", dpi=120)
#    plt.close()

#################################################################################################################

def rms_over_time(ax, path, outfile, star, regionMatrix):
    simCount = 0
    simTracker = []
    rmsTracker = []
    
    for regionList in regionMatrix:
        simCount += 1
        for region in regionList:
            rmsTracker.append(region.rms())
            simTracker.append(simCount)
            
    ax.scatter(simTracker, rmsTracker, marker='o', color='b')
    ax.set_xlabel('Simulation ID')
    ax.set_ylabel('RMS')
    #ax.set_title('RMS vs Simulation ID')
    
#    plt.savefig(path + outfile + ".png", dpi=120)
#    plt.close()

#################################################################################################################

def twoD_brights_over_time(ax, path, outfile, star, regionList, stripes=False):
    '''
    Creates the box and stripe plots
    path: path to output
    outfile: base name of file, do not include the path or the extension name
    star: a Star object from above
    regionList: a list of Region objects from above
    stripes: Do you want a stripe or a box plot? The default is box plot
    '''
    DIMENSION_SCALE = 50
    
    if stripes:
        startIndex = star.nBoxes
        stopIndex = star.nsb
        cm = 'hot'
        #cm = 'RdGy'
    else:
        startIndex = 0
        stopIndex = star.nBoxes
        cm = 'hot'
        #cm = 'RdGy'
    nRegions = stopIndex - startIndex
    
    #Set up arrays for brightnesses and positions
    nx = (len(regionList[0].get_brights()) + 1) * DIMENSION_SCALE + 4 #Number of files * 50 + 1 goal column + 4 pixels to separate
    ny = len(regionList[startIndex:stopIndex]) * DIMENSION_SCALE #Number of boxes and stripes * 50
    img = np.zeros((ny,nx))

    #Populate the pixel values
    idy1 = 0
    idy2 = DIMENSION_SCALE
    idx1 = 0
    idx2 = DIMENSION_SCALE
    
    
    for region in regionList[startIndex:stopIndex]:
        idx1 = 0
        idx2 = DIMENSION_SCALE
        for brightness in region.get_brights():
            img[idy1:idy2, idx1:idx2] = brightness
            idx1 += DIMENSION_SCALE
            idx2 += DIMENSION_SCALE
        img[idy1:idy2, nx - DIMENSION_SCALE:nx] = region.goalVal
        idy1 += DIMENSION_SCALE
        idy2 += DIMENSION_SCALE


    tempImg = ax.imshow(img, cmap=cm, vmin=0.88, vmax=1.05)
    if stripes:
        plt.colorbar(tempImg, shrink=0.75, ax=ax)

    ax.set_xlabel('Window Number')
    ax.set_xticks(np.arange(0, len(regionList[0].get_brights()), 5) * 50)
    ax.set_xticklabels(range(0,len(regionList[0].get_brights()), 5))
    
    ax.set_yticks(np.arange(0, nRegions * 50, 100))
    ax.set_yticklabels(np.arange(startIndex, stopIndex, 2))
    ax.set_ylabel('Region')
    
    #ax.set_title('Brightness values over time: ' + path)
    ax.set_title('Brightness values')
    
#    plt.savefig(path + outfile + ".png", dpi=180)
#    plt.close()
    
#################################################################################################################

def make_bright_image(star, regionList, currentWindow, goal=False):
    '''
    Utility that makes the matrix for the brightness map. Shouldn't really be called on its own
    
    For real data: change the brights array to the commented list comprehension
    rather than the if/else statements
    '''
    #brights = np.array([region.get_bright_at_time(currentWindow) for region in regionList])
    
    if goal:
        brights = np.array([region.goalVal for region in regionList])
    else:
        brights = np.array([region.get_mean() for region in regionList])
        
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
    
    idy1_lat = int((star.lat1/math.pi) * ny/2 + ny/2 + 0.5)
    idy2_lat = int((star.lat2/math.pi) * ny/2 + ny/2 + 0.5) - 1

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
    
#################################################################################################################

def plot_lc(ax, path, outfile, modelLCList, binnedLC, star, regionList, transFileBaseName):
    '''
    Components of this routine:
        Overall binned data in one color
        Windowed binned data in another color
        Model lightcurve in another color
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
        
        ###############################
        # Make the overall lightcurve #
        ###############################
        
        #Order matters here. MPL arranges things such that the most recent is on the top
        #ax.scatter(binnedLC.time, binnedLC.flux, c='green', marker = '+', label="Data")
        ax.scatter(binnedLC.time[start:stop], binnedLC.flux[start:stop], c='red', marker='o', s=5, lw=0, label="Brightness Map Region")
        ax.plot(modelLC.time, modelLC.flux, c='black', label="Model")
        
        ax.set_ylim(ymax=binMax + .05 * binDiff)
        ax.set_ylim(ymin=binMin - .02 * binDiff)
        ax.set_xlim(xmax=binnedLC.time[start])
        ax.set_xlim(xmax=binnedLC.time[stop] + (binnedLC.time[1] - binnedLC.time[0]) * 10)
        
        ax.set_xlabel("Orbital Phase")
        ax.set_ylabel("Relative Flux")
        ax.set_title("Light Curve Fit")
#        plt.savefig(path + "model_fit_w%d.png" % currentWindow) 
#        plt.close()
            
        #Increment the window count for the brightness map maker
        currentWindow += 1
        
################################################################################################################# 

def plot_brights(ax, path, star, regionList, goal=False):
    '''
    Components of this routine:
        Projected brightness map
        
    Please note that this has been modified for use in diagnostic plots, 
    there should really be a way to specify a windowNumber for real data
    '''
    currentWindow = 0

    ###########################
    # Make the brightness map #
    ###########################
    img = make_bright_image(star, regionList, currentWindow, goal=goal)
    
    plt.imsave(path + "temp.jpg", img, cmap='hot', vmin=0.85, vmax=1.15)
    plt.imshow(img, cmap='hot')
    #Create the plot
    bmap = Basemap(projection='moll', lon_0 = 0, ax=ax)
    bmap.warpimage(path + "temp.jpg", ax=ax)
    
    if goal:
        ax.set_title("Desired Map")
    else:
        ax.set_title("Average Map")

#################################################################################################################
def transit_plots(axList, path, modelLCList, binnedLC, transFileBaseName):
    '''
    Return a list of axes objects containing transit plots for each window
    TODO: Add a way to keep track of axes.
        Either pass in a list of axes or create one in here.
        Maybe an automatic gridspec type of thing?
    '''
    transInfo = path + transFileBaseName + ".out"
    try:
        information = csv.reader(open(transInfo))
    except:
        information = []
        print "Warning, no transit information provided\nFile %s not opened or contains no information.\n" % transInfo

    currentAx = 0
    transits = []
    for row in information:
        transits.append(Transit_Container(int(row[0]), int(row[1])))
        for modelLC in modelLCList:
            #Determine the current window start and stop indices
            start, stop = get_start_and_stop(binnedLC, modelLC)
            if transit_in_window(start, stop, row):
                #Add this window's transit to the appropriate container
                transits[-1].add_transit(modelLC, start)
            else:
                continue
            '''
            Just a note here on what is done so that I don't get lost.
            I now have made a list of Transit_Containers that is nTransits long (length of the transit file)
            Each element of that list should contain time (redundant), flux, and error values for each window
            You can recover them individually by hacking into the list, or you can average them
            Now we should be able to plot
            '''
        #Now we are out of the modelLC loop, but still in the row loop. Every appropriate transit should be in transits[-1] right now
        #axList[currentAx].plot(transits[-1].get_mean_time(), transits[-1].get_mean_flux(), c='black', lw=1.5, label="Average")
        for i in range(len(transits[-1].fluxList)):
            axList[currentAx].plot(transits[-1].timeList[i], transits[-1].fluxList[i], c='gray', lw=.75)
            
        axList[currentAx].scatter(binnedLC.time[int(row[0]):int(row[1])], binnedLC.flux[int(row[0]):int(row[1])], c='red', marker='+', label='Data')

        currentAx += 1

def transit_in_window(start, stop, (trans_start, trans_stop)):
    '''
    Helper Function for transit plots
    Returns True if the given transit limits are in between start and stop
    Returns False otherwise
    '''
    return (int(trans_start) <= stop) and (int(trans_stop) >= start)

def get_start_and_stop(binnedLC, modelLC):
    '''
    Helper Function for transit plots
    Returns the start and stop indices for the binned LC of a given model LC
    '''
    start = 0
    while binnedLC.time[start] < modelLC.time[0]:
        start += 1
    stop = start
    while binnedLC.time[stop] < modelLC.time[-1]:
        stop += 1  
    return start, stop

#################################################################################################################

def noise_comp(input_vector):
    '''
    Stupid little method that I wrote to plot out the different levels of noise for a paper. 
    It has no scientific value really and is only in here in case we write a paper about a different object.
    '''
    LC_List = []
    scale = 0
    for inp in input_vector:
        LC_List.append(Light_Curve(inp))
        
    colorArray = ['b','r','g']
    iterator = 0
    for LC in LC_List:
        plt.scatter(LC.get_time(), LC.get_flux() + scale, marker='.', color=colorArray[iterator])
        iterator += 1
        scale += .033
        
    plt.title("Noise levels for synthetic light curves")
    plt.xlabel("Phase")
    plt.ylabel("Relative Flux\n(Offset for each curve to visualize comparative noise)")
    
    plt.xlim(xmin=LC_List[0].get_time()[0], xmax=LC_List[0].get_time()[-1])
    
    plt.savefig("noise_levels.png", dpi=180)
    
#################################################################################################################

def make_average_file(path, output, star, regionList):
    '''
    Some quantifiable measure of the average of a region over time, not actually a plot
    '''
    boxAve = 0
    stripeAve = 0
    darkBoxTotal = 0
    darkStripeTotal = 0
    darkStripe = 0
    darkBox = 0
    
    e_region = ((math.sin(math.pi/2) - math.sin(-math.pi/2)) * (star.lat2 - star.lat1 - math.sin(star.lat2)*math.cos(star.lat2) + math.sin(star.lat1)*math.cos(star.lat1)))/(2 * math.pi);
    
    outfile = open(path + output + ".txt", "w")
    
    for region in regionList[:star.nBoxes]:
        if abs(1.0 - region.goalVal) < .00000001:
            boxAve += region.get_mean()
            outfile.write("%f\n" % region.get_mean())
        else:
            darkBox += 1
            darkBoxTotal += region.get_mean()
            outfile.write("%f\n" % region.get_mean())
            
    for region in regionList[star.nBoxes:]:
        if abs(1.0 - region.goalVal) < .00000001:
            stripeAve += region.get_mean()
            outfile.write("%f\n" % region.get_mean())
        else:
            darkStripe += 1
            darkStripeTotal += region.get_mean()
            outfile.write("%f\n" % region.get_mean())
    
    nBoxAve = boxAve/float(star.nBoxes - darkBox)
    boxAve = (boxAve + darkBoxTotal)/float(star.nBoxes)
    
    nStripeAve = stripeAve/float(star.nStripes - darkStripe)
    stripeAve = (stripeAve + darkStripeTotal)/float(star.nStripes)
    
    nTotalAve = nBoxAve * e_region + nStripeAve * (1 - e_region)
    totalAve = boxAve * e_region + stripeAve * (1 - e_region)
    
    outfile.write("\n------------ Totals ------------\nTotal Average: %f\nBox Average: %f\nStripe Average: %f\n\nNon-Darkened Total Average: %f\nNon-Darkened Box Average: %f\nNon-Darkened Stripe Average: %f\n" % (totalAve, boxAve, stripeAve, nTotalAve, nBoxAve, nStripeAve))
        
    outfile.close()
