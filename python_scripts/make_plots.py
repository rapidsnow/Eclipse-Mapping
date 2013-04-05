from plots_scratch import *
import matplotlib.pyplot as plt
import numpy as np

###########################
## Edit These Every Time ##
###########################
star = Star(.0179, .12953, 11, 22)
nFiles = 15
###########################
###########################

fullPathList = ["1b/noise/", "1b_1s/noise/", "2b/noise/", "2b_1s/noise/", "2b_2s/noise/", "3b_1s/noise/", "3b_2s/noise/", "3b_3s/noise/", "4b_3s/noise/", "all_b_2s/noise/", "all_varied/noise/", "1b/14_noise/", "1b_1s/14_noise/", "2b/14_noise/", "2b_1s/14_noise/", "2b_2s/14_noise/", "3b_1s/14_noise/", "3b_2s/14_noise/", "3b_3s/14_noise/", "4b_3s/14_noise/", "all_b_2s/14_noise/", "all_varied/14_noise/", "1b/no_noise/", "1b_1s/no_noise/", "2b/no_noise/", "2b_1s/no_noise/", "2b_2s/no_noise/", "3b_1s/no_noise/", "3b_2s/no_noise/", "3b_3s/no_noise/", "4b_3s/no_noise/", "all_b_2s/no_noise/", "all_varied/no_noise/"]

pathList = ["./"]

for path in pathList:
    fig = plt.figure()
    axModelMap = fig.add_subplot(3,2,1) #nrows, ncolumns, plot number (index at 1)
    axAverageMap = fig.add_subplot(3,2,2)
    axBoxPlot = fig.add_subplot(3,2,3)
    axStripePlot = fig.add_subplot(3,2,4)
    axLC = fig.add_subplot(3,2,5)
    axRMS = fig.add_subplot(3,2,6)

    regionList = make_brightness_structures(path + "description.txt", [path + "brights_%d.out" % fileNumber for fileNumber in range(nFiles)])
    twoD_brights_over_time(axBoxPlot, path, 'box_plot', star, regionList)
    twoD_brights_over_time(axStripePlot, path, 'stripe_plot', star, regionList, stripes=True)
    
    rms_vs_region(axRMS, path, "rms_over_time", star, regionList)
    make_average_file(path, "average", star, regionList)
    
    modelLCList = [Light_Curve(path + "model_%d.out" % num) for num in range(nFiles)]
    binnedLC = Light_Curve(path + "binned.out")
    
    plot_lc(axLC, path, "lc_and_brights", modelLCList, binnedLC, star, regionList, "trans")
    plot_brights(axModelMap, path, "lc_and_brights", modelLCList, binnedLC, star, regionList, "trans")
    plot_brights(axAverageMap, path, "lc_and_brights", modelLCList, binnedLC, star, regionList, "trans")
#    transit_plots(path, "stacked_transits", modelLCList, binnedLC, star, regionList, "trans")
    
    plt.show()
    plt.close()
    
    print path + " Done!"
