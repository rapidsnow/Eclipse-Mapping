from plots_scratch import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
    fig = plt.figure(figsize=(10,12), dpi=500)
#    axBoxPlot = fig.add_subplot(2,2,1)
#    axStripePlot = fig.add_subplot(2,2,2)
#    axLC = fig.add_subplot(2,2,3)
#    axRMS = fig.add_subplot(2,2,4)
#    axModelMap = fig.add_subplot(3,2,1) #nrows, ncolumns, plot number (index at 1)
#    axAverageMap = fig.add_subplot(3,2,2)

    gs = gridspec.GridSpec(4,2, hspace=.35)#, height_ratios=[1,5,5,5,5])

    #axTitle = fig.add_subplot(gs[0,:])
    axModelMap = fig.add_subplot(gs[0,0])
    axAverageMap = fig.add_subplot(gs[0,1])

    regionList = make_brightness_structures(path + "description.txt", [path + "brights_%d.out" % fileNumber for fileNumber in range(nFiles)])
    
    modelLCList = [Light_Curve(path + "model_%d.out" % num) for num in range(nFiles)]
    binnedLC = Light_Curve(path + "binned.out")
    
    #Just a note here. The Basemap package is really poorly designed to interface with the rest of MPL
    #You need to do this part first and then define the rest of the axes in order for it to not mess anything up.
    plot_brights(axModelMap, path, star, regionList, goal=True)
    plot_brights(axAverageMap, path, star, regionList, goal=False)
    
    axBoxPlot = fig.add_subplot(gs[1, 0])
    axStripePlot = fig.add_subplot(gs[1,1])
    axRMS = fig.add_subplot(gs[2,:])
    axLC = fig.add_subplot(gs[3,:])
    
    plot_lc(axLC, path, "lc_and_brights", modelLCList, binnedLC, star, regionList, "trans")
    
    twoD_brights_over_time(axBoxPlot, path, 'box_plot', star, regionList)
    twoD_brights_over_time(axStripePlot, path, 'stripe_plot', star, regionList, stripes=True)
    
    rms_vs_region(axRMS, path, "rms_over_time", star, regionList)
    
#    axTitle.set_title("Diagnostic Plots: %s" % path)
#    axTitle.get_xaxis().set_visible(False)
#    axTitle.get_yaxis().set_visible(False)
#    for sp in axTitle.spines.values():
#        sp.set_visible(False)
    plt.savefig("page_plot.eps", format='eps')
    plt.close()

    axList = [0]*20
    fig2 = plt.figure(figsize=(8.5, 11), dpi=500,)
    gs2 = gridspec.GridSpec(5,4,wspace=0,hspace=0)
    for i in range(5):
        for j in range(4):
            axList[j + i * 4] = fig2.add_subplot(gs2[i, j])

    transit_plots(axList, path, modelLCList, binnedLC, "trans")
    for ax in axList:
        ax.set_xticks([])
        ax.set_yticks([])
    plt.savefig("transit_page.eps", format='eps')
    
    make_average_file(path, "average", star, regionList)
    
    print path + " Done!"
