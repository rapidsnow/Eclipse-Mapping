from plots_scratch import *

###########################
## Edit These Every Time ##
###########################

star = Star(.0179, .12953, 11, 22)
outfile = "regions_over_time"
nFiles = 24

fullPathList = ["1b/noise/", "1b_1s/noise/", "2b/noise/", "2b_1s/noise/", "2b_2s/noise/", "3b_1s/noise/", "3b_2s/noise/", "3b_3s/noise/", "4b_3s/noise/", "all_b_2s/noise/", "all_varied/noise/", "1b/14_noise/", "1b_1s/14_noise/", "2b/14_noise/", "2b_1s/14_noise/", "2b_2s/14_noise/", "3b_1s/14_noise/", "3b_2s/14_noise/", "3b_3s/14_noise/", "4b_3s/14_noise/", "all_b_2s/14_noise/", "all_varied/14_noise/", "1b/no_noise/", "1b_1s/no_noise/", "2b/no_noise/", "2b_1s/no_noise/", "2b_2s/no_noise/", "3b_1s/no_noise/", "3b_2s/no_noise/", "3b_3s/no_noise/", "4b_3s/no_noise/", "all_b_2s/no_noise/", "all_varied/no_noise/"]

pathList = ["1b/no_noise/", "all_1s/no_noise/", "darkspot/no_noise/", "mediumspot/no_noise/"]

#regionMatrix = []

for path in pathList:
    regionList = make_brightness_structures(path + "../description.txt", [path + "brights_%d.out" % fileNumber for fileNumber in range(nFiles)])
    #twoD_brights_over_time(path, 'box_plot', star, regionList)
    #twoD_brights_over_time(path, 'stripe_plot', star, regionList, stripes=True)
    
    #regionMatrix.append(regionList)
    #regions_over_time(path, outfile, star, regionList)
    
    #rms_vs_region(path, "rms_over_time", star, regionList)
    make_average_file(path, "average", star, regionList)
    
    #modelLCList = [Light_Curve(path + "model_%d.out" % num) for num in range(nFiles)]
    #binnedLC = Light_Curve(path + "binned.out")
    #plot_lc_and_brights(path, "lc_and_brights", modelLCList, binnedLC, star, regionList, "trans")
    #transit_plots(path, "stacked_transits", modelLCList, binnedLC, star, regionList, "trans")
    
    print path + " Done!"
    
#rms_over_time("./", "rms_over_complexity", star, regionMatrix)

'''
For Kepler-17
path = "../"
nFiles = 41

regionList = make_brightness_structures("1b/description.txt", [path + "brights_%d.out" % fileNumber for fileNumber in range(nFiles)])
regions_over_time(path, outfile, star, regionList)
twoD_brights_over_time(path, 'region_plot', star, regionList)
'''
