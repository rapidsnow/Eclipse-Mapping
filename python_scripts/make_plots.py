#import data_plots as dp
from plots_scratch import *

#for ii in range(39):
#	dp.plot_lc_and_brights("1b/noise/", "brights_%d.out" % ii, "chi_%d.out" % ii, "model_%d.out" % ii, "spot_model.out", .268, .13031, 7, 21, "brights_lc_test_%02d" % ii, ii, "trans_%d.out" % ii)

###########################
## Edit These Every Time ##
###########################

star = Star(.268, .13031, 6, 18)
outfile = "regions_over_time"
nFiles = 27

pathList = ["1b/noise/", "1b_1s/noise/", "2b/noise/", "2b_1s/noise/", "2b_2s/noise/", "3b_1s/noise/", "3b_2s/noise/", "3b_3s/noise/", "4b_3s/noise/", "all_b_2s/noise/", "all_varied/noise/", "1b/14_noise/", "1b_1s/14_noise/", "2b/14_noise/", "2b_1s/14_noise/", "2b_2s/14_noise/", "3b_1s/14_noise/", "3b_2s/14_noise/", "3b_3s/14_noise/", "4b_3s/14_noise/", "all_b_2s/14_noise/", "all_varied/14_noise/", "1b/no_noise/", "1b_1s/no_noise/", "2b/no_noise/", "2b_1s/no_noise/", "2b_2s/no_noise/", "3b_1s/no_noise/", "3b_2s/no_noise/", "3b_3s/no_noise/", "4b_3s/no_noise/", "all_b_2s/no_noise/", "all_varied/no_noise/"]

#regionMatrix = []

for path in pathList:
    regionList = make_brightness_structures(path + "../description.txt", [path + "brights_%d.out" % fileNumber for fileNumber in range(nFiles)])
    #twoD_brights_over_time(path, 'box_plot', star, regionList)
    #twoD_brights_over_time(path, 'stripe_plot', star, regionList, stripes=True)
    #regionMatrix.append(regionList)
    #regions_over_time(path, outfile, star, regionList)
    #rms_vs_region(path, "rms_over_time", star, regionList)
    
    modelLCList = [Light_Curve(path + "model_%d.out" % num) for num in range(nFiles)]
    binnedLC = Light_Curve(path + "binned.out")
    plot_lc_and_brights(path, "lc_and_brights", modelLCList, binnedLC, star, regionList, "trans")
    print path + " Done!"
#rms_over_time("./", "rms_over_time", star, regionMatrix)

#path = "../"
#nFiles = 41

#regionList = make_brightness_structures("1b/description.txt", [path + "brights_%d.out" % fileNumber for fileNumber in range(nFiles)])
#regions_over_time(path, outfile, star, regionList)
#twoD_brights_over_time(path, 'region_plot', star, regionList)
