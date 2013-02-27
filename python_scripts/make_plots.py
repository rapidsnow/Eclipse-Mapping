from plots_scratch import *

###########################
## Edit These Every Time ##
###########################

star = Star(.268, .13031, 6, 18)
outfile = "regions_over_time"
nFiles = 47

pathList = ["./"]

regionMatrix = []

for path in pathList:
    regionList = make_brightness_structures(path + "../description.txt", [path + "brights_%d.out" % fileNumber for fileNumber in range(nFiles)])
    #twoD_brights_over_time(path, 'region_plot', star, regionList)
    #regions_over_time(path, outfile, star, regionList)
    #rms_vs_region(path, "rms_over_time", star, regionList)
    
    binnedLC = Light_Curve(path + "binned.out")
    modelLCList = [Light_Curve(path + "model_%d.out" % number) for number in range(nFiles)]
    plot_lc_and_brights(path, "lc_and_brights", modelLCList, binnedLC, star, regionList, "trans")
