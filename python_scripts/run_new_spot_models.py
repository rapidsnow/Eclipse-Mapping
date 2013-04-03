import subprocess

#pathList = ["1b/noise/", "1b_1s/noise/", "2b/noise/", "2b_1s/noise/", "2b_2s/noise/", "3b_1s/noise/", "3b_2s/noise/", "3b_3s/noise/", "4b_3s/noise/", "all_b_2s/noise/", "all_varied/noise/", "1b/14_noise/", "1b_1s/14_noise/", "2b/14_noise/", "2b_1s/14_noise/", "2b_2s/14_noise/", "3b_1s/14_noise/", "3b_2s/14_noise/", "3b_3s/14_noise/", "4b_3s/14_noise/", "all_b_2s/14_noise/", "all_varied/14_noise/", "1b/no_noise/", "1b_1s/no_noise/", "2b/no_noise/", "2b_1s/no_noise/", "2b_2s/no_noise/", "3b_1s/no_noise/", "3b_2s/no_noise/", "3b_3s/no_noise/", "4b_3s/no_noise/", "all_b_2s/no_noise/", "all_varied/no_noise/"]

pathList = ["1b/no_noise/", "darkspot/no_noise/", "mediumspot/no_noise/"]
for path in pathList:
    subprocess.call(["./main", path + "test.in", path])
    print path + " Done!"
    
subprocess.call(["python", "make_plots.py"])
