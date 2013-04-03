import os
import errno
import shutil
import numpy as np
from plots_scratch import *

pathList = ["1b/", "1b_1s/", "2b/", "2b_1s/", "2b_2s/", "3b_1s/", "3b_2s/", "3b_3s/", "4b_3s/", "all_b_2s/", "all_varied/"]
star = Star(.0179, .12953, 11, 22)
outfile = open("description.txt", "w")
[outfile.write("%2f" % 1.00) for a in range(star.nsb)]
outfile.close()

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            
def make_test_file(path, star):
    inputfile = open(path + "test.in", "w")
    
    inputfile.write(os.getcwd() + path + "input_lc.in\n")
    inputfile.write("%f\n" % 8)
    inputfile.write("%d\n" % 3)
    inputfile.write("%d\n" % 95)
    inputfile.write("%f\n" % 0)
    inputfile.write("%f\n" % 0.06519)
    inputfile.write("%f\n" % 1.0)
    inputfile.write("%f\n" % 0.0179)
    inputfile.write("%f\n" % 0.12953)
    inputfile.write("%f\n" % 5.67)
    inputfile.write("%f\n" % 0)
    inputfile.write("%f\n" % 0)
    inputfile.write("%d\n" % star.nStripes)
    inputfile.write("%d\n" % star.nBoxes)
    inputfile.write("%d\n" % 1)
    inputfile.write("%f\n" % 7.5)
    inputfile.write("%f\n" % 0.5)
    inputfile.write("%f\n" % 0.420868)
    inputfile.write("%f\n" % 0.244728)

    inputfile.close()
            
if __name__ == '__main__':
    for path in pathList:
        make_sure_path_exists(path)
        
        make_sure_path_exists(path + "no_noise/")
        make_test_file(path + "no_noise/", star)
        
        make_sure_path_exists(path + "12_noise/")
        make_test_file(path + "12_noise/", star)
        
        make_sure_path_exists(path + "14_noise/")
        make_test_file(path + "14_noise/", star)
        
        shutil.copyfile("description.txt", path + "description.txt")

