import sys
import shutil
import subprocess
import data_plots

inputFile = sys.argv[1]
outfile = open('test.in', 'r')
output = []

for line in open(inputFile):
    if line[0] != '#':
        output.append(line)
        
for line in output:
    outfile.write(line)
    
subprocess.call('main', 'test.in', '.')
data_plots.plot_lc_and_brights('STUFF HERE!!!!!')

subprocess.call('leslieCode', 'test.in', '.')
data_plots.plot_lc_and_brights('MORE STUFF HERE!!!!')

shutil.copy('ALLTHEFILES!")
