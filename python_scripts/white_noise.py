'''
Python white noise generator
Inputs: File with three columns: Time Flux Error
Outputs: File with the same three columns with white noise added
'''
import random
import sys


def mag_to_index(mag):
    '''
    Converts a given magnitude value to a corresponding 
    index in the sigmas array
    '''
    return int((float(mag) - 7) * 2)

def read_three_or_two(filename):
    '''
    This will try to read a three-column file and return phase, flux, and error lists.
    If it is a two column file, the error list will be empty.
    '''
    error = []
    flux = []
    phase = []

    #Open the file and read in the arrays
    for line in open(filename, 'r'):
        temp = line.split(' ')
    
        phase.append(float(temp[0]))
        flux.append(float(temp[1]))
        #If it is only two columns, that's okay
        #This is assuming that the user puts in a magnitude value for two columned data.
        #Either way, this won't raise an error now
        try:
            error.append(float(temp[2]))
        except IndexError:
            error.append(0)
    return phase, flux, error

def generate_noise(flux, error, mag=None):
    '''
    Uses the error counts provided by Kepler to simulate Gaussian noise on a noiseless synthetic lightcurve.
    '''
    #Count errors for Kepler light curves.
    #Magnitudes start at 7.0 for index 0 and go up by .5 for each successive index
    sigmas = [float(x) for x in [39, 49, 62, 78, 99, 125, 159, 202, \
                259, 336, 441, 589, 804, 1127, 1620, 2388, 3594, 5495,\
                8503, 13270, 20810, 32770, 51720, 81760, 129400, 204800,\
                324400, 513900, 814300]]
                
    #If a magnitude is given, the argv should have length > 2
    if mag == None:
        #Create noise based on the errors given at each datapoint
        newFlux = [flux[x] + random.gauss(0, error[x]) for x in range(len(flux))]
    else:
        #Find which sigma we are using
        sigma = sigmas[mag_to_index(mag)]/(1e6)
        #Generate the random noise
        newFlux = [x + random.gauss(0, sigma) for x in flux]
        #Create the error array
        error = [sigma]*len(newFlux)
    return newFlux, error

def write_noise(output, phase, flux, error):
    #Write the noise file
    outfile = open(output, 'w')
    for i in range(len(flux)):
        outfile.write(str(phase[i]) + ' ' + str(flux[i]) + ' ' + str(error[i]) + '\n')
    
if __name__ == '__main__':
    filename = sys.argv[1]
    try:
        mag = sys.argv[2]
    except IndexError:
        mag=None
    
    phase, flux, error = read_three_or_two(filename)
    flux, error = generate_noise(flux, error, mag)
    write_noise('noise.out', phase, flux, error)
    
