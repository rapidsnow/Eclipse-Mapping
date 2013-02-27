'''
This code provides utilities to take a set of visibility 
curves and a list of brightness values and produce a model lightcurve
with and without desired gaussian noise
'''

import white_noise as wn

class Visibility(object):
    def __init__(self, path, nBoxes, nStripes):
        self.path = path
        self.nBoxes = nBoxes
        self.nStripes = nStripes
        self.q = nBoxes/nStripes
        self.visMatrix, self.phase = self.open_files()
        self.make_S_vals()
    def open_files(self):
        visMatrix = []
        for i in range(self.nBoxes):
            phase, vis, dummy = wn.read_three_or_two(self.path + "b%d.out" % i)
            visMatrix.append(vis)
        for i in range(self.nStripes):
            phase, vis, dummy = wn.read_three_or_two(self.path + "z%d.out" % i)
            visMatrix.append(vis)
        return visMatrix, phase
    def make_S_vals(self):
        '''
        This is ported from my C Code, sorry for poor readability
        Theoretically, it should take the visibility matrix, which
        is in B-Z space and turn it into B-S space
        '''
        for i in range(len(self.visMatrix[0])):
            for j in range(self.nBoxes, self.nBoxes + self.nStripes):
                tempVis = 0
                for k in range(self.q):
                    tempVis += self.visMatrix[(j - self.nBoxes) * self.q + k][i]
                self.visMatrix[j][i] -= tempVis
    def get_visMatrix(self):
        return self.visMatrix[:][:]
    def get_phase(self):
        return self.phase[:]

class Brightness(object):
    '''
    Just an automatic container for a list of brightness values
    '''
    def __init__(self, path, filename):
        self.brights = [float(bright.rstrip('\r\n')) for bright in open(path + filename)]
    def get_brights(self):
        return self.brights[:]

def calculate_model(vis, brights):
    visMatrix = vis.get_visMatrix()
    brights = brights.get_brights()
    
    flux = []
    
    for i in range(len(visMatrix[0])):
        tempFlux = 0
        for j in range(len(brights)):
            tempFlux += visMatrix[j][i] * brights[j]
        flux.append(tempFlux)
    return flux
    
def make_noise_curve(path, outfile, mag, phase, flux, error=[]):
    flux, error = wn.generate_noise(flux, error, mag)
    make_clean_curve(path, outfile, phase, flux, error)
    
def make_clean_curve(path, outfile, phase, flux, error=[]):
    writeTo = open(path + outfile, "w")
    for i in range(len(phase)):
        writeTo.write("%lf %lf %lf\n" % (phase[i], flux[i], error[i]))
        
if __name__ == '__main__':
    import sys
    path = sys.argv[1]
    
    error = [float(temp.split()[2]) for temp in open("../binned.out")]
    
    visMatrix = Visibility("../Vis_Plots/", 18, 6)
    brights = Brightness(path, "description.txt")
    
    flux = calculate_model(visMatrix, brights)
    phase = visMatrix.get_phase()
    
    make_noise_curve(path, "12_noise_model.out", 12, phase, flux, error)
    make_clean_curve(path, "model.out", phase, flux, error)
