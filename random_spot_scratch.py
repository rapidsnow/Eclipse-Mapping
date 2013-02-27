import numpy as np
import math
import random
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

PIXEL_PRECISION = 1000

'''
Assuming the input is a list of spots, let's make a spot class
'''
class Spot(object):
    '''
    Stores information about a given starspot
    Can output a rectangular spot model representation
    Inputs:
        theta: angle from 0 -> pi
        phi:   angle from 0 -> 2pi
        r:     radius of a spot from 0 -> 1. This is actually r_{spot}/r{star}
        val:   relative brightness value compared to normal
    '''
    def __init__(self, theta, phi, r, val):
        self.theta = theta
        self.phi = phi
        self.r = r
        self.val = val
    def to_rectangle(self, rectangle):
        '''
        Convert from spherical to rectangular projection coordinates and add them to a rectangular field
        '''
        #Center
        cX = int((self.phi/(2 * math.pi)) * rectangle.get_width())
        cY = int((self.theta/(math.pi)) * rectangle.get_height())
        
        #Limits
        x1 = cX - int(self.r * rectangle.get_width())
        x2 = cX + int(self.r * rectangle.get_width())
        y1 = cY - int(self.r * rectangle.get_height())
        y2 = cY + int(self.r * rectangle.get_height())
        
        rectangle.add_spot(x1, x2, y1, y2, cX, cY, self.r * PIXEL_PRECISION, self.val)
        
        
class Rectangle(object):
    '''
    The field on which our spots are set
    '''
    def __init__(self):
        self.width = 1
        self.height = 2
        self.matrix = np.array([random.gauss(1, .04) for i in range(self.width * PIXEL_PRECISION * self.height * PIXEL_PRECISION)]).reshape(self.width * PIXEL_PRECISION, self.height * PIXEL_PRECISION)
        
#        self.matrix = np.ones(self.width * PIXEL_PRECISION * self.height * PIXEL_PRECISION).reshape(self.width * PIXEL_PRECISION, self.height * PIXEL_PRECISION)
    def get_width(self):
        return self.width * PIXEL_PRECISION
    def get_height(self):
        return self.height * PIXEL_PRECISION
    def add_spot(self, x1, x2, y1, y2, x0, y0, r, val):
        '''
        Darkens the appropriate pixels in the rectangular representation
        '''
        for i in range(x1, x2):
            for j in range(y1, y2):
                if within_r((x0, y0), (i, j), r):
                    if i > 0:
                        i = i % self.get_width()
                    if j > 0:
                        j = j % self.get_height()
                    self.matrix[i][j] -= (1 - val) #Did this rather than pure assignment so that spots can overlap
    def proj_plot(self):
        '''
        Creates a Mollweide projection of Rectangle and saves it
        '''
        #ax = plt.axes([.15, .6, .32, .32])
        plt.imsave("/temp.png", self.matrix, cmap='hot', vmin=0.85, vmax=1.15)
        plt.imshow(self.matrix, cmap='hot')    
        plt.clim(0.85,1.15)
        plt.colorbar(shrink=0.5)

        #Create the plot
        bmap = Basemap(projection='robin', lon_0 = 0)
        bmap.warpimage("/temp.png")
        
        fig = plt.gcf()
        fig.set_size_inches(18.5,11.5)
        plt.savefig("spot_proj.png", dpi=120)

def within_r(p1, p2, r):
    '''
    Takes two points p1 and p2 and a radius value r. Returns True if |p2 - p1| <= r
    '''
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2) <= r
    
if __name__ == '__main__':
    import random
    field = Rectangle()
    spots = [Spot(random.uniform(0, math.pi), random.uniform(0, 2 * math.pi), random.uniform(0.03, .15), random.uniform(.80, 1.05)) for i in range (12)]
    for spot in spots:
        spot.to_rectangle(field)
    
    field.proj_plot()
