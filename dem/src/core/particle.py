# Generic imports
import math
import numpy as np

### ************************************************
### Class defining array of particles
class particles:
    ### ************************************************
    ### Constructor
    def __init__(self,
                 n       = 10,
                 density = 1.0,
                 radius  = 1.0):

        self.n       = n
        self.density = density
        self.radius  = radius
        self.vol     = (4.0/3.0)*math.pi*radius**3
        self.mass    = self.vol*self.density

        self.reset()

    ### ************************************************
    ### Reset arrays
    def reset(self):

        self.m = np.ones((self.n))*self.mass   # masses
        self.r = np.ones((self.n))*self.radius # radii
        self.x = np.zeros((self.n,2))          # positions
        self.v = np.zeros((self.n,2))          # velocities
        self.a = np.zeros((self.n,2))          # accelerations

    ### ************************************************
    ### Compute collisions between particles
    ### Returns a list of colliding particles
    def collisions(self):

        pass

    ### ************************************************
    ### Compute forces
    def forces(self):

        pass

    ### ************************************************
    ### Update positions
    def update(self):

        pass
