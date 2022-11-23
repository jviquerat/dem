# Generic imports
import math
import numpy as np

### ************************************************
### Class defining array of particles
class particles:
    ### ************************************************
    ### Constructor
    def __init__(self,
                 n           = 10,
                 nt          = 1000,
                 density     = 1.0,
                 radius      = 1.0,
                 restitution = 0.8,
                 young       = 100.0,
                 poisson     = 0.1,
                 color       = "b",
                 store       = False):

        self.n           = n
        self.nt          = nt
        self.density     = density
        self.radius      = radius
        self.restitution = restitution
        self.young       = young
        self.poisson     = poisson
        self.sigma       = (1.0-poisson**2)/young
        self.vol         = (4.0/3.0)*math.pi*radius**3
        self.mass        = self.vol*self.density
        self.store       = store
        self.color       = color

        self.reset()

    ### ************************************************
    ### Reset arrays
    def reset(self):

        self.m = np.ones((self.n), np.float32)*self.mass        # masses
        self.r = np.ones((self.n), np.float32)*self.radius      # radii
        self.x = np.zeros((self.n,2), np.float32)               # positions
        self.d = np.zeros((self.n,2), np.float32)               # displacements
        self.v = np.zeros((self.n,2), np.float32)               # velocities
        self.a = np.zeros((self.n,2), np.float32)               # accelerations
        self.f = np.zeros((self.n,2), np.float32)               # forces
        self.e = np.ones((self.n), np.float32)*self.restitution # restitution coeff
        self.y = np.ones((self.n), np.float32)*self.young       # young modulus
        self.p = np.ones((self.n), np.float32)*self.poisson     # poisson ratio
        self.c = [self.color]*self.n                # colors

        # Optional storage
        if self.store:
            self.history = np.zeros((self.nt, self.n, 2), np.float32)

        self.set_particles()

    ### ************************************************
    ### Set particle coefficients from inputs
    def set_particles(self):

        # Compute Eb coefficient
        # Eb = (1.0 - p**2)/y
        self.Eb = np.ones((self.n))*(1.0-np.square(self.p[:]))/self.y[:]

        # Compute Gb coefficient
        # Gb = 2*(2 + p)*(1 - p)/y
        self.Gb = np.ones((self.n))*(2.0*(2.0+self.p[:])*(1.0-self.p[:])/self.y[:])

        # Pre-compute g coefficients
        # g =-sqrt(5/6)*2*ln(e)/sqrt(pi**2 + ln(e)**2)
        self.g  = np.ones((self.n))*(-2.0*np.log(self.e[:]))
        self.g *= np.sqrt(1.0/(math.pi**2 + np.square(np.log(self.e[:]))))
        self.g *= math.sqrt(5.0/6.0)

    ### ************************************************
    ### Reset forces
    def reset_forces(self):

        self.f[:,:] = 0.0

    ### ************************************************
    ### Compute maximal velocity
    def max_velocity(self):

        return np.max(np.linalg.norm(self.v))

    ### ************************************************
    ### Compute minimal mass
    def min_mass(self):

        return np.min(self.m)

    ### ************************************************
    ### Compute maximal stiffness
    def max_stiffness(self):

        return np.max((4.0/3.0)*self.y*np.sqrt(self.r))

    ### ************************************************
    ### Compute minimal radius
    def min_radius(self):

        return np.min(self.r)

    ### ************************************************
    ### Compute collisions between particles
    def collisions(self):

        pass

    ### ************************************************
    ### Add gravity
    def gravity(self, g):

        self.a[:,1] -= self.m[:]*g

    ### ************************************************
    ### Update positions using verlet method
    def update(self, dt, it):

        self.f[:,0] /= self.m[:]
        self.f[:,1] /= self.m[:]
        self.v[:,:] += 0.5*dt*(self.a[:,:] + self.f[:,:])
        self.d[:,:]  = dt*self.v[:,:] + 0.5*dt*dt*self.f[:,:]
        self.x[:,:] += self.d[:,:]
        self.a[:,:]  = self.f[:,:]

        if (self.store):
            self.history[it,:,:] = self.x[:,:]
