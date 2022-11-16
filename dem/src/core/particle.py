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
                 density     = 1.0,
                 radius      = 1.0,
                 restitution = 0.8,
                 young       = 100.0,
                 poisson     = 0.1,
                 color       = "b",
                 store       = False):

        self.n           = n
        self.density     = density
        self.radius      = radius
        self.restitution = restitution
        self.young       = young
        self.poisson     = poisson
        self.sigma       = (1.0-poisson**2)/young
        self.vol         = (4.0/3.0)*math.pi*radius**3
        self.mass        = self.vol*self.density
        self.g           = 9.81
        self.store       = store

        self.reset()

    ### ************************************************
    ### Reset arrays
    def reset(self):

        self.m = np.ones((self.n))*self.mass        # masses
        self.r = np.ones((self.n))*self.radius      # radii
        self.x = np.zeros((self.n,2))               # positions
        self.d = np.zeros((self.n,2))               # displacements
        self.v = np.zeros((self.n,2))               # velocities
        self.a = np.zeros((self.n,2))               # accelerations
        self.f = np.zeros((self.n,2))               # forces
        self.e = np.ones((self.n))*self.restitution # restitution coeff
        self.y = np.ones((self.n))*self.young       # young modulus
        self.p = np.ones((self.n))*self.poisson     # poisson ratio

        # Optional storage
        if self.store:
            self.history = np.empty((2*self.n), float)

        self.set_particles()

    ### ************************************************
    ### Set particle coefficients from inputs
    def set_particles(self):

        # Compute sigma coefficient
        self.sigma = np.ones((self.n))*(1.0-np.square(self.p[:]))/self.y[:]

        # Compute kappa coefficient
        self.kappa = np.ones((self.n))*(2.0*(2.0+self.p[:])*(1.0-self.p[:])/self.y[:])

        # Compute alpha coefficient
        self.alpha  = np.ones((self.n))*(-2.0*np.log(self.e[:]))
        self.alpha *= np.sqrt(1.0/(math.pi**2 + np.square(np.log(self.e[:]))))
        self.alpha *= math.sqrt(5.0/6.0)

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
    def gravity(self):

        self.a[:,1] -= self.m[:]*self.g

    ### ************************************************
    ### Update positions using verlet method
    def update(self, dt):

        self.v[:] += 0.5*dt*(self.a[:] + self.f[:]/self.m[:])
        self.d[:]  = dt*self.v[:] + 0.5*dt*dt*self.f[:]/self.m[:]
        self.x[:] += self.d[:]
        self.a[:]  = self.f[:]/self.m[:]

        if (self.store):
            self.history = np.vstack((self.history, self.x.reshape((1,-1))))
