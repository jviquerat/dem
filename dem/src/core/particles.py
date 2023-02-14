# Generic imports
import math
import numpy as np

# Custom imports
from dem.src.core.collision import *

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

        self.m = np.ones((self.n),    np.float32)*self.mass        # masses
        self.r = np.ones((self.n),    np.float32)*self.radius      # radii
        self.x = np.zeros((self.n,2), np.float32)                  # positions
        self.d = np.zeros((self.n,2), np.float32)                  # displacements
        self.v = np.zeros((self.n,2), np.float32)                  # velocities
        self.a = np.zeros((self.n,2), np.float32)                  # accelerations
        self.f = np.zeros((self.n,2), np.float32)                  # forces
        self.e = np.ones((self.n),    np.float32)*self.restitution # restitution coeff
        self.y = np.ones((self.n),    np.float32)*self.young       # young modulus
        self.p = np.ones((self.n),    np.float32)*self.poisson     # poisson ratio
        self.c = [self.color]*self.n                               # colors

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

        ci, cj, cd = particles_distances(self.n, self.x, self.r)
        particles_collisions(ci, cj, cd, self.x, self.r, self.m,
                             self.v, self.g, self.Eb, self.Gb, self.f)

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

### ************************************************
### Compute inter-particles distances and return those
### which need to be accounted for
@nb.njit(cache=True)
def particles_distances(n, x, r):

    ci = np.empty((0), np.uint16)
    cj = np.empty((0), np.uint16)
    cd = np.empty((0), np.float32)

    for i in range(n):
        for j in range(i+1,n):
            dist = (x[i,0]-x[j,0])**2 + (x[i,1]-x[j,1])**2
            dist = math.sqrt(dist)
            dx   = dist - r[i] - r[j]

            if (dx < 0.0):
                ci = np.append(ci, np.uint16(i))
                cj = np.append(cj, np.uint16(j))
                cd = np.append(cd, np.float32(abs(dx)))

    return ci, cj, cd

### ************************************************
### Compute collisions between particles
@nb.njit(cache=True)
def particles_collisions(ci, cj, cd, x, r, m, v, g, Eb, Gb, f):

    n_coll = len(ci)

    for k in range(n_coll):
        i = ci[k]
        j = cj[k]

        # Compute normal and tangent
        x_ij = x[j,:] - x[i,:]
        n    = x_ij/(cd[k] + r[i] + r[j])
        t    = np.zeros(2)
        t[0] =-n[1]
        t[1] = n[0]

        # Return forces from collision parameters
        # - normal elastic
        # - normal damping,
        # - tangential elastic
        # - tangential damping
        fne, fnd, fte, ftd = hertz(cd[k],          # penetration
                                   r[i],   r[j],   # radii
                                   m[i],   m[j],   # masses
                                   v[i,:], v[j,:], # velocities
                                   n,      t,      # normal and tangent
                                   g[i],   g[j],   # restitution coeffs
                                   Eb[i],  Eb[j],  # Eb coeffs
                                   Gb[i],  Gb[j])  # Gb coeffs

        # normal elastic force
        f[i,:] -= fne[:]
        f[j,:] += fne[:]

        # normal damping force
        f[i,:] -= fnd[:]
        f[j,:] += fnd[:]

        # # tangential elastic force
        # f[i,:] += fte[:]
        # f[j,:] -= fte[:]

        # # tangential damping force
        # f[i,:] += ftd[:]
        # f[j,:] -= ftd[:]
