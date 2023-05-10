# Generic imports
import math
import numpy as np

# Custom imports
from dem.src.core.collision    import *
from dem.src.material.material import *
from dem.src.utils.buff        import *

### ************************************************
### Class defining array of particles
class particles:
    ### ************************************************
    ### Constructor
    def __init__(self,
                 np          = 10,
                 nt          = 1000,
                 material    = "steel",
                 radius      = 1.0,
                 color       = "b",
                 store       = False):

        self.np          = np
        self.nt          = nt
        self.mtr         = material_factory.create(material)
        self.radius      = radius
        self.vol         = (4.0/3.0)*math.pi*radius**3
        self.mass        = self.vol*self.mtr.density
        self.store       = store
        self.color       = color

        self.reset()

    ### ************************************************
    ### Reset arrays
    def reset(self):

        self.m   = np.ones( (self.np),   np.float32)*self.mass   # masses
        self.r   = np.ones( (self.np),   np.float32)*self.radius # radii
        self.x   = np.zeros((self.np,2), np.float32)             # positions
        self.d   = np.zeros((self.np,2), np.float32)             # displacements
        self.v   = np.zeros((self.np,2), np.float32)             # velocities
        self.a   = np.zeros((self.np,2), np.float32)             # accelerations
        self.f   = np.zeros((self.np,2), np.float32)             # forces
        self.mat = [self.mtr  ]*self.np                          # materials
        self.c   = [self.color]*self.np                          # colors

        # Optional storage
        if self.store:
            self.names = ["x", "y", "vx", "vy"]
            self.buff  = buff(self.np, self.nt, self.names)

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
    def collisions(self, dt):

        # Compute list of collisions
        ci, cj, cd = list_collisions(self.np, self.x, self.r)

        # Check if there are collisions
        n_coll = len(ci)
        if (n_coll == 0): return

        # Compute forces
        collide(self, dt, cd, ci, cj, n_coll)

    ### ************************************************
    ### Add gravity
    def gravity(self, g):

        self.f[:,1] -= self.m[:]*g

    ### ************************************************
    ### Update positions using verlet method
    def update(self, t, dt, it):

        if (self.store):
            self.buff.store( t, self.names,
                            [self.x[:,0], self.x[:,1],
                             self.v[:,0], self.v[:,1]])

        self.f[:,0] /= self.m[:]
        self.f[:,1] /= self.m[:]
        self.v[:,:] += 0.5*dt*(self.a[:,:] + self.f[:,:])
        self.d[:,:]  = dt*self.v[:,:] + 0.5*dt*dt*self.f[:,:]
        self.x[:,:] += self.d[:,:]
        self.a[:,:]  = self.f[:,:]

### ************************************************
### Compute inter-particles distances and return data
### to compute collisions
@nb.njit(cache=True)
def list_collisions(n, x, r):

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
#@nb.njit(cache=True)
def collide(p, dt, dx, ci, cj, n_coll):

    for k in range(n_coll):
        i = ci[k]
        j = cj[k]

        # Compute normal and tangent
        n    = np.zeros(2)
        x_ij = p.x[j,:] - p.x[i,:]
        n    = x_ij/(dx[k] + p.r[i] + p.r[j])

        # Return forces from collision parameters
        # - normal elastic
        # - normal damping,
        # - tangential elastic
        # - tangential damping
        fn, ft = hertz(dx[k],            # penetration
                       dt,               # timestep
                       p.r[i],           # radius 1
                       p.r[j],           # radius 2
                       p.m[i],           # mass 1
                       p.m[j],           # mass 2
                       p.v[i,:],         # velocity 1
                       p.v[j,:],         # velocity 2
                       n[:],             # normal from 1 to 2
                       p.mat[i].e_part,  # restitution 1
                       p.mat[j].e_part,  # restitution 2
                       p.mat[i].Y,       # effective young modulus 1
                       p.mat[j].Y,       # effective young modulus 2
                       p.mat[i].G,       # effective shear modulus 1
                       p.mat[j].G,       # effective shear modulus 2
                       p.mat[i].mu_part, # static friction 1
                       p.mat[j].mu_part) # static friction 2

        # normal force
        p.f[i,:] -= fn[:]
        p.f[j,:] += fn[:]

        # tangential force
        p.f[i,:] -= ft[:]
        p.f[j,:] += ft[:]
