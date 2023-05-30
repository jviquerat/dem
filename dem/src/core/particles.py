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
                 store       = False,
                 search      = "linear",
                 rad_coeff   = 5.0):

        self.np          = np
        self.nt          = nt
        self.mtr         = material_factory.create(material)
        self.radius      = radius
        self.vol         = (4.0/3.0)*math.pi*radius**3
        self.mass        = self.vol*self.mtr.density
        self.store       = store
        self.color       = color
        self.search      = search
        self.rad_coeff   = rad_coeff

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
        self.ngb = None                                          # neighbors

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
    ### Compute maximal radius
    def max_radius(self):

        return np.max(self.r)

    ### ************************************************
    ### Compute collisions between particles
    def collisions(self, dt):

        if (self.search == "linear"):

            # Compute list of collisions
            ci, cj, cd = linear_search(self.np, self.x, self.r, 0.0)

            # Check if there are collisions
            n_coll = len(ci)
            if (n_coll == 0): return

            # Compute forces
            collide(self, dt, cd, ci, cj, n_coll)

        if (self.search == "nearest"):

            # If the list is not set yet
            #if (self.ngb == None):
            self.reset_ngb()

            # Update if list is not valid anymore

            # Compute collisions
            for i in range(self.np):
                for j in self.ngb[i]:
                    collide_single(self, dt, i, j)

    ### ************************************************
    ### Reset neighbor particles
    def reset_ngb(self):

        self.ngb = [None]*self.np
        for i in range(self.np):
            self.ngb[i] = np.empty((0), np.uint16)

        r_ref      = self.rad_coeff*self.max_radius()
        ci, cj, cd = linear_search(self.np, self.x, self.r, r_ref)

        for k in range(len(ci)):
            i = ci[k]
            j = cj[k]
            self.ngb[i] = np.append(self.ngb[i], np.uint(j))
            self.ngb[j] = np.append(self.ngb[j], np.uint(i))

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
### to compute collisions. r_ref is a reference radius
### that is equal to 0 for regular linear search detection,
### and that is equal to alpha*max_radius when generating
### nearest neighbor lists
@nb.njit(cache=True)
def linear_search(n, x, r, r_ref):

    ci = np.empty((0), np.uint16)
    cj = np.empty((0), np.uint16)
    cd = np.empty((0), np.float32)

    for i in range(n):
        for j in range(i+1,n):
            dist = (x[i,0]-x[j,0])**2 + (x[i,1]-x[j,1])**2
            dist = math.sqrt(dist)
            dx   = dist - r[i] - r[j]

            if (dx < r_ref):
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

### ************************************************
### Compute collisions between particles
#@nb.njit(cache=True)
def collide_single(p, dt, i, j):

    dx = (p.x[i,0]-p.x[j,0])**2 + (p.x[i,1]-p.x[j,1])**2
    dx = math.sqrt(dx)
    if (dx - p.r[i] - p.r[j] > 0.0): return
    dx = abs(dx)

    # Compute normal and tangent
    n    = np.zeros(2)
    x_ij = p.x[j,:] - p.x[i,:]
    n    = x_ij/(dx + p.r[i] + p.r[j])

    # Return forces from collision parameters
    # - normal elastic
    # - normal damping,
    # - tangential elastic
    # - tangential damping
    fn, ft = hertz(dx,            # penetration
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
    #p.f[j,:] += fn[:]

    # tangential force
    p.f[i,:] -= ft[:]
    #p.f[j,:] += ft[:]
