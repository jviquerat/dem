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
                 rad_coeff   = 1.0):

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
        self.ngbs        = None
        self.d_ngb       = None

        self.reset()

    ### ************************************************
    ### Reset arrays
    def reset(self):

        self.m = np.ones( (self.np), )*self.mass   # masses
        self.r = np.ones( (self.np), )*self.radius # radii
        self.x = np.zeros((self.np,2))             # positions
        self.d = np.zeros((self.np,2))             # displacements
        self.v = np.zeros((self.np,2))             # velocities
        self.a = np.zeros((self.np,2))             # accelerations
        self.f = np.zeros((self.np,2))             # forces

        # default material parameters
        self.e_wall  = np.ones((self.np))*self.mtr.e_wall
        self.mu_wall = np.ones((self.np))*self.mtr.mu_wall
        self.e_part  = np.ones((self.np))*self.mtr.e_part
        self.mu_part = np.ones((self.np))*self.mtr.mu_part
        self.Y       = np.ones((self.np))*self.mtr.Y
        self.G       = np.ones((self.np))*self.mtr.G

        self.c       = [self.color]*self.np # colors
        self.m_rad   = self.rad_coeff*self.max_radius()

        # Optional storage
        if self.store:
            self.names = ["x", "y", "vx", "vy"]
            self.buff  = buff(self.np, self.nt, self.names)

    ### ************************************************
    ### Reset forces
    def reset_forces(self):

        self.f[:,:] = 0.0

    ### ************************************************
    ### Set particle i with material m
    def set_material(self, i, m):

        self.e_wall[i]  = m.e_wall
        self.mu_wall[i] = m.mu_wall
        self.e_part[i]  = m.e_part
        self.mu_part[i] = m.mu_part
        self.Y[i]       = m.Y
        self.G[i]       = m.G

    ### ************************************************
    ### Compute maximal velocity
    def max_velocity(self):

        if (len(self.v) == 0):
            return 0.0
        else:
            return np.max(np.linalg.norm(self.v, axis=1))

    ### ************************************************
    ### Compute maximal radius
    def max_radius(self):

        if (len(self.r) == 0):
            return 0.0
        else:
            return np.max(self.r)

    ### ************************************************
    ### Remove a list of particles
    def delete(self, lst):

        self.np     -= len(lst)
        self.m       = np.delete(self.m,       lst, axis=0)
        self.r       = np.delete(self.r,       lst, axis=0)
        self.x       = np.delete(self.x,       lst, axis=0)
        self.d       = np.delete(self.d,       lst, axis=0)
        self.v       = np.delete(self.v,       lst, axis=0)
        self.a       = np.delete(self.a,       lst, axis=0)
        self.f       = np.delete(self.f,       lst, axis=0)
        self.e_wall  = np.delete(self.e_wall,  lst, axis=0)
        self.mu_wall = np.delete(self.mu_wall, lst, axis=0)
        self.e_part  = np.delete(self.e_part,  lst, axis=0)
        self.mu_part = np.delete(self.mu_part, lst, axis=0)
        self.Y       = np.delete(self.Y,       lst, axis=0)
        self.G       = np.delete(self.G,       lst, axis=0)
        self.c       = np.delete(np.array(self.c), lst, axis=0).tolist()

    ### ************************************************
    ### Compute collisions between particles
    def collisions(self, dt, force_nearest=False):

        if (self.search == "linear"):

            # Search for collisions linearly and compute forces
            linear_search(self.x,      self.r,  self.m, self.v,
                          self.e_part, self.Y,  self.G, self.mu_part,
                          self.f,      self.np, dt)

        if (self.search == "nearest"):

            # If the list of neighbors is not set yet
            if (self.d_ngb is None):
                self.d_ngb = 0.0
                self.m_rad = self.rad_coeff*self.max_radius()
                self.ngbi, self.ngbj = list_ngbs(self.np, self.x,
                                                 self.r,  self.m_rad)

            # If list is not valid anymore
            v_max       = self.max_velocity()
            self.d_ngb += 2.0*v_max*dt
            if ((self.d_ngb > 0.99*self.m_rad) or force_nearest):
                self.d_ngb = 0.0
                self.ngbi, self.ngbj = list_ngbs(self.np, self.x,
                                                 self.r,  self.m_rad)

            # Check if there are collisions
            n_coll = len(self.ngbi)
            if (n_coll == 0): return

            # Compute forces
            collide_ngbs(self.x,      self.r,    self.m,    self.v,
                         self.e_part, self.Y,    self.G,    self.mu_part,
                         self.f,      self.ngbi, self.ngbj, n_coll, dt)

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
def linear_search(x, r, m, v, e_part, Y, G, mu_part, f, n, dt):

    # Loop on particles twice
    for i in range(n):
        for j in range(i+1,n):
            dx = (x[i,0]-x[j,0])**2 + (x[i,1]-x[j,1])**2
            dx = math.sqrt(dx)
            dx = dx - r[i] - r[j]

            # If particles intersect
            if (dx < 0.0):

                # Compute normal
                dx   = abs(dx)
                nrm  = np.zeros(2)
                x_ij = x[j,:] - x[i,:]
                nrm  = x_ij/(dx + r[i] + r[j])

                # Return forces from collision parameters
                # - normal elastic
                # - normal damping,
                # - tangential elastic
                # - tangential damping
                fn, ft = hertz(dx,         # penetration
                               dt,         # timestep
                               r[i],       # radius 1
                               r[j],       # radius 2
                               m[i],       # mass 1
                               m[j],       # mass 2
                               v[i,:],     # velocity 1
                               v[j,:],     # velocity 2
                               nrm[:],     # normal from 1 to 2
                               e_part[i],  # restitution 1
                               e_part[j],  # restitution 2
                               Y[i],       # effective young modulus 1
                               Y[j],       # effective young modulus 2
                               G[i],       # effective shear modulus 1
                               G[j],       # effective shear modulus 2
                               mu_part[i], # static friction 1
                               mu_part[j]) # static friction 2

                # normal force
                f[i,:] -= fn[:]
                f[j,:] += fn[:]

                # tangential force
                f[i,:] -= ft[:]
                f[j,:] += ft[:]

### ************************************************
### Compute inter-particles distances and return data
### to compute collisions. r_ref is a reference radius
### that is equal to 0 for regular linear search detection,
### and that is equal to alpha*max_radius when generating
### nearest neighbor lists
@nb.njit(cache=True)
def list_ngbs(n, x, r, r_ref):

    ngbi = np.empty((0), np.uint16)
    ngbj = np.empty((0), np.uint16)
    for i in range(n):
        for j in range(i+1,n):
            dist = (x[i,0]-x[j,0])**2 + (x[i,1]-x[j,1])**2
            dist = math.sqrt(dist)
            dx   = dist - r[i] - r[j]

            if (dx < r_ref):
                ngbi = np.append(ngbi, np.uint16(i))
                ngbj = np.append(ngbj, np.uint16(j))

    return ngbi, ngbj

### ************************************************
### Compute inter-particles distances and return data
### to compute collisions. r_ref is a reference radius
### that is equal to 0 for regular linear search detection,
### and that is equal to alpha*max_radius when generating
### nearest neighbor lists
@nb.njit(cache=True)
def collide_ngbs(x, r, m, v, e_part, Y, G, mu_part, f, ngbi, ngbj, n, dt):

    for k in range(n):
        i = ngbi[k]
        j = ngbj[k]

        dx = (x[i,0]-x[j,0])**2 + (x[i,1]-x[j,1])**2
        dx = math.sqrt(dx)
        dx = dx - r[i] - r[j]

        # If particles intersect
        if (dx < 0.0):

            # Compute normal
            dx   = abs(dx)
            nrm  = np.zeros(2)
            x_ij = x[j,:] - x[i,:]
            nrm  = x_ij/(dx + r[i] + r[j])

            # Return forces from collision parameters
            # - normal elastic
            # - normal damping,
            # - tangential elastic
            # - tangential damping
            fn, ft = hertz(dx,         # penetration
                           dt,         # timestep
                           r[i],       # radius 1
                           r[j],       # radius 2
                           m[i],       # mass 1
                           m[j],       # mass 2
                           v[i,:],     # velocity 1
                           v[j,:],     # velocity 2
                           nrm[:],     # normal from 1 to 2
                           e_part[i],  # restitution 1
                           e_part[j],  # restitution 2
                           Y[i],       # effective young modulus 1
                           Y[j],       # effective young modulus 2
                           G[i],       # effective shear modulus 1
                           G[j],       # effective shear modulus 2
                           mu_part[i], # static friction 1
                           mu_part[j]) # static friction 2

            # normal force
            f[i,:] -= fn[:]
            f[j,:] += fn[:]

            # tangential force
            f[i,:] -= ft[:]
            f[j,:] += ft[:]
