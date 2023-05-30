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
                 rad_coeff   = 2.0):

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

        self.m = np.ones( (self.np),   np.float32)*self.mass   # masses
        self.r = np.ones( (self.np),   np.float32)*self.radius # radii
        self.x = np.zeros((self.np,2), np.float32)             # positions
        self.d = np.zeros((self.np,2), np.float32)             # displacements
        self.v = np.zeros((self.np,2), np.float32)             # velocities
        self.a = np.zeros((self.np,2), np.float32)             # accelerations
        self.f = np.zeros((self.np,2), np.float32)             # forces

        # default material parameters
        self.e_wall  = np.ones((self.np), np.float32)*self.mtr.e_wall
        self.mu_wall = np.ones((self.np), np.float32)*self.mtr.mu_wall
        self.e_part  = np.ones((self.np), np.float32)*self.mtr.e_part
        self.mu_part = np.ones((self.np), np.float32)*self.mtr.mu_part
        self.Y       = np.ones((self.np), np.float32)*self.mtr.Y
        self.G       = np.ones((self.np), np.float32)*self.mtr.G

        self.c       = [self.color]*self.np # colors

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
            linear_search(self.x,      self.r,  self.m, self.v,
                          self.e_part, self.Y,  self.G, self.mu_part,
                          self.f,      self.np, dt)

            # # Check if there are collisions
            # n_coll = len(ci)
            # if (n_coll == 0): return

            # # Compute forces
            # collide(self, dt, ci, cj, n_coll)

        if (self.search == "nearest"):

            pass

            # # If the list is not set yet
            # if (self.ngb == None):
            #     self.d_ngb = 0.0
            #     self.m_rad = self.rad_coeff*self.max_radius()
            #     ci, cj, cd = linear_search(self.np, self.x, self.r, self.m_rad)

            # # If list is not valid anymore
            # v_max = self.max_velocity()
            # self.d_ngb += 2.0*v_max*dt
            # if (self.d_ngb > 0.95*self.m_rad):
            #     ci, cj = linear_search(self.np, self.x, self.r, self.m_rad)
            #     self.d_ngb = 0.0

            # # Check if there are collisions
            # n_coll = len(ci)
            # if (n_coll == 0): return

            # # Compute forces
            # collide(self, dt, ci, cj, n_coll)

            # # # Compute collisions
            # # for i in range(self.np):
            # #     for j in self.ngb[i]:
            # #         dist = (self.x[i,0]-self.x[j,0])**2 + (self.x[i,1]-self.x[j,1])**2
            # #         dist = math.sqrt(dist)
            # #         dx   = dist - self.r[i] - self.r[j]
            # #         if (dx < 0.0):
            # #             dx = abs(dx)
            # #             collide_single(self, dx, dt, i, j)

    ### ************************************************
    ### Reset neighbor particles
    def reset_ngb(self):

        self.ngb = [None]*self.np
        for i in range(self.np):
            self.ngb[i] = np.empty((0), np.uint16)

        ci, cj, cd = linear_search(self.np, self.x, self.r, self.m_rad)

        for k in range(len(ci)):
            i = ci[k]
            j = cj[k]
            if (i > j):
                self.ngb[i] = np.append(self.ngb[i], np.uint(j))
            else:
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
#@nb.njit(cache=True)
# def linear_search(n, x, r, r_ref):

#     #ci = np.empty((0), np.uint16)
#     #cj = np.empty((0), np.uint16)
#     #cd = np.empty((0), np.float32)

#     for i in range(n):
#         for j in range(i+1,n):
#             dist = (x[i,0]-x[j,0])**2 + (x[i,1]-x[j,1])**2
#             dist = math.sqrt(dist)
#             dx   = dist - r[i] - r[j]

#             if (dx < r_ref):
#                 ci = np.append(ci, np.uint16(i))
#                 cj = np.append(cj, np.uint16(j))
#                 #cd = np.append(cd, np.float32(abs(dx)))

#     return ci, cj

### ************************************************
### Compute collisions between particles
#@nb.njit(cache=True)
def collide(p, dt, ci, cj, n_coll):

    for k in range(n_coll):
        i = ci[k]
        j = cj[k]

        dist = (p.x[i,0]-p.x[j,0])**2 + (p.x[i,1]-p.x[j,1])**2
        dist = math.sqrt(dist)
        dx   = abs(dist - p.r[i] - p.r[j])

        # Compute normal and tangent
        n    = np.zeros(2)
        x_ij = p.x[j,:] - p.x[i,:]
        n    = x_ij/(dx + p.r[i] + p.r[j])

        # Return forces from collision parameters
        # - normal elastic
        # - normal damping,
        # - tangential elastic
        # - tangential damping
        fn, ft = hertz(dx,               # penetration
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
def collide_single(p, dx, dt, i, j):

    # dx = (p.x[i,0]-p.x[j,0])**2 + (p.x[i,1]-p.x[j,1])**2
    # dx = math.sqrt(dx)
    # if (dx - p.r[i] - p.r[j] > 0.0): return
    # dx = abs(dx)

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
    p.f[j,:] += fn[:]

    # tangential force
    p.f[i,:] -= ft[:]
    p.f[j,:] += ft[:]
