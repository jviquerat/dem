# Generic imports
import math
import numpy as np

# Custom imports
from dem.src.core.collision    import *
from dem.src.core.grid         import *
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
        self.search      = "linear"
        self.grid        = None

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
    ### Set grid
    def set_grid(self, domain, nx, ny):

        if (self.search == "grid"):
            self.grid = grid(domain, nx, ny)
            self.grid.set(self.x)

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
    ### Compute maximal radius
    def max_radius(self):

        return np.max(self.r)

    ### ************************************************
    ### Compute collisions between particles
    def collisions(self, dt):

        if ((self.grid is None) and (self.search == "grid")):
            self.grid.set(self.x)

        # Linear search type
        if (self.search == "linear"):

            idx = np.arange(0, self.np, dtype=int)
            ci, cj, cd = linear_search(self.np, idx, self.x, self.r)

            n_coll = len(ci)
            if (n_coll == 0): return

            self.collide(dt, cd, ci, cj)

        # Grid search type
        if (self.search == "grid"):

            # Correct grid
            self.grid.correct(self.x)

            # Set checked array
            checked = np.zeros((self.grid.n, self.grid.n), bool)

            # Loop on grid cells
            for i in range(self.grid.nx):
                for j in range(self.grid.ny):

                    # Loop on neighbors (including local cell)
                    ngb = self.grid.ngb[i,j]
                    for n in range(len(ngb)):

                        # Retrieve indices of current neighbor and id
                        k = ngb[n][0]
                        l = ngb[n][1]

                        # Check if interactions were already computed
                        id_ij = self.grid.id(i,j)
                        id_kl = self.grid.id(k,l)
                        if (checked[id_ij, id_kl]): continue

                        # Retrieve particles in neighbor cellOA
                        parts = self.grid.parts[k,l]
                        if (len(parts) == 0): continue

                        # Loop on local particles
                        for p in self.grid.parts[i,j]:

                            # Check collisions
                            ci, cj, cd = single_search(p,
                                                       self.x[p,:],
                                                       self.r[p],
                                                       len(parts), parts,
                                                       self.x[parts,:],
                                                       self.r[parts])

                            n_coll = len(ci)
                            if (n_coll == 0): continue

                            # Compute forces
                            self.collide(dt, cd, ci, cj)

                        # Note computed interactions
                        checked[id_ij, id_kl] = True
                        checked[id_kl, id_ij] = True

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
    ### Compute collisions between particles
    def collide(self, dt, dx, ci, cj):

        n_coll = len(ci)

        for k in range(n_coll):
            i = ci[k]
            j = cj[k]

            if (i == j): continue

            # Compute normal and tangent
            n    = np.zeros(2)
            x_ij = self.x[j,:] - self.x[i,:]
            n    = x_ij/(dx[k] + self.r[i] + self.r[j])

            # Return forces from collision parameters
            # - normal elastic
            # - normal damping,
            # - tangential elastic
            # - tangential damping
            fn, ft = hertz(dx[k],               # penetration
                           dt,                  # timestep
                           self.r[i],           # radius 1
                           self.r[j],           # radius 2
                           self.m[i],           # mass 1
                           self.m[j],           # mass 2
                           self.v[i,:],         # velocity 1
                           self.v[j,:],         # velocity 2
                           n[:],                # normal from 1 to 2
                           self.mat[i].e_part,  # restitution 1
                           self.mat[j].e_part,  # restitution 2
                           self.mat[i].Y,       # effective young modulus 1
                           self.mat[j].Y,       # effective young modulus 2
                           self.mat[i].G,       # effective shear modulus 1
                           self.mat[j].G,       # effective shear modulus 2
                           self.mat[i].mu_part, # static friction 1
                           self.mat[j].mu_part) # static friction 2

            # normal force
            self.f[i,:] -= fn[:]
            self.f[j,:] += fn[:]

            # tangential force
            self.f[i,:] -= ft[:]
            self.f[j,:] += ft[:]

### ************************************************
### Linear search for collisions (all particles together)
@nb.njit(cache=True)
def linear_search(n, idx, x, r):

    # Set arrays
    ci = np.empty((0), np.uint16)
    cj = np.empty((0), np.uint16)
    cd = np.empty((0), np.float32)

    for ii in range(n):
        for jj in range(ii+1,n):
            i = idx[ii]
            j = idx[jj]

            dist = (x[i,0]-x[j,0])**2 + (x[i,1]-x[j,1])**2
            dist = math.sqrt(dist)
            dx   = dist - r[i] - r[j]

            if (dx < 0.0):
                ci = np.append(ci, np.uint16(i))
                cj = np.append(cj, np.uint16(j))
                cd = np.append(cd, np.float32(abs(dx)))

    return ci, cj, cd

### ************************************************
### Linear search for collisions (one particle against all others)
@nb.njit(cache=True)
def single_search(i, x0, r0, n, idx, x, r):

    # Set arrays
    ci = np.empty((0), np.uint16)
    cj = np.empty((0), np.uint16)
    cd = np.empty((0), np.float32)

    for jj in range(n):
        j = idx[jj]
        if (i==j): continue

        dist = (x0[0]-x[jj,0])**2 + (x0[1]-x[jj,1])**2
        dist = math.sqrt(dist)
        dx   = dist - r0 - r[jj]

        if (dx < 0.0):
            ci = np.append(ci, np.uint16(i))
            cj = np.append(cj, np.uint16(j))
            cd = np.append(cd, np.float32(abs(dx)))

    return ci, cj, cd
