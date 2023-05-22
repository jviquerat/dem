# Generic imports
import matplotlib
import numpy as np
import numba as nb

from matplotlib.patches import Rectangle

# Custom imports
from dem.src.domain.base       import *
from dem.src.core.collision    import *
from dem.src.material.material import *

### ************************************************
### Class defining rectangle domain
class rectangle(base_domain):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 x_min    = 0.0,
                 x_max    = 1.0,
                 y_min    = 0.0,
                 y_max    = 1.0,
                 material = "steel"):

        # External boundaries
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        # Material
        self.mat = material_factory.create(material)

        # Ficticious parameters for domain
        self.r     = 1.0e8
        self.m     = 1.0e8
        self.g     = 1.0
        self.v     = np.zeros((2))
        self.d     = np.zeros((2))

        # Define a*x + b*y + c = 0 for all four borders
        # Ridge 0 is the bottom one, then we pursue in
        # trigonometric order. d array corresponds to the
        # normalization parameter sqrt(a*a + b*b).
        # n and t are the normal and tangential vectors
        self.a = np.array([0.0, 1.0, 0.0, 1.0])
        self.b = np.array([1.0, 0.0, 1.0, 0.0])
        self.c = np.array([-self.y_min,-self.x_max,
                           -self.y_max,-self.x_min])
        self.d = np.array([1.0, 1.0, 1.0, 1.0])

        # Inward normals
        self.n = np.array([[ 0.0, 1.0],
                           [-1.0, 0.0],
                           [ 0.0,-1.0],
                           [ 1.0, 0.0]], np.float32)

        # Trigonometric direction tangents
        self.t = np.array([[ 1.0, 0.0],
                           [ 0.0, 1.0],
                           [-1.0, 0.0],
                           [ 0.0,-1.0]], np.float32)

    ### ************************************************
    ### Plot domain
    def plot(self, ax):

        ax.add_patch(Rectangle((self.x_min, self.y_min),
                                self.x_max-self.x_min,
                                self.y_max-self.y_min,
                                fill=False, color='r'))

    ### ************************************************
    ### Compute collisions with a set of particles
    def collisions(self, p, dt):

        # Compute distances to domain boundaries
        ci, cj, cd = rectangle_distance(self.a, self.b, self.c, self.d,
                                        p.x, p.r, p.np)

        # Check if there are collisions
        n_coll = len(ci)
        if (n_coll == 0): return

        # Compute forces
        self.collide(p, dt, cd, ci, cj, n_coll)

    ### ************************************************
    ### Collision forces on particle in rectangle domain
    def collide(self, p, dt, dx, ci, cj, n_coll):

        # Loop on collisions with domain
        for k in range(n_coll):
            i = ci[k]
            j = cj[k]

            # Compute normal and tangent
            n    = np.zeros(2, np.float32)
            n[:] =-self.n[j,:]

            # Return forces from collision parameters
            # - normal elastic
            # - normal damping,
            # - tangential elastic
            # - tangential damping
            fn, ft = hertz(dx[k],            # penetration
                           dt,               # timestep
                           p.r[i],           # radius 1
                           self.r,           # radius 2
                           p.m[i],           # mass 1
                           self.m,           # mass 2
                           p.v[i,:],         # velocity 1
                           self.v,           # velocity 2
                           n[:],             # normal from 1 to 2
                           p.mat[i].e_wall,  # restitution 1
                           p.mat[i].e_wall,  # restitution 2
                           p.mat[i].Y,       # effective young modulus 1
                           self.mat.Y,       # effective young modulus 2
                           p.mat[i].G,       # effective shear modulus 1
                           self.mat.G,       # effective shear modulus 2
                           p.mat[i].mu_wall, # static friction 1
                           p.mat[i].mu_wall) # static friction 2

            # normal force
            p.f[i,:] -= fn[:]

            # tangential force
            p.f[i,:] -= ft[:]

### ************************************************
### Distance from rectangle domain to given coordinates
@nb.njit(cache=True)
def rectangle_distance(a, b, c, d, x, r, n):

    ci = np.empty((0), np.uint16)
    cj = np.empty((0), np.uint16)
    cd = np.empty((0), np.float32)

    for i in range(n):
        for j in range(4):
            dist = abs(a[j]*x[i,0] + b[j]*x[i,1] + c[j])/d[j]
            dx   = dist - r[i]

            if (dx < 0.0):
                ci = np.append(ci, np.uint16(i))
                cj = np.append(cj, np.uint16(j))
                cd = np.append(cd, np.float32(abs(dx)))

    return ci, cj, cd
