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
        self.v     = np.zeros((2))

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

        # Search for collisions linearly and compute forces
        linear_search(self.a, self.b, self.c, self.d, self.r, self.m,
                      self.v, self.mat.Y, self.mat.G, self.n,
                      p.x, p.r, p.m, p.v, p.e_wall, p.mu_wall,
                      p.Y, p.G, p.f, p.np, dt)

### ************************************************
### Distance from rectangle domain to given coordinates
### Prefix d_ corresponds to domain
### Prefix p_ corresponds to particle
@nb.njit(cache=True)
def linear_search(d_a, d_b, d_c, d_d, d_r, d_m,
                  d_v, d_mat_Y, d_mat_G, d_n,
                  p_x, p_r, p_m, p_v, p_e_wall, p_mu_wall,
                  p_Y, p_G, p_f, n, dt):

    # Loop on particles
    for i in range(n):

        # Loop on rectangle sides
        for j in range(4):
            dx = abs(d_a[j]*p_x[i,0] + d_b[j]*p_x[i,1] + d_c[j])/d_d[j]
            dx = dx - p_r[i]

            # If particle intersects boundary
            if (dx < 0.0):

                # Compute normal
                dx     = abs(dx)
                nrm    = np.zeros(2, np.float32)
                nrm[:] =-d_n[j,:]

                # Return forces from collision parameters
                # - normal elastic
                # - normal damping,
                # - tangential elastic
                # - tangential damping
                fn, ft = hertz(dx,           # penetration
                               dt,           # timestep
                               p_r[i],       # radius 1
                               d_r,          # radius 2
                               p_m[i],       # mass 1
                               d_m,          # mass 2
                               p_v[i,:],     # velocity 1
                               d_v,          # velocity 2
                               nrm[:],       # normal from 1 to 2
                               p_e_wall[i],  # restitution 1
                               p_e_wall[i],  # restitution 2
                               p_Y[i],       # effective young modulus 1
                               d_mat_Y,      # effective young modulus 2
                               p_G[i],       # effective shear modulus 1
                               d_mat_G,      # effective shear modulus 2
                               p_mu_wall[i], # static friction 1
                               p_mu_wall[i]) # static friction 2

                # normal force
                p_f[i,:] -= fn[:]

                # tangential force
                p_f[i,:] -= ft[:]
