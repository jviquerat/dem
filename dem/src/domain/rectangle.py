# Generic imports
import matplotlib
import numpy as np
import numba as nb

from matplotlib.patches import Rectangle

# Custom imports
from dem.src.domain.base_domain import *
from dem.src.core.collision     import *

### ************************************************
### Class defining rectangle domain
class rectangle(base_domain):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 x_min   = 0.0,
                 x_max   = 1.0,
                 y_min   = 0.0,
                 y_max   = 1.0,
                 young   = 210.0e9,
                 poisson = 0.25):

        # External boundaries
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        # Material (steel)
        self.y  = young   # young modulus
        self.p  = poisson # poisson ratio
        self.Eb = np.ones((4))*(1.0 - self.p**2)/self.y
        self.Gb = np.ones((4))*2.0*(2.0 + self.p)*(1.0 - self.p)/self.y

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
        self.n = np.array([[ 0.0, 1.0],
                           [-1.0, 0.0],
                           [ 0.0,-1.0],
                           [ 1.0, 0.0]])
        self.t = np.array([[ 1.0, 0.0],
                           [ 0.0, 1.0],
                           [-1.0, 0.0],
                           [ 0.0,-1.0]])

    ### ************************************************
    ### Plot domain
    def plot(self, ax):

        ax.add_patch(Rectangle((self.x_min, self.y_min),
                                self.x_max-self.x_min,
                                self.y_max-self.y_min,
                                fill=False, color='r'))

    ### ************************************************
    ### Compute collisions with a set of particles
    def collisions(self, p):

        # Compute distances to domain boundaries
        ci, cj, cd = rectangle_distance(self.a, self.b, self.c, self.d,
                                        p.x, p.r, p.n)

        # Check if there are collisions
        n_coll = len(ci)
        if (n_coll == 0): return

        # Compute forces
        rectangle_forces(p.f, p.r, p.m, p.v, p.g, p.Eb, p.Gb, cd,
                         self.Eb, self.Gb, self.n, self.t, ci, cj, n_coll)

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

### ************************************************
### Collision forces on particle in rectangle domain
@nb.njit(cache=True)
def rectangle_forces(f, p_r, p_m, p_v, p_g, p_Eb, p_Gb, dx,
                     d_Eb, d_Gb, n, t, ci, cj, n_coll):

    # Ficticious parameters for domain
    d_r     = 1.0e8
    d_m     = 1.0e8
    d_g     = 1.0
    d_v     = np.zeros((2))

    # Loop on collisions with domain
    for k in range(n_coll):
        i = ci[k]
        j = cj[k]

        # Return forces from collision parameters
        # - normal elastic
        # - normal damping,
        # - tangential elastic
        # - tangential damping
        fne, fnd, fte, ftd = hertz(dx[k],             # penetration
                                   p_r[i],   d_r,     # radii
                                   p_m[i],   d_m,     # masses
                                   p_v[i,:], d_v,     # velocities
                                   n[j,:],   t[j,:],  # normal and tangent
                                   p_g[i],   d_g,     # restitution coeffs
                                   p_Eb[i],  d_Eb[j], # Eb coeffs
                                   p_Gb[i],  d_Gb[j]) # Gb coeffs

        # normal elastic force
        f[i,:] += fne[:]

        # normal damping force
        f[i,:] += fnd[:]

        # tangential elastic force
        f[i,:] += fte[:]

        # tangential damping force
        f[i,:] += ftd[:]
