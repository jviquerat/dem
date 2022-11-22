# Generic imports
import math
import numpy as np
import matplotlib

import numba as nb
from   numba              import jit, njit
from   matplotlib.patches import Rectangle

# Custom imports
from dem.src.domain.base_domain import *

### ************************************************
### Distance from domain to given coordinates
@jit(cache=True, nopython=True, fastmath=True)
def domain_distance(a, b, c, d, x):

    dist    = np.zeros((4), np.float32)
    dist[:] = np.absolute(a[:]*x[0] +
                          b[:]*x[1] +
                          c[:])/d[:]

    return dist

### ************************************************
### Collision forces
@jit(cache=True, fastmath=True)
def forces(f, dx, r, alpha, p_sigma, p_kappa, d_sigma, d_kappa, m, v, n, t, ci, cj, n_coll):

    k_n  = np.zeros((n_coll))
    nu_n = np.zeros((n_coll))
    k_t  = np.zeros((n_coll))
    nu_t = np.zeros((n_coll))
    vn   = np.zeros((n_coll), np.float32)
    vt   = np.zeros((n_coll), np.float32)

    # normal stiffness
    k_n[:]  = (4.0/3.0)*np.sqrt(r[ci[:]])/(p_sigma[ci[:]] + d_sigma)

    # normal damping
    nu_n[:] = alpha[ci[:]]*np.sqrt(1.5*k_n*m[ci[:]])

    # tangential stiffness
    k_t  = 8.0*np.sqrt(r[ci[:]])/(p_kappa[ci[:]] + d_kappa)

    # tangential damping
    nu_t = alpha[ci[:]]*np.sqrt(k_t*m[ci[:]])

    vn[:]   = v[ci[:],0]*n[cj[:],0] + v[ci[:],1]*n[cj[:],1]
    vt[:]   = v[ci[:],0]*t[cj[:],0] + v[ci[:],1]*t[cj[:],1]

    # normal elastic force
    f[ci[:],0] += np.power(dx[:],1.5)*k_n[:]*n[cj[:],0]
    f[ci[:],1] += np.power(dx[:],1.5)*k_n[:]*n[cj[:],1]

    # normal damping force
    f[ci[:],0] -= np.power(dx[:],0.25)*nu_n[:]*vn[:]*n[cj[:],0]
    f[ci[:],1] -= np.power(dx[:],0.25)*nu_n[:]*vn[:]*n[cj[:],1]

    # tangential elastic force
    #p.f[i,:] -= pow(dx,0.5)*k_t*vt*0.00001*t[:]

    # tangential damping force
    f[ci[:],0] -= np.power(dx[:],0.25)*nu_t[:]*vt[:]*t[cj[:],0]
    f[ci[:],1] -= np.power(dx[:],0.25)*nu_t[:]*vt[:]*t[cj[:],1]


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
        self.y     = young   # young modulus
        self.p     = poisson # poisson ratio
        self.sigma = (1.0 - self.p**2)/self.y
        self.kappa = 2.0*(2.0 + self.p)*(1.0 - self.p)/self.y

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
    ### Compute collisions with a particle
    def collisions(self, p):

        ci = np.empty((0), dtype=int)
        cj = np.empty((0), dtype=int)
        cd = np.empty((0), dtype=float)

        for i in range(p.n):
            dist = domain_distance(self.a, self.b, self.c, self.d, p.x[i,:])

            for j in range(4):
                dx = dist[j] - p.r[i] # relative distance

                if (dx < 0.0):
                    ci = np.append(ci, i)
                    cj = np.append(cj, j)
                    cd = np.append(cd, dx)

        n_coll = len(ci)
        if (n_coll == 0): return

        cd = np.abs(cd)

        forces(p.f, cd, p.r, p.alpha, p.sigma, p.kappa, self.sigma, self.kappa, p.m, p.v, self.n, self.t, ci, cj, n_coll)
