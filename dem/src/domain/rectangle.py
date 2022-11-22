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
@nb.njit(cache=True)
def domain_distance(a, b, c, d, x, r, n):

    ci = np.empty((0), np.uint16)
    cj = np.empty((0), np.uint16)
    cd = np.empty((0), np.float32)

    for i in range(n):
        for j in range(4):
            dist = np.absolute(a[j]*x[i,0] + b[j]*x[i,1] + c[j])/d[j]
            dx   = dist - r[i] # relative distance

            if (dx < 0.0):
                ci = np.append(ci, np.uint16(i))
                cj = np.append(cj, np.uint16(j))
                cd = np.append(cd, np.float32(abs(dx)))

    return ci, cj, cd

### ************************************************
### Collision forces
@njit(cache=True)
def forces(f, dx, r, alpha, p_sigma, p_kappa, d_sigma, d_kappa, m, v, n, t, ci, cj, n_coll):

    for k in range(n_coll):
        i = ci[k]
        j = cj[k]

        # normal stiffness
        k_n  = (4.0/3.0)*np.sqrt(r[i])/(p_sigma[i] + d_sigma[i])

        # normal damping
        nu_n = alpha[i]*np.sqrt(1.5*k_n*m[i])

        # tangential stiffness
        k_t  = 8.0*np.sqrt(r[i])/(p_kappa[i] + d_kappa[i])

        # tangential damping
        nu_t = alpha[i]*np.sqrt(k_t*m[i])

        vn   = v[i,0]*n[j,0] + v[i,1]*n[j,1]
        vt   = v[i,0]*t[j,0] + v[i,1]*t[j,1]

        # normal elastic force
        f[i,0] += np.power(dx[k],1.5)*k_n*n[j,0]
        f[i,1] += np.power(dx[k],1.5)*k_n*n[j,1]

        # normal damping force
        f[i,0] -= np.power(dx[k],0.25)*nu_n*vn*n[j,0]
        f[i,1] -= np.power(dx[k],0.25)*nu_n*vn*n[j,1]

        # tangential elastic force
        #p.f[i,:] -= pow(dx,0.5)*k_t*vt*0.00001*t[:]

        # tangential damping force
        f[i,0] -= np.power(dx[k],0.25)*nu_t*vt*t[j,0]
        f[i,1] -= np.power(dx[k],0.25)*nu_t*vt*t[j,1]


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

        # Compute distances to domain boundaries
        ci, cj, cd = domain_distance(self.a, self.b, self.c, self.d, p.x, p.r, p.n)

        # Check if there are collisions
        n_coll = len(ci)
        if (n_coll == 0): return

        # Generate arrays for force computation
        sigma = np.ones((n_coll))*self.sigma
        kappa = np.ones((n_coll))*self.kappa

        # Compute forces
        forces(p.f, cd, p.r, p.alpha, p.sigma, p.kappa, sigma, kappa, p.m, p.v, self.n, self.t, ci, cj, n_coll)
