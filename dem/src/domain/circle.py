# Generic imports
import math
import matplotlib
import numpy as np
import numba as nb

from matplotlib.patches import Rectangle

# Custom imports
from dem.src.domain.base       import *
from dem.src.core.collision    import *
from dem.src.material.material import *

### ************************************************
### Class defining circular domain
class circle(base_domain):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 x_c       = 0.0,
                 y_c       = 0.0,
                 rad       = 1.0,
                 plot_fill = False,
                 velocity  = 0.0,
                 material  = "steel"):

        # Type
        self.type = "circle"

        # Center coordinates and radius
        self.x_c = x_c
        self.y_c = y_c
        self.rad = rad

        # Rolling velocity
        self.velocity = velocity

        # Plot filling
        self.plot_fill = plot_fill

        # Material
        self.mat = material_factory.create(material)

        # Ficticious parameters for domain
        self.r = 1.0e8
        self.m = 1.0e8
        self.v = np.array([velocity, 0.0])#np.zeros((2))

        # External boundaries
        self.x_min = x_c - rad
        self.x_max = x_c + rad
        self.y_min = y_c - rad
        self.y_max = y_c + rad

    ### ************************************************
    ### Compute collisions with a set of particles
    def collisions(self, p, dt):

        # Search for collisions linearly and compute forces
        linear_search(self.x_c, self.y_c, self.rad, self.m,
                      self.v, self.mat.Y, self.mat.G,
                      p.x, p.r, p.m, p.v, p.e_wall, p.mu_wall,
                      p.Y, p.G, p.f, p.np, dt)

### ************************************************
### Distance from rectangle domain to given coordinates
### Prefix d_ corresponds to domain
### Prefix p_ corresponds to particle
@nb.njit(cache=True)
def linear_search(d_xc, d_yc, d_r, d_m,
                  d_v, d_mat_Y, d_mat_G,
                  p_x, p_r, p_m, p_v, p_e_wall, p_mu_wall,
                  p_Y, p_G, p_f, n, dt):

    # Loop on particles
    for i in range(n):

        # Check distance to circle border
        dx, nrm = p_to_p(p_x[i], [d_xc, d_yc])
        if ((dx > d_r + p_r[i]) or (dx < d_r - p_r[i])): continue

        if (dx < d_r):
            dist   = dx + p_r[i] - d_r
            nrm[:] = -nrm[:]
            d_rad  = -d_r
        if (dx > d_r):
            dist  = dx - p_r[i] - d_r
            d_rad = d_r

        # Compute approximate velocity at impact point
        cost = p_x[i,0]/d_r
        sint = p_x[i,1]/d_r
        nv   = np.linalg.norm(d_v)
        v    = nv*np.array([-sint, cost])

        # Return forces from collision parameters
        # - normal elastic
        # - normal damping,
        # - tangential elastic
        # - tangential damping
        fn, ft = hertz(abs(dist),    # penetration
                       dt,           # timestep
                       p_r[i],       # radius 1
                       d_rad,        # radius 2
                       p_m[i],       # mass 1
                       d_m,          # mass 2
                       p_v[i,:],     # velocity 1
                       v[:],       # velocity 2
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

### ************************************************
### Compute distance from point p to point pc
@nb.njit(cache=True)
def p_to_p(p, pc):

    x,  y  =  p[0],  p[1]
    xc, yc = pc[0], pc[1]

    dx   = xc - x
    dy   = yc - y
    dist = math.sqrt(dx*dx + dy*dy)
    nrm  = np.array([dx, dy])/dist # outward normal

    return dist, nrm

### ************************************************
### Compute domain angular velocity at contact point
@nb.njit(cache=True)
def v_at_p(p, pc):

    x,  y  =  p[0],  p[1]
    xc, yc = pc[0], pc[1]

    dx   = xc - x
    dy   = yc - y
    dist = math.sqrt(dx*dx + dy*dy)
    nrm  = np.array([dx, dy])/dist # outward normal

    return dist, nrm
