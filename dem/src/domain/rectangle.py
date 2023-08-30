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
### Class defining rectangle domain
class rectangle(base_domain):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 x_min         = 0.0,
                 x_max         = 1.0,
                 y_min         = 0.0,
                 y_max         = 1.0,
                 angle         = 0.0, # in degrees
                 plot_fill     = False,
                 material      = "steel",
                 open_boundary = [False, False, False, False]):

        # Type
        self.type = "rectangle"

        # Angle
        self.angle = angle

        # Lateral sizes
        self.dx = x_max - x_min
        self.dy = y_max - y_min

        # Plot filling
        self.plot_fill = plot_fill

        # Define center of domain
        self.pc = np.array([0.5*(x_max+x_min), 0.5*(y_max+y_min)])

        # Material
        self.mat = material_factory.create(material)

        # Open boundaries
        self.open_boundary = np.array(open_boundary)

        # Ficticious parameters for domain
        self.r     = 1.0e8
        self.m     = 1.0e8
        self.v     = np.zeros((2))

        # Define boundary points of the 4 segments
        # First segment is bottom one, then right one,
        # then top one, then left one
        self.p1 = np.array([x_min, y_min]) # bottom left
        self.p2 = np.array([x_max, y_min]) # bottom right
        self.p3 = np.array([x_max, y_max]) # top right
        self.p4 = np.array([x_min, y_max]) # top left

        # Rotate points around center
        self.rotate(self.p1, self.pc, self.angle)
        self.rotate(self.p2, self.pc, self.angle)
        self.rotate(self.p3, self.pc, self.angle)
        self.rotate(self.p4, self.pc, self.angle)

        # External boundaries
        self.x_min = min(self.p1[0], self.p2[0], self.p3[0], self.p4[0])
        self.x_max = max(self.p1[0], self.p2[0], self.p3[0], self.p4[0])
        self.y_min = min(self.p1[1], self.p2[1], self.p3[1], self.p4[1])
        self.y_max = max(self.p1[1], self.p2[1], self.p3[1], self.p4[1])

        # Define segments
        self.seg_pts = np.array([ [self.p1, self.p2],
                                  [self.p2, self.p3],
                                  [self.p3, self.p4],
                                  [self.p4, self.p1] ])

    ### ************************************************
    ### Compute collisions with a set of particles
    def collisions(self, p, dt):

        # Search for collisions linearly and compute forces
        linear_search(self.seg_pts, self.r, self.m,
                      self.v, self.mat.Y, self.mat.G, self.open_boundary,
                      p.x, p.r, p.m, p.v, p.e_wall, p.mu_wall,
                      p.Y, p.G, p.f, p.np, dt)

    ### ************************************************
    ### Rotate point around another point from given angle
    def rotate(self, p, pc, angle):

        cost  = math.cos(math.radians(angle))
        sint  = math.sin(math.radians(angle))
        dp    = np.zeros((2))
        dp[0] = (p[0]-pc[0])*cost - (p[1]-pc[1])*sint + pc[0]
        dp[1] = (p[0]-pc[0])*sint + (p[1]-pc[1])*cost + pc[1]
        p[:]  = dp[:]

    ### ************************************************
    ### Check if point is in domain
    ### pt is assumed to be an np array of size 2
    def is_in(self, pm):

        p1p2 = self.p2 - self.p1
        p1pm =      pm - self.p1
        p1p4 = self.p4 - self.p1

        p1p2p1pm = np.dot(p1p2, p1pm)
        p1p2p1p2 = np.dot(p1p2, p1p2)
        p1p4p1pm = np.dot(p1p4, p1pm)
        p1p4p1p4 = np.dot(p1p4, p1p4)

        if ((p1p2p1pm > 0.0)      and
            (p1p2p1p2 > p1p2p1pm) and
            (p1p4p1pm > 0.0)      and
            (p1p4p1p4 > p1p4p1pm)): return True

        return False

### ************************************************
### Distance from rectangle domain to given coordinates
### Prefix d_ corresponds to domain
### Prefix p_ corresponds to particle
@nb.njit(cache=True)
def linear_search(d_pts, d_r, d_m,
                  d_v, d_mat_Y, d_mat_G, open_boundary,
                  p_x, p_r, p_m, p_v, p_e_wall, p_mu_wall,
                  p_Y, p_G, p_f, n, dt):

    # Loop on particles
    for i in range(n):

        # Loop on rectangle sides
        for j in range(4):

            # If this boundary is open, loop
            if (open_boundary[j]): continue

            dx, nrm = p_to_segment(p_x[i], d_pts[j,0], d_pts[j,1])
            dx      = dx - p_r[i]

            # If particle intersects boundary
            if (dx < 0.0):

                # Return forces from collision parameters
                # - normal elastic
                # - normal damping,
                # - tangential elastic
                # - tangential damping
                fn, ft = hertz(abs(dx),      # penetration
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

### ************************************************
### Compute distance from point p to segment [p1,p2]
@nb.njit(cache=True)
def p_to_segment(p, p1, p2):

    x,  y  =  p[0],  p[1]
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]

    a = x  - x1
    b = y  - y1
    c = x2 - x1
    d = y2 - y1

    dot   = a*c + b*d
    lgt   = c*c + d*d
    ratio = dot/lgt

    if (ratio < 0.0):
        px = x1;
        py = y1;
    elif (ratio > 1.0):
        px = x2;
        py = y2;
    else:
        px = x1 + ratio*c
        py = y1 + ratio*d

    dx   = px - x
    dy   = py - y
    dist = math.sqrt(dx*dx + dy*dy)
    nrm  = np.array([dx, dy])/dist # outward normal

    return dist, nrm
