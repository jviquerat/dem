# Generic imports
import math
import numpy as np
import matplotlib

from   matplotlib.patches import Rectangle

# Custom imports
from dem.src.domain.base_domain import *

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
        # normalization parameter sqrt(a*a + b*b)
        self.a = np.array([0.0, 1.0, 0.0, 1.0])
        self.b = np.array([1.0, 0.0, 1.0, 0.0])
        self.c = np.array([-self.y_min,-self.x_max,
                           -self.y_max,-self.x_min])
        self.d = np.array([1.0, 1.0, 1.0, 1.0])
        self.n = np.array([[ 0.0, 1.0],
                           [-1.0, 0.0],
                           [ 0.0,-1.0],
                           [ 1.0, 0.0]])

    ### ************************************************
    ### Distance to given coordinates
    def distance(self, x):

        dist    = np.zeros((4))
        dist[:] = np.absolute(self.a[:]*x[0] +
                              self.b[:]*x[1] +
                              self.c[:])/self.d[:]

        return dist

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

        d_min = 1.0e8

        n  = np.zeros((2)) # normal  to border
        t  = np.zeros((2)) # tangent to border
        x  = np.zeros((2)) # position
        d  = np.zeros((2)) # displacement
        v  = np.zeros((2)) # velocity
        vn = np.zeros((2)) # normal     velocity
        vt = np.zeros((2)) # tangential velocity

        coll = []
        for i in range(p.n):

            x[:] = p.x[i,:] # position
            r    = p.r[i]   # radius
            dist = self.distance(x)

            for j in range(4):
                dx = dist[j] - r   # relative distance

                if (dx < 0.0):
                    coll.append([i,j,dx])

        # n_coll = len(coll)
        # k_n = np.zeros((n_coll))
        # nu_n = np.zeros((n_coll))
        # k_t = np.zeros((n_coll))
        # nu_t = np.zeros((n_coll))

        for k in range(len(coll)):
            i  = coll[k][0]
            j  = coll[k][1]
            dx = coll[k][2]

            # normal stiffness
            k_n  = (4.0/3.0)*math.sqrt(p.r[i])/(p.sigma[i] + self.sigma)

            # normal damping
            nu_n = p.alpha[i]*math.sqrt(1.5*k_n*p.m[i])

            # tangential stiffness
            k_t  = 8.0*math.sqrt(r)/(p.kappa[i] + self.kappa)

            # tangential damping
            nu_t = p.alpha[i]*math.sqrt(k_t*p.m[i])

            n[:] = self.n[j,:] # normal  to boundary
            t[0] = n[1]        # tangent to boundary
            t[1] =-n[0]        # tangent to boundary
            vn   = np.dot(p.v[i],n) # normal     velocity
            vt   = np.dot(p.v[i],t) # tangential velocity
            dn   = np.dot(p.d[i],n)
            dt   = np.dot(p.d[i],t)
            dx   = abs(dx)

            # normal elastic force
            p.f[i,:] += pow(dx,1.5)*k_n*n[:]

            # normal damping force
            p.f[i,:] -= pow(dx,0.25)*nu_n*vn*n[:]

            # tangential elastic force
            #p.f[i,:] -= pow(dx,0.5)*k_t*vt*0.00001*t[:]

            # tangential damping force
            p.f[i,:] -= pow(dx,0.25)*nu_t*vt*t[:]
