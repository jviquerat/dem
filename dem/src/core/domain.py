# Generic imports
import math
import numpy as np

### ************************************************
### Class defining domain
class domain:
    ### ************************************************
    ### Constructor
    def __init__(self,
                 dtype   = "rectangle",
                 x_min   = 0.0,
                 x_max   = 1.0,
                 y_min   = 0.0,
                 y_max   = 1.0,
                 young   = 210.0e9,
                 poisson = 0.25):

        # Type of domain
        self.dtype = dtype

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
        # Ridge 0 is the bottom one, then we pursue
        # in trigonometric order
        if (dtype == "rectangle"):
            self.a = np.zeros((4))
            self.b = np.zeros((4))
            self.c = np.zeros((4))
            self.n = np.zeros((4,2))

            # Bottom
            self.a[0]   = 0.0
            self.b[0]   = 1.0
            self.c[0]   =-self.y_min
            self.n[0,0] = 0.0
            self.n[0,1] = 1.0

            # Right
            self.a[1]   = 1.0
            self.b[1]   = 0.0
            self.c[1]   =-self.x_max
            self.n[1,0] =-1.0
            self.n[1,1] = 0.0

            # Top
            self.a[2]   = 0.0
            self.b[2]   = 1.0
            self.c[2]   =-self.y_max
            self.n[2,0] = 0.0
            self.n[2,1] =-1.0

            # Left
            self.a[3]   = 1.0
            self.b[3]   = 0.0
            self.c[3]   =-self.x_min
            self.n[3,0] = 1.0
            self.n[3,1] = 0.0

    ### ************************************************
    ### Compute collisions with a particle
    def collisions(self, p):

        d_min = 1.0e8

        if (self.dtype == "rectangle"):
            n  = np.zeros((2)) # normal  to border
            t  = np.zeros((2)) # tangent to border
            x  = np.zeros((2)) # position
            d  = np.zeros((2)) # displacement
            v  = np.zeros((2)) # velocity
            vn = np.zeros((2)) # normal     velocity
            vt = np.zeros((2)) # tangential velocity
            for i in range(p.n):
                m    = p.m[i]     # mass
                r    = p.r[i]     # radius
                s    = p.sigma[i] # sigma
                k    = p.kappa[i] # kappa
                x[:] = p.x[i,:]   # position
                v[:] = p.v[i,:]   # velocity
                d[:] = p.d[i,:]   # displacement

                # normal stiffness
                k_n  = (4.0/3.0)*math.sqrt(r)/(s + self.sigma)

                # normal damping
                nu_n = p.alpha[i]*math.sqrt(1.5*k_n*m)

                # tangential stiffness
                k_t  = 8.0*math.sqrt(r)/(k + self.kappa)

                # tangential damping
                nu_t = p.alpha[i]*math.sqrt(k_t*m)

                for j in range(4):
                    a    = self.a[j]   # coefficients
                    b    = self.b[j]   # coefficients
                    c    = self.c[j]   # coefficients
                    n[:] = self.n[j,:] # normal  to boundary
                    t[0] = n[1]        # tangent to boundary
                    t[1] =-n[0]        # tangent to boundary

                    dist  = abs(a*x[0] + b*x[1] + c) # distance to border
                    dist /= math.sqrt(a*a + b*b)     # distance to border
                    dx = dist - r                    # relative distance

                    if (dx < 0.0):
                        vn   = np.dot(v,n) # normal     velocity
                        vt   = np.dot(v,t) # tangential velocity
                        dn   = np.dot(d,n)
                        dt   = np.dot(d,t)
                        dx   = abs(dx)

                        # normal elastic force
                        p.a[i,:] += pow(dx,1.5)*k_n*n[:]
                        #print(pow(dx,1.5)*k_n*n[:])

                        # normal damping force
                        p.a[i,:] -= pow(dx,0.25)*nu_n*vn*n[:]
                        #print(-pow(dx,0.25)*nu_n*vn*n[:])

                        # tangential elastic force
                        #p.a[i,:] += pow(dx,0.5)*k_t*dt*t[:]
                        #print(pow(dx,0.5)*k_t*dt*t[:])

                        # tangential damping force
                        p.a[i,:] -= pow(dx,0.25)*nu_t*vt*t[:]
                        #print(-pow(dx,0.25)*nu_t*vt*t[:])

                        #print("")

        if (self.dtype == "circle"):
            pass
