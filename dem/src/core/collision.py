# Generic imports
import math
import numpy as np
import numba as nb

### ************************************************
### Compute forces given input parameters of collision
@nb.njit(cache=True)
def hertz(dx, r1, r2, m1, m2, v1, v2, n, t,
          alpha1, alpha2, sigma1, sigma2, kappa1, kappa2):

    # averaged values
    r = 1.0/(1.0/r1 + 1.0/r2)
    m = 1.0/(1.0/m1 + 1.0/m2)

    # normal stiffness
    k_n  = (4.0/3.0)*np.sqrt(r)/(sigma1 + sigma2)

    # normal damping
    nu_n = alpha1*np.sqrt(1.5*k_n*m)

    # tangential stiffness
    k_t  = 8.0*np.sqrt(r)/(kappa1 + kappa2)

    # tangential damping
    nu_t = alpha1*np.sqrt(k_t*m)

    vn   = v1[0]*n[0] + v1[1]*n[1]
    vt   = v1[0]*t[0] + v1[1]*t[1]

    # forces
    fne = np.zeros((2))
    fnd = np.zeros((2))
    fte = np.zeros((2))
    ftd = np.zeros((2))

    dx15  = pow(dx, 1.5)
    dx025 = pow(dx, 0.25)

    # normal elastic force
    fne[0] = dx15*k_n*n[0]
    fne[1] = dx15*k_n*n[1]

    # normal damping force
    fnd[0] =-dx025*nu_n*vn*n[0]
    fnd[1] =-dx025*nu_n*vn*n[1]

    # tangential elastic force
    #p.f[i,:] -= pow(dx,0.5)*k_t*vt*0.00001*t[:]

    # tangential damping force
    ftd[0] =-dx025*nu_t*vt*t[0]
    ftd[1] =-dx025*nu_t*vt*t[1]

    return fne, fnd, fte, ftd
