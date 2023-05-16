# Generic imports
import math
import numpy as np
import numba as nb

### ************************************************
### Compute forces given input parameters of collision
### Penetration dx is an absolute distance
### Normals must be from 1 towards 2
### Relative velocity is computed as v12 = v2-v1
#@nb.njit(cache=True)
def hertz(dx, dt, r1, r2, m1, m2, v1, v2, n,
          e1, e2, Y1, Y2, G1, G2, mu1, mu2):

    # averaged values
    R = 1.0/(1.0/r1 + 1.0/r2)
    M = 1.0/(1.0/m1 + 1.0/m2)
    Y = 1.0/(Y1 + Y2)
    G = 1.0/(G1 + G2)
    E = 0.5*(e1 + e2) + 1.0e-8
    L = np.log(E)/np.sqrt(math.pi**2+np.log(E)**2)
    C = 0.5*(mu1 + mu2)

    # relative velocity
    v    = np.zeros((2), np.float32)
    v[:] = v2[:] - v1[:]

    # normal stiffness
    k_n = (4.0/3.0)*Y*np.sqrt(R)

    # normal damping
    g_n =-L*np.sqrt(5.0*k_n*M)

    # tangential stiffness
    k_t = 8.0*G*np.sqrt(R)

    # tangential damping
    g_t =-L*np.sqrt(5.0*k_t*M/6.0)

    # normal and tangential velocities
    t     = np.zeros((2), np.float32)
    vn    = np.zeros((2), np.float32)
    vt    = np.zeros((2), np.float32)
    vn[:] = np.dot(v, n)*n[:]
    vt[:] = v[:] - vn[:]
    t[:]  = vt[:]/np.sqrt(np.dot(vt[:],vt[:]) + 1.0e-8)

    # tangential displacement
    # very simple approximation
    dxt    = np.zeros((2))
    dxt[:] = vt[:]*dt

    # forces
    fn = np.zeros((2), np.float32)
    ft = np.zeros((2), np.float32)

    dx15  = pow(dx, 1.5)
    dx05  = pow(dx, 0.5)
    dx025 = pow(dx, 0.25)

    # normal force
    fn[:] = dx15*k_n*n[:]   - dx025*g_n*vn[:]

    # tangential force
    ft[:] =-dx05*k_t*dxt[:] - dx025*g_t*vt[:]

    # apply coulomb law
    vfn = np.sqrt(np.dot(fn, fn))
    vft = np.sqrt(np.dot(ft, ft))
    if (vft > C*vfn):
        ft[:] =-C*vfn*t[:]

    return fn, ft
