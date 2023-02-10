# Generic imports
import os
import math
import random

# Custom imports
from dem.src.app.base_app import *
from dem.src.plot.plot    import *

### ************************************************
### Dropping of several spheres to check inter-particle contacts
class drop(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self):
        super().__init__()

        self.name      = 'drop'
        self.t_max     = 10.0
        self.dt        = 0.00002
        self.nt        = int(self.t_max/self.dt)+1
        self.plot_freq = 1000
        self.plot_it   = 0
        self.plot_show = True
        self.plot_png  = False

        density = 2200    # steel
        young   = 210.0e9 # steel
        poisson = 0.25    # steel

        self.p = particles(n           = 10,
                           nt          = self.nt,
                           density     = density,
                           radius      = 0.05,
                           restitution = 1.0,
                           young       = young,
                           poisson     = poisson,
                           color       = "b",
                           store       = True)

        # Restitution ratios
        self.p.e[:] = 0.8
        self.p.set_particles()

        # Colors
        self.p.c[0] = 'b'
        self.p.c[1] = 'r'
        self.p.c[2] = 'y'
        self.p.c[3] = 'b'
        self.p.c[4] = 'r'
        self.p.c[5] = 'y'
        self.p.c[6] = 'b'
        self.p.c[7] = 'r'
        self.p.c[8] = 'y'
        self.p.c[9] = 'b'

        self.d = domain_factory.create("rectangle",
                                       x_min      = 0.0,
                                       x_max      = 2.0,
                                       y_min      = 0.0,
                                       y_max      = 5.0,
                                       young      = young,
                                       poisson    = poisson)

        self.path = self.name
        os.makedirs(self.path, exist_ok=True)

        self.reset()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.t = 0.0

        for i in range(2):
            for j in range(5):
                self.p.x[5*i+j,0] = 0.8 + 0.2*i
                self.p.x[5*i+j,1] = 1.0 + 0.5*j
                self.p.v[5*i+j,0] = random.random()

    ### ************************************************
    ### Compute forces
    def forces(self):

        self.p.reset_forces()
        self.p.collisions()
        self.d.collisions(self.p)
        self.p.gravity(self.g)

    ### ************************************************
    ### Update positions
    def update(self, it):

        self.p.update(self.dt, it)
        self.t += self.dt

    ### ************************************************
    ### Plot
    def plot(self, it):

        if (it%self.plot_freq == 0):
            plot(self.d, self.p, self.path, self.plot_it,
                 show=self.plot_show, png=self.plot_png)
            self.plot_it += 1

    ### ************************************************
    ### Finalize
    def finalize(self):

        plot_history(self.p.n, self.p.history, self.p.c)
