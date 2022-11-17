# Generic imports
import os
import math

# Custom imports
from dem.src.app.base_app import *
from dem.src.plot.plot    import *

### ************************************************
### Gravity fall
class gravity(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self):
        super().__init__()

        self.name      = 'gravity'
        self.t_max     = 1.0
        self.dt        = 0.000025
        self.plot_freq = 1000
        self.plot_it   = 0
        self.plot_show = True
        self.plot_png  = False

        density = 2200    # steel
        young   = 210.0e9 # steel
        poisson = 0.25    # steel

        self.p = particles(n           = 3,
                           density     = density,
                           radius      = 0.05,
                           restitution = 1.0,
                           young       = young,
                           poisson     = poisson,
                           color       = "b",
                           store       = True)

        # Restitution ratios
        self.p.e[0] = 1.0
        self.p.e[1] = 0.8
        self.p.e[2] = 0.5
        self.p.set_particles()

        # Colors
        self.p.c[0] = 'b'
        self.p.c[1] = 'r'
        self.p.c[2] = 'y'

        self.d = domain(dtype      = "rectangle",
                        x_min      = 0.0,
                        x_max      = 1.0,
                        y_min      = 0.0,
                        y_max      = 0.5,
                        young      = young,
                        poisson    = poisson)

        self.path = self.name
        os.makedirs(self.path, exist_ok=True)

        self.reset()
        self.r_min = self.p.min_radius()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.t = 0.0

        self.p.x[:,0] = 0.5
        self.p.x[:,1] = 0.4
        self.p.v[:,0] = 1.0

    ### ************************************************
    ### Compute forces
    def forces(self):

        self.p.reset_forces()
        self.p.collisions()
        self.d.collisions(self.p)
        self.p.gravity()

    ### ************************************************
    ### Update positions
    def update(self):

        self.p.update(self.dt)
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
