# Generic imports
import os
import math

# Custom imports
from dem.src.app.base_app import *

### ************************************************
### Single sphere under gravity with perfect restitution
class gravity(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self):
        super().__init__()

        self.name      = 'gravity'
        self.t_max     = 3.0
        self.dt        = 0.000025
        self.nt        = int(self.t_max/self.dt)+1
        self.plot_freq = 1000
        self.plot_it   = 0
        self.plot_show = True
        self.plot_png  = False

        density = 2200    # steel
        young   = 210.0e9 # steel
        poisson = 0.25    # steel

        self.p = particles(n           = 1,
                           nt          = self.nt,
                           density     = density,
                           radius      = 0.05,
                           restitution = 1.0,
                           young       = young,
                           poisson     = poisson,
                           color       = "b",
                           store       = True)

        self.d = domain_factory.create("rectangle",
                                       x_min      = 0.0,
                                       x_max      = 1.0,
                                       y_min      = 0.0,
                                       y_max      = 0.5,
                                       young      = young,
                                       poisson    = poisson)

        self.path = self.name
        os.makedirs(self.path, exist_ok=True)

        self.reset()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.t = 0.0

        self.p.x[:,0] = 0.5
        self.p.x[:,1] = 0.4

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
