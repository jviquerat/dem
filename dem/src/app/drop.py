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
        self.plot_show = False
        self.plot_png  = True

        density = 2200    # steel
        young   = 210.0e9 # steel
        poisson = 0.25    # steel

        self.n_row = 25 # nb of particles on a row at start
        self.n_col = 35 # nb of particles on a col at start
        self.radius = 0.025

        self.p = particles(n           = self.n_row*self.n_col,
                           nt          = self.nt,
                           density     = density,
                           radius      = self.radius,
                           restitution = 0.9,
                           young       = young,
                           poisson     = poisson,
                           color       = "b",
                           store       = False)

        # Colors
        colors = np.array(['r', 'g', 'b', 'c', 'm', 'y', 'k'])
        self.p.c = colors[np.random.randint(0,len(colors),size=self.p.n)]

        self.d = domain_factory.create("rectangle",
                                       x_min      = 0.0,
                                       x_max      = 5.0,
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
        sep = 4.0*self.radius

        for i in range(self.n_row):
            for j in range(self.n_col):
                self.p.x[self.n_col*i+j,0] = sep + sep*i + self.radius*random.random()
                self.p.x[self.n_col*i+j,1] = 10*sep + sep*j

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

        if (self.p.store):
            plot_history(self.p.n, self.p.history, self.p.c)
