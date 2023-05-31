# Generic imports
import os
import math
import random

# Custom imports
from dem.app.base_app import *

### ************************************************
### Dropping of several spheres to check inter-particle contacts
class drop(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 name            = 'drop',
                 t_max           = 2.0,
                 dt              = 2.5e-5,
                 plot_freq       = 1000,
                 plot_show       = False,
                 plot_trajectory = False,
                 plot_png        = True):
        super().__init__()

        self.name            = name
        self.t_max           = t_max
        self.dt              = dt
        self.plot_freq       = plot_freq
        self.plot_show       = plot_show
        self.plot_trajectory = plot_trajectory
        self.plot_png        = plot_png

        self.nt      = int(self.t_max/self.dt)
        self.plot_it = 0

        self.n_row  = 40 # nb of particles on a row at start
        self.n_col  = 40 # nb of particles on a col at start
        self.radius = 0.025

        self.p = particles(np       = self.n_row*self.n_col,
                           nt       = self.nt,
                           material = "steel",
                           radius   = self.radius,
                           color    = "b",
                           store    = False,
                           search   = "nearest")

        # Set different resitution ratio
        self.p.e_wall[:] = 0.7

        # Colors
        colors = np.array(['r', 'g', 'b', 'c', 'm', 'y', 'k'])
        self.p.c = colors[np.random.randint(0,len(colors),size=self.p.np)]

        self.d = domain_factory.create("rectangle",
                                       x_min      = 0.0,
                                       x_max      = 5.0,
                                       y_min      = 0.0,
                                       y_max      = 5.0,
                                       material   = "steel")

        self.path = self.base_path+'/'+self.name
        os.makedirs(self.path, exist_ok=True)

        self.reset()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.it = 0
        self.t  = 0.0
        sep = 4.0*self.radius

        for i in range(self.n_row):
            for j in range(self.n_col):
                self.p.x[self.n_col*i+j,0] = sep + sep*i + self.radius*random.random()
                self.p.x[self.n_col*i+j,1] = 10*sep + sep*j

    ### ************************************************
    ### Compute forces
    def forces(self):

        self.p.reset_forces()
        self.p.collisions(self.dt)
        self.d.collisions(self.p, self.dt)
        self.p.gravity(self.g)
