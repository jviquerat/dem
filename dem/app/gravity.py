# Generic imports
import os
import math

# Custom imports
from dem.app.base_app import *

### ************************************************
### Single sphere under gravity with perfect restitution
class gravity(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 name            = 'gravity',
                 t_max           = 1.0,
                 dt              = 2.5e-5,
                 plot_freq       = 500,
                 plot_show       = True,
                 plot_trajectory = True,
                 plot_png        = True):
        super().__init__()

        self.name            = name
        self.t_max           = t_max
        self.dt              = dt
        self.plot_freq       = plot_freq
        self.plot_show       = plot_show
        self.plot_trajectory = plot_trajectory
        self.plot_png        = plot_png

        self.nt        = int(self.t_max/self.dt)
        self.plot_it   = 0

        self.p = particles(np          = 1,
                           nt          = self.nt,
                           material    = "steel",
                           radius      = 0.05,
                           color       = "b",
                           store       = True)

        self.d = domain_factory.create("rectangle",
                                       x_min      = 0.0,
                                       x_max      = 0.3,
                                       y_min      = 0.0,
                                       y_max      = 0.5,
                                       angle      = 0.0,
                                       material   = "steel")

        # Set perfect restitution
        self.p.e_wall[:] = 1.0

        self.path = self.base_path+'/'+self.name
        os.makedirs(self.path, exist_ok=True)

        self.reset()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.it = 0
        self.t  = 0.0

        self.p.x[:,0] = 0.15
        self.p.x[:,1] = 0.4

    ### ************************************************
    ### Compute forces
    def forces(self):

        self.p.reset_forces()
        self.d.collisions(self.p, self.dt)
        self.p.gravity(self.g)
