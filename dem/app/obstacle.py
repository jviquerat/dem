# Generic imports
import os
import math

# Custom imports
from dem.app.base_app import *

### ************************************************
### Single sphere under gravity with obstacle in the domain
class obstacle(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 name            = 'obstacle',
                 t_max           = 2.5,
                 dt              = 2.5e-5,
                 angle           = 0.0,
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
                           radius      = 0.02,
                           color       = "b",
                           store       = True)

        self.p.e_wall[:] = 0.8

        self.d = domain_factory.create("rectangle",
                                       x_min      = 0.0,
                                       x_max      = 0.3,
                                       y_min      = 0.0,
                                       y_max      = 0.5,
                                       angle      = 0.0,
                                       material   = "steel")

        self.o1 = domain_factory.create("rectangle",
                                        x_min      =-0.05,
                                        x_max      = 0.15,
                                        y_min      = 0.2,
                                        y_max      = 0.22,
                                        angle      =-20.0,
                                        plot_fill  = True,
                                        material   = "steel")
        self.o2 = domain_factory.create("rectangle",
                                        x_min      = 0.15,
                                        x_max      = 0.35,
                                        y_min      = 0.1,
                                        y_max      = 0.12,
                                        angle      = 20.0,
                                        plot_fill  = True,
                                        material   = "steel")

        self.d_lst = [self.d, self.o1, self.o2]

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
        self.o1.collisions(self.p, self.dt)
        self.o2.collisions(self.p, self.dt)
        self.p.gravity(self.g)
