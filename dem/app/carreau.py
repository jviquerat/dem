# Generic imports
import os
import math

# Custom imports
from dem.app.base_app import *

### ************************************************
### Coup type carreau entre deux billes
class carreau(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 name            = 'carreau',
                 t_max           = 1.0,
                 dt              = 2.5e-5,
                 plot_freq       = 200,
                 plot_show       = True,
                 plot_trajectory = False,
                 plot_png        = False):
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

        self.p = particles(np          = 2,
                           nt          = self.nt,
                           material    = "steel",
                           radius      = 0.05,
                           color       = "b",
                           store       = True,
                           search      = "nearest")

        self.s1 = material_factory.create("steel")
        self.s1.e_part = 1.0
        self.s1.e_wall = 1.0
        self.p.mat[0] = self.s1
        self.p.mat[1] = self.s1

        # Colors
        self.p.c[0] = 'b'
        self.p.c[1] = 'r'

        self.d = domain_factory.create("rectangle",
                                       x_min      = 0.0,
                                       x_max      = 0.5,
                                       y_min      = 0.0,
                                       y_max      = 0.5,
                                       material   = "steel")

        self.path = self.base_path+'/'+self.name
        os.makedirs(self.path, exist_ok=True)

        self.reset()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.it = 0
        self.t  = 0.0

        self.p.x[0,0] = 0.1
        self.p.x[1,0] = 0.25
        self.p.x[:,1] = 0.25
        self.p.v[0,0] = 2.0

    ### ************************************************
    ### Compute forces
    def forces(self):

        self.p.reset_forces()
        self.p.collisions(self.dt)
        self.d.collisions(self.p, self.dt)
