# Generic imports
import os
import math
import random

# Custom imports
from dem.app.base_app import *

### ************************************************
### Rolling mill test case
class mill(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 name            = 'mill',
                 t_max           = 4.0,
                 dt              = 2.5e-5,
                 plot_freq       = 500,
                 plot_show       = True,
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

        self.nt        = int(self.t_max/self.dt)
        self.plot_it   = 0

        self.n_row  = 20 # nb of particles on a row at start
        self.n_col  = 20 # nb of particles on a col at start
        self.radius = 0.025

        self.p = particles(np          = self.n_row*self.n_col,
                           nt          = self.nt,
                           material    = "glass",
                           radius      = self.radius,
                           color       = "b",
                           store       = False,
                           search      = "nearest",
                           rad_coeff   = 2.0)

        self.p.e_wall[:] = 0.5
        self.p.e_part[:] = 0.5

        colors = np.array(['r', 'g', 'b', 'c', 'm', 'y', 'k'])
        self.p.c = colors[np.random.randint(0,len(colors),size=self.p.np)]

        self.d = domain_factory.create("circle",
                                       x_c      = 0.0,
                                       y_c      = 0.0,
                                       rad      = 1.0,
                                       velocity = 2.0,
                                       material = "steel")

        self.d_lst = [self.d]

        self.path = self.base_path+'/'+self.name
        os.makedirs(self.path, exist_ok=True)

        self.reset()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.it = 0
        self.t  = 0.0

        sep    = 2.5*self.radius
        mid_x  = self.d.x_c
        mid_y  = self.d.y_c
        half_x = 0.5*self.n_row*(self.radius + 0.5*sep)
        half_y = 0.5*self.n_col*(self.radius + 0.5*sep)

        for i in range(self.n_row):
            for j in range(self.n_col):
                self.p.x[self.n_col*i+j,0] = mid_x - half_x + sep*i + 0.5*self.radius*random.random()
                self.p.x[self.n_col*i+j,1] = mid_y - half_y + sep*j

    ### ************************************************
    ### Compute forces
    def forces(self):

        self.p.reset_forces()
        self.p.collisions(self.dt)
        self.d.collisions(self.p, self.dt)
        self.p.gravity(self.g)
