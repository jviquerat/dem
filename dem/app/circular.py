# Generic imports
import os
import math

# Custom imports
from dem.app.base_app import *

### ************************************************
### Several spheres falling on circular obstacles
class circular(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self,
                 name             = 'circular',
                 t_max            = 2.5,
                 dt               = 2.5e-5,
                 angle            = 0.0,
                 plot_freq        = 500,
                 plot_show        = True,
                 plot_trajectory  = True,
                 plot_png         = True,
                 plot_coll_radius = False):
        super().__init__()

        self.name             = name
        self.t_max            = t_max
        self.dt               = dt
        self.plot_freq        = plot_freq
        self.plot_show        = plot_show
        self.plot_trajectory  = plot_trajectory
        self.plot_png         = plot_png
        self.plot_coll_radius = plot_coll_radius

        self.nt        = int(self.t_max/self.dt)
        self.plot_it   = 0

        self.p = particles(np          = 16,
                           nt          = self.nt,
                           material    = "steel",
                           radius      = 0.025,
                           color       = "b",
                           store       = True,
                           search      = "nearest",
                           rad_coeff   = 2.0)

        colors = np.array(['r', 'g', 'b', 'c', 'm', 'y', 'k'])
        self.p.c = colors[np.random.randint(0,len(colors),size=self.p.np)]

        self.d = domain_factory.create("rectangle",
                                       x_min      = 0.0,
                                       x_max      = 1.0,
                                       y_min      = 0.0,
                                       y_max      = 2.0,
                                       angle      = 0.0,
                                       material   = "steel")

        self.o0 = domain_factory.create("circle",
                                        x_c       = 0.25,
                                        y_c       = 0.5,
                                        rad       = 0.2,
                                        plot_fill = True,
                                        material  = "steel")
        self.o1 = domain_factory.create("circle",
                                        x_c       = 0.75,
                                        y_c       = 0.5,
                                        rad       = 0.2,
                                        plot_fill = True,
                                        material  = "steel")

        self.d_lst = [self.d, self.o0, self.o1]

        self.path = self.base_path+'/'+self.name
        os.makedirs(self.path, exist_ok=True)

        self.reset()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.it = 0
        self.t  = 0.0

        for i in range(self.p.np):
            self.p.x[i,0] = 0.025*(1.0 + 2.5*i)
            self.p.x[i,1] = 1.5

    ### ************************************************
    ### Compute forces
    def forces(self):

        self.p.reset_forces()
        self.p.collisions(self.dt)
        for d in self.d_lst:
            d.collisions(self.p, self.dt)
        self.p.gravity(self.g)
