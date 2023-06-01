# Generic imports
import math

# Custom imports
from dem.src.core.particles    import *
from dem.src.domain.domain     import *
from dem.src.material.material import *
from dem.src.utils.plot        import *

### ************************************************
### Base app
class base_app():

    ### ************************************************
    ### Constructor
    def __init__(self):

        self.plot_coll_radius = False
        self.base_path        = 'results'
        self.g                = 9.81

    ### ************************************************
    ### Iteration printings
    def printings(self):

        print("# it = "+str(self.it)+", t = {0:.3f}".format(self.t)+" / {0:.3f}".format(self.t_max), end='\r')

    ### ************************************************
    ### Check stopping criterion
    def check_stop(self):

        compute = True
        if (self.it >= self.nt):
            compute = False
            print('\n')
            print('# Computation ended: it>=nt')

        return compute

    ### ************************************************
    ### Compute forces
    def forces(self):

        raise NotImplementedError

    ### ************************************************
    ### Update positions
    def update(self):

        self.p.update(self.t, self.dt, self.it)
        self.it += 1
        self.t  += self.dt

    ### ************************************************
    ### Plot
    def plot(self):

        if (self.plot_show or self.plot_png):
            if (self.it%self.plot_freq == 0):
                plot(self)
                self.plot_it += 1

    ### ************************************************
    ### Finalize
    def finalize(self):

        if (self.plot_trajectory):
            plot_trajectory(self.p.np, self.p.buff, self.p.c)
