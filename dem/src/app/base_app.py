# Generic imports
import math

# Custom imports
from dem.src.core.particle import *
from dem.src.core.domain   import *

### ************************************************
### Base app
class base_app():

    ### ************************************************
    ### Constructor
    def __init__(self):

        self.plt_show = False
        self.plt_png  = False

    ### ************************************************
    ### Iteration printings
    def printings(self, it):

        print("# it = "+str(it)+", t = {0:.3f}".format(self.t)+" / {0:.3f}".format(self.t_max), end='\r')

    ### ************************************************
    ### Check stopping criterion
    def check_stop(self):

        compute = True
        if (self.t >= self.t_max):
            compute = False
            print('\n')
            print('# Computation ended: t>t_max')

        return compute

    ### ************************************************
    ### Compute forces
    def forces(self):

        raise NotImplementedError

    ### ************************************************
    ### Update positions
    def update(self):

        raise NotImplementedError

    ### ************************************************
    ### Plot
    def plot(self, it):

        raise NotImplementedError

    ### ************************************************
    ### Finalize
    def finalize(self):

        raise NotImplementedError
