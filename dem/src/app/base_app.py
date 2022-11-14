# Generic imports
import math

# Custom imports
from dem.src.core.particle import *
from dem.src.core.domain   import *

###############################################
### Base app
class base_app():

    ### Iteration printings
    def printings(self, it):

        print('# it = '+str(it)+' / '+str(self.it_max), end='\r')

    ### Check stopping criterion
    def check_stop(self, it):

        compute = True
        if (it >= self.it_max):
            compute = False
            print('\n')
            print('# Computation ended: it>it_max')

        return compute
