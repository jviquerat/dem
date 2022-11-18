# Generic imports
import math

### ************************************************
### Base domain
class base_domain():

    ### ************************************************
    ### Constructor
    def __init__(self):

        pass

    ### ************************************************
    ### Compute distance to particle
    def distance(self, p):

        raise NotImplementedError

    ### ************************************************
    ### Detect collision with particle
    def collide(self, p):

        raise NotImplementedError
