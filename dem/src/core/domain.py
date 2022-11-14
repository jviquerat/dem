# Generic imports
import numpy as np

### ************************************************
### Class defining domain
class domain:
    ### ************************************************
    ### Constructor
    def __init__(self, dtype, x_min, x_max, y_min, y_max):

        self.dtype = dtype
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

    ### ************************************************
    ### Compute collision with a particle
    def collision(self, p):

        if (dtype == "rectangle"):
            pass

        if (dtype == "circle"):
            pass
