# Generic imports
import math

# Custom imports
from dem.src.app.base_app import *

### ************************************************
### Gravity fall
class gravity(base_app):
    ### ************************************************
    ### Constructor
    def __init__(self):

        self.name        = 'gravity'
        self.t_max       = 1.0
        self.dt          = 0.1
        self.it_max      = math.floor(self.t_max/self.dt)
        self.p = particles(1, 1.0, 0.01)
        self.d = domain("rectangle", 0.0, 1.0, 0.0, 1.0)

        self.reset()

    ### ************************************************
    ### Reset app
    def reset(self):

        self.p.x[0,0] = 0.5
        self.p.x[0,1] = 0.8

