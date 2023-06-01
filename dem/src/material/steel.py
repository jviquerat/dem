from numba import float64
from numba.experimental import jitclass

### ************************************************
### Class defining steel material
### Taken from:
###   "DEM Investigation of the Influence of Particulate Properties
###    and Operating Conditions on the Mixing Process in Rotary Drums:
###    Part 1—Determination of the DEM Parameters and Calibration Process",
###    J. Hlosta et al, Processes, 8, 222 (2020)
spec = [
    ('density', float64),
    ('young',   float64),
    ('poisson', float64),
    ('mu_wall', float64),
    ('mu_part', float64),
    ('e_wall',  float64),
    ('e_part',  float64),
    ('Y',       float64),
    ('G',       float64)
]
@jitclass(spec)
class steel():
    ### ************************************************
    ### Constructor
    def __init__(self):

        # Density (kg/m3)
        self.density = 7850.0

        # Young modulus (N/m2)
        self.young   = 210.0e9

        # Poisson ratio (unitless)
        self.poisson = 0.3

        # Static friction coeff. on a wall of same material (unitless)
        self.mu_wall = 0.13

        # Static friction coeff. on a particle of same material (unitless)
        self.mu_part = 0.17

        # Restitution coeff. on a wall of same material (unitless)
        self.e_wall = 0.26

        # Restitution coeff. on a particle of same material (unitless)
        self.e_part = 0.85

        # Effective young modulus
        self.Y = (1.0 - self.poisson**2)/self.young

        # Effective shear modulus
        self.G = 2.0*(2.0 + self.poisson)*(1.0 - self.poisson)/self.young
