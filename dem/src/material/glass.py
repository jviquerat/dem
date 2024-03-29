from numba import float64
from numba.experimental import jitclass

### ************************************************
### Class defining glass material
### Taken from:
###   "Validating N-body code CHRONO for granular DEM simulations
###    in reduced-gravity environments",
###    C. Sunday et al, Monthly Notices of the Royal Astronomical Society,
###    498, 1062-1079 (2020)
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
class glass():
    ### ************************************************
    ### Constructor
    def __init__(self):

        # Density (kg/m3)
        self.density = 2500

        # Young modulus (N/m2)
        self.young   = 70.0e9

        # Poisson ratio (unitless)
        self.poisson = 0.24

        # Static friction coeff. on a wall of same material (unitless)
        self.mu_wall = 0.45

        # Static friction coeff. on a particle of same material (unitless)
        self.mu_part = 0.16

        # Restitution coeff. on a wall of same material (unitless)
        self.e_wall = 0.82

        # Restitution coeff. on a particle of same material (unitless)
        self.e_wall = 0.97

        # Effective young modulus
        self.Y = (1.0 - self.poisson**2)/self.young

        # Effective shear modulus
        self.G = 2.0*(2.0 + self.poisson)*(1.0 - self.poisson)/self.young
