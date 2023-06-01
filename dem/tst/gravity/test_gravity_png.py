# Generic imports
import os
import pytest

# Custom imports
from dem.app.gravity  import *
from dem.src.core.run import *

###############################################
### Test that the png plotting does not fail
def test_gravity_png():

    # Initial space
    print("")

    # Run gravity app
    app = gravity(name            = 'gravity_test',
                  t_max           = 1.0,
                  dt              = 2.5e-5,
                  angle           = 40.0,
                  plot_show       = False,
                  plot_png        = True,
                  plot_freq       = 1000,
                  plot_trajectory = False)
    app.p.e_wall[:] = 0.5
    run(app)
