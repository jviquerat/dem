# Generic imports
import os
import pytest

# Custom imports
from dem.app.carreau  import *
from dem.src.core.run import *

###############################################
### Test carreau case between two particles
def test_carreau_nearest():

    # Initial space
    print("")

    ##################
    ### Run carreau app with nearest neighbor
    app = carreau(t_max           = 0.2,
                  dt              = 2.5e-5,
                  search          = "nearest",
                  plot_show       = False,
                  plot_png        = False,
                  plot_trajectory = False)
    app.p.e_wall[:] = 1.0
    app.p.e_part[:] = 1.0
    run(app)

    # Check that final velocity of first particle is zero
    vx0 = app.p.buff.get("vx", 0)
    vy0 = app.p.buff.get("vy", 0)
    assert(abs(vx0[-1]) < 1.0e-3)
    assert(abs(vy0[-1]) < 1.0e-8)

    # Check that final velocity of second particle is equal to first one
    vx1 = app.p.buff.get("vx", 1)
    assert(abs(abs(vx1[-1]) - abs(vx0[0])) < 3.0e-3)
