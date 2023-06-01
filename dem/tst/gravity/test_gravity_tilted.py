# Generic imports
import os
import pytest

# Custom imports
from dem.app.gravity  import *
from dem.src.core.run import *

###############################################
### Test gravity fall of one ball in a tilted rectangle
def test_gravity_tilted():

    # Initial space
    print("")

    # Run gravity app with tilted domain
    app = gravity(t_max           = 1.0,
                  dt              = 2.5e-5,
                  angle           = 40.0,
                  plot_show       = False,
                  plot_png        = False,
                  plot_trajectory = False)
    app.p.e_wall[:] = 0.5
    run(app)

    # Retrieve history
    x  = app.p.buff.get("x",  0)
    y  = app.p.buff.get("y",  0)
    vx = app.p.buff.get("vx", 0)
    vy = app.p.buff.get("vy", 0)
    t  = app.p.buff.get("t")

    # Check that final velocity is zero
    assert(abs(vx[-1]) < 1.0e-3)
    assert(abs(vy[-1]) < 1.0e-3)

    # Check final position
    assert(abs(x[-1] - 0.20195) < 1.0e-3)
    assert(abs(y[-1] - 0.03251) < 1.0e-3)
