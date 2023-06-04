# Generic imports
import os
import pytest

# Custom imports
from dem.app.obstacle import *
from dem.src.core.run import *

###############################################
### Test presence of multiple obstacles
def test_multiple_obstacles():

    # Initial space
    print("")

    # Run gravity app with tilted domain
    app = obstacle(t_max            = 1.2,
                   dt               = 2.5e-5,
                   angle            = 0.0,
                   plot_show        = False,
                   plot_png         = False,
                   plot_trajectory  = False,
                   plot_coll_radius = False)

    app.p = particles(np          = 1,
                      nt          = app.nt,
                      material    = "steel",
                      radius      = 0.02,
                      color       = "b",
                      store       = True,
                      search      = "nearest",
                      rad_coeff   = 2.0)
    app.p.x[0,0]    = 0.07
    app.p.x[0,1]    = 0.4
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
    assert(abs(x[-1] - 0.059229) < 1.0e-3)
    assert(abs(y[-1] - 0.019999) < 1.0e-3)
