# Generic imports
import os
import pytest

# Custom imports
from dem.app.obstacle  import *
from dem.src.core.run import *

###############################################
### Test collision with obstacle corner
def test_obstacle_corner():

    # Initial space
    print("")

    # Run gravity app with tilted domain
    app = obstacle(t_max            = 0.8,
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
    app.o0 = domain_factory.create("rectangle",
                                   x_min      = 0.0,
                                   x_max      = 0.14,
                                   y_min      = 0.1,
                                   y_max      = 0.24,
                                   angle      = 45.0,
                                   plot_fill  = True,
                                   material   = "steel")
    app.d_lst = [app.d, app.o0]
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
    assert(abs(x[-1] - 0.07) < 1.0e-3)
    assert(abs(y[-1] - 0.28899) < 1.0e-3)
