# Generic imports
import os
import pytest

# Custom imports
from dem.app.gravity  import *
from dem.src.core.run import *

###############################################
### Test gravity fall of one ball with full restitution
def test_gravity_full_restitution():

    # Initial space
    print("")

    # Run gravity app with restitution equal to 1
    app = gravity(t_max           = 0.6,
                  dt              = 2.5e-5,
                  angle           = 0.0,
                  plot_show       = False,
                  plot_png        = False,
                  plot_trajectory = False)
    app.p.e_wall[:] = 1.0
    run(app)

    # Retrieve history
    x  = app.p.buff.get("x",  0)
    y  = app.p.buff.get("y",  0)
    vx = app.p.buff.get("vx", 0)
    vy = app.p.buff.get("vy", 0)
    t  = app.p.buff.get("t")

    # Check height of ball after first bounce
    # Index values where the ball is close to the ground
    yc = np.where(y < 1.05*app.p.r[0])[0][0]

    # Index values where the ball is close to initial height
    yt = np.where(abs(y - y[0]) < 1.0e-3)[0]

    # Exclude initial steps where the ball is at initial position
    it = yt[yt > yc][0]

    # Find max altitudes for i>it
    ht = np.amax(y[it:])

    assert(abs(ht - y[0]) < 1.0e-3)
