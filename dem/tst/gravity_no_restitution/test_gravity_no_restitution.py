# Generic imports
import os
import pytest

# Custom imports
from dem.app.gravity  import *
from dem.src.core.run import *

###############################################
### Test gravity fall of one ball with no restitution
def test_gravity_no_restitution():

    # Initial space
    print("")

    # Run gravity app with restitution equal to 0
    app = gravity(t_max=0.4, dt=2.5e-5, plot_show=False, plot_trajectory=False)
    app.p.mtr.e_wall = 0.0
    run(app)

    # Retrieve history
    x  = app.p.buff.get("x",  0)
    y  = app.p.buff.get("y",  0)
    vx = app.p.buff.get("vx", 0)
    vy = app.p.buff.get("vy", 0)
    t  = app.p.buff.get("t")

    # Check final position
    assert(abs( x[-1] - 0.15      ) < 1.0e-5)
    assert(abs( y[-1] - app.p.r[0]) < 1.0e-5)

    # Check final velocity
    assert(abs(vx[-1]             ) < 1.0e-5)
    assert(abs(vy[-1]             ) < 5.0e-5)

    # Check time at which ball hit the ground
    # Index values where the ball is close to the ground
    yc = np.where(y < 1.05*app.p.r[0])[0][0]

    # Index values where the velocity is close to 0
    vc = np.where(abs(vy) < 1.0e-2)[0]

    # Exclude initial steps where the ball is at initial position
    ic = vc[vc > yc][0]

    # Get contact index and compute time
    tc     = t[ic]
    tc_ref = np.sqrt(2.0*(y[0] - app.p.r[0])/app.g)

    assert(abs(tc - tc_ref) < 1.0e-3)
