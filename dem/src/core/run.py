# Generic imports
import os
import time

# Custom imports
from dem.src.app.app   import *


########################
# Run dem simulation
########################
def run(app):

    # Timer and loop data
    start_time = time.time()
    it         = 0
    compute    = True

    # Solve
    print('### Solving')
    while (compute):

        # Printings and plot
        app.printings(it)
        app.plot(it)

        # Compute forces
        app.forces()

        # Update positions
        app.update(it)

        # Check stopping criterion
        compute = app.check_stop()

        # Increment iteration
        it += 1

    # Count time
    end_time = time.time()
    print("# Loop time = {:f}".format(end_time - start_time))

    # Perform final operations and outputs
    app.finalize()

