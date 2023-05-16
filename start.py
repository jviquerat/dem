# Generic imports
import os
import sys

# Custom imports
from dem.app.app      import *
from dem.src.core.run import *

########################
# Run dem simulation
########################
if __name__ == '__main__':

    # Check command-line input
    if (len(sys.argv) == 2):
        app_name = sys.argv[1]
    else:
        print('Command line error, please use as follows:')
        print('python3 start.py app_name')

    # Instanciate app
    app = app_factory.create(app_name)

    # Run
    run(app)
