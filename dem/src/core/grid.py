# Generic imports
import math
import numpy as np

# Custom imports
from dem.src.core.collision    import *
from dem.src.material.material import *
from dem.src.utils.buff        import *

### ************************************************
### Class defining grid for neighboring cell algorithm
class grid:
    ### ************************************************
    ### Constructor
    def __init__(self,
                 domain,   # reference domain
                 nx = 10,  # number of bins in x direction
                 ny = 10): # number of bins in y direction

        self.x_min = domain.x_min
        self.x_max = domain.x_max
        self.y_min = domain.y_min
        self.y_max = domain.y_max
        self.nx    = nx
        self.ny    = ny
        self.n     = nx*ny

        # Minimal number of bins per direction is 2
        if (nx < 2): nx = 2
        if (ny < 2): ny = 2

        self.x_bins    = np.linspace(self.x_min, self.x_max, self.nx+1)
        self.y_bins    = np.linspace(self.y_min, self.y_max, self.ny+1)
        self.dx        = (self.x_max-self.x_min)/float(self.nx+1)
        self.dy        = (self.y_max-self.y_min)/float(self.ny+1)
        self.dist      = min(self.dx, self.dy)
        self.counter   = 0
        self.frequency = 200

        # Create list of cell neighbors
        self.build_neighbors()

    ### ************************************************
    ### Set particles array
    def set(self, x):

        # Sort particles in bins
        dig_x = np.digitize(x[:,0], self.x_bins)
        dig_y = np.digitize(x[:,1], self.y_bins)

        # Build arrays:
        # - parts[i][j] returns the list of particle
        #   indexes contained in cell (i,j)
        # - cells[k,:] returns the indexes (i,j) of
        #   the cell containing particle k
        self.parts = np.empty((self.nx, self.ny), object)
        self.cells = np.zeros((len(x),  2),       np.int16)
        for i in range(self.nx):
            for j in range(self.ny):
                self.parts[i,j] = []

        for k in range(len(x)):
            i = dig_x[k]-1
            j = dig_y[k]-1
            self.parts[i,j].append(k)
            self.cells[k,0] = i
            self.cells[k,1] = j

            # Reset counter
            self.counter = 0

    ### ************************************************
    ### Correct particles array
    def correct(self, x):

        if (self.counter%self.frequency == 0):

            for i in range(self.nx):
                for j in range(self.ny):
                    for p in self.parts[i,j]:
                        is_in = self.is_in_cell(i, j, x[p,:])
                        if (not is_in):
                            self.parts[i,j].remove(p)
                            ngb = self.ngb[i,j]
                            for n in range(len(ngb)):
                                k = ngb[n][0]
                                l = ngb[n][1]
                                is_in = self.is_in_cell(k, l, x[p,:])
                                if (is_in):
                                    self.parts[k,l].append(p)
                                    self.cells[p,0] = k
                                    self.cells[p,1] = l

            self.counter = 0

        else:

            # Increment counter
            self.counter += 1

    ### ************************************************
    ### Check if coordinate is in cell
    def is_in_cell(self, i, j, x):

        if ((x[0] > i*self.dx)     and
            (x[0] < (i+1)*self.dx) and
            (x[1] > j*self.dy)     and
            (x[1] < (j+1)*self.dy)): return True

        return False

    ### ************************************************
    ### Return list of neighbors of given cell
    ### The list also includes current cell for direct
    ### research of colliding particles
    def build_neighbors(self):

        self.ngb = np.empty((self.nx,self.ny), object)

        for i in range(self.nx):
            for j in range(self.ny):
                self.ngb[i,j] = []

                # Bottom left
                if (i==0 and j==0):
                    self.append_ngb(i, j, [0,1], [0,1])

                # Left
                elif (i==0 and j>0 and j<self.ny-1):
                    self.append_ngb(i, j, [0,1], [-1,0,1])

                # Top left
                elif (i==0 and j==self.ny-1):
                    self.append_ngb(i, j, [0,1], [-1,0])

                # Top
                elif (i>0 and i<self.nx-1 and j==self.ny-1):
                    self.append_ngb(i, j, [-1,0,1], [-1,0])

                # Top right
                elif (i==self.nx-1 and j==self.ny-1):
                    self.append_ngb(i, j, [-1,0], [-1,0])

                # Right
                elif (i==self.nx-1 and j<self.ny-1 and j>0):
                    self.append_ngb(i, j, [-1,0], [-1,0,1])

                # Bottom right
                elif (i==self.nx-1 and j==0):
                    self.append_ngb(i, j, [-1,0], [0,1])

                # Bottom
                elif (i<self.nx-1 and i>0 and j==0):
                    self.append_ngb(i, j, [-1,0,1], [0,1])

                # General case
                else:
                    self.append_ngb(i, j, [-1,0,1], [-1,0,1])

    ### ************************************************
    ### Append neighbors to list
    def append_ngb(self, i, j, lst_i, lst_j):

        for k in lst_i:
            for l in lst_j:
                self.ngb[i,j].append([i+k, j+l])

    ### ************************************************
    ### Check if grid needs to be reset
    def reset(self):

        self.counter += 1
        if (self.counter%self.frequency == 0):
            self.counter = 0
            return True
        else:
            return False

    ### ************************************************
    ### Return cell id from coordinates
    def id(self, i, j):

        return j*self.nx+i
