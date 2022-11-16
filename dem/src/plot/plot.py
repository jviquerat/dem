# Generic imports
import os
import math
import numpy             as np
import matplotlib.pyplot as plt

from   matplotlib.patches import Rectangle, Circle

### ************************************************
### Output plot given domain and particles
def plot(d, p, path, it, show=False, png=False):

    # Plot counting and initialization
    if (it == 0):
        plt.ion()

    # Plot domain
    ax  = plt.gca()
    fig = plt.gcf()
    ax.set_xlim([d.x_min, d.x_max])
    ax.set_ylim([d.y_min, d.y_max])
    if (d.dtype == "rectangle"):
        ax.add_patch(Rectangle((d.x_min, d.y_min),
                                d.x_max-d.x_min,
                                d.y_max-d.y_min,
                                fill=False, color='r'))
    if (d.dtype == "circle"):
        ax.add_patch(Circle((0.5*(d.x_max+d.x_min),
                             0.5*(d.y_max+d.y_min)),
                             0.5*(d.x_max-d.x_min),
                             fill=False, color='r'))

    # Plot particles
    for i in range(p.n):
        ax.add_patch(Circle((p.x[i,0], p.x[i,1]), p.r[i],
                            fill=True, color='b'))

    ax.set_aspect('equal')
    fig.tight_layout()
    plt.grid()
    if png: fig.savefig(path+'/'+str(it)+'.png',
                        bbox_inches='tight')
    if show: plt.pause(0.0001)
    plt.clf()

### ************************************************
### Plot history of positions
def plot_history(n, h):

    plt.ioff()
    ax  = plt.gca()
    fig = plt.gcf()
    plt.plot(h[:,1])
    plt.grid()
    plt.show()
