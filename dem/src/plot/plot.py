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
                            fill=True, color=p.c[i]))

    ax.set_aspect('equal')
    fig.tight_layout()
    plt.grid()
    if png: fig.savefig(path+'/'+str(it)+'.png',
                        bbox_inches='tight')
    if show: plt.pause(0.0001)
    plt.clf()

### ************************************************
### Plot history of positions
def plot_history(n, h, c):

    h = h.reshape((-1,2*n))
    plt.ioff()
    ax  = plt.gca()
    fig = plt.gcf()
    for i in range(n):
        ax.add_patch(Circle((h[0,2*i],h[0,2*i+1]), 0.01, color=c[i]))
        plt.plot(h[:,2*i],h[:,2*i+1], color=c[i], linestyle='dashed')
    ax.set_aspect('equal')
    fig.tight_layout()
    plt.grid()
    plt.show()
