# Generic imports
import os
import math
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.patches     import Rectangle, Circle
from matplotlib.collections import PatchCollection

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

    # Plot particles
    patches = []
    for i in range(p.np):
        patches.append(Circle((p.x[i,0], p.x[i,1]), p.r[i],
                              fill=True, color=p.c[i]))
        #if (hasattr(p, 'm_rad')):
        #    patches.append(Circle((p.x[i,0], p.x[i,1]), p.m_rad,
        #                          fill=False, color=p.c[i]))

    col = PatchCollection(patches, match_original=True)
    ax.add_collection(col)

    ax.set_aspect('equal')
    fig.tight_layout()
    #plt.grid()
    if png: fig.savefig(path+'/'+str(it)+'.png',
                        bbox_inches='tight')
    if show: plt.pause(0.0001)
    plt.clf()

### ************************************************
### Plot trajectory
def plot_trajectory(np, h, c):

    hx = h.data["x"]
    hy = h.data["y"]

    plt.ioff()
    ax  = plt.gca()
    fig = plt.gcf()

    for i in range(np):
        ax.add_patch(Circle((hx[i,0], hy[i,0]), 0.01, color=c[i]))
        plt.plot(hx[i,:], hy[i,:], color=c[i], linestyle='dashed')

    ax.set_aspect('equal')
    fig.tight_layout()
    plt.grid()
    plt.show()
