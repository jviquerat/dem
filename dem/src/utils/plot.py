# Generic imports
import os
import math
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.patches     import Rectangle, Circle
from matplotlib.collections import PatchCollection

### ************************************************
### Output plot of an app
def plot(app):

    # Plot counting and initialization
    if (app.plot_it == 0):
        plt.ion()

    # Plot domain
    ax  = plt.gca()
    fig = plt.gcf()
    ax.set_xlim([app.d_lst[0].x_min, app.d_lst[0].x_max])
    ax.set_ylim([app.d_lst[0].y_min, app.d_lst[0].y_max])
    ax.set_axis_off()
    fig.tight_layout()
    plt.margins(0,0)
    plt.subplots_adjust(0,0,1,1,0,0)
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())

    patches = []

    # Plot main domain
    for d in app.d_lst:
        patches.append(Rectangle(d.p1, d.dx, d.dy,
                                 angle          = d.angle,
                                 rotation_point = 'xy',
                                 fill           = d.plot_fill))

    # Plot particles
    for i in range(app.p.np):
        patches.append(Circle((app.p.x[i,0], app.p.x[i,1]), app.p.r[i],
                              fill=True, color=app.p.c[i]))
        #if (hasattr(app.p, 'm_rad')):
        #    patches.append(Circle((app.p.x[i,0], app.p.x[i,1]), app.p.m_rad,
        #                          fill=False, color=app.p.c[i]))

    col = PatchCollection(patches, match_original=True)
    ax.add_collection(col)
    ax.set_aspect('equal')

    #plt.grid()
    if app.plot_png: fig.savefig(app.path+'/'+str(app.plot_it)+'.png',
                                 bbox_inches='tight', pad_inches=0)
    if app.plot_show: plt.pause(0.01)
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
