# Generic imports
import numpy as np

###############################################
### Buffer class to store data during computation
class buff:
    def __init__(self, np, nt, names):

        self.np    = np
        self.nt    = nt
        self.names = names
        self.reset()

    def reset(self):

        self.data = {}
        for name in self.names:
            self.data[name] = np.zeros((self.np, self.nt))

        self.t = np.zeros((self.nt))
        self.k = 0

    def store(self, time, names, fields):

        n,f = zip(*[(i,j) for i,j in zip(names,fields) if i in self.names])
        for name, field in zip(n, f):
            self.data[name][:,self.k] = field[:]

        self.t[self.k] = time
        self.k        += 1

    def get(self, name, particle=None):

        if (name != "t"):
            return self.data[name][particle]
        else:
            return self.t
