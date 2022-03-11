import numpy as np
from math import pi, cos
import json

def initialize(params):
    return SourceTerm(json.loads(params))

class SourceTerm:
    def __init__(self, params):
        self._scale = params['scale']

    def setParams(self, mu):
        _mu = mu

    def eval(self, X):
        x = X[0]
        y = X[1]
        z = X[2]
        t = X[3]
        if z == 0.0:
            return self._scale * t * np.array([1.0*x, 2.0*x, 1.0*y, 2.0*y])
        else:
            return self._scale * t * np.array([1.0*x, 2.0*x, 3.0*x, 1.0*y, 2.0*y, 3.0*y, 1.0*z, 2.0*z, 3.0*z])
