import numpy as np
from math import pi, cos
import json

def initialize(params):
    return SourceTerm(json.loads(params))

class SourceTerm:
    def __init__(self, params):
        self._scale = params['scale']
        try:
            self._with33 = params["withtt"]
        except:
            self._with33 = False

    def setParams(self, mu):
        _mu = mu

    def eval(self, X):
        x = X[0]
        y = X[1]
        z = X[2]
        t = X[3]
        if z == 0.0:
            if self._with33:
                return self._scale * t * np.array([x*x, y*y, 3.0, x*y])
            else:
                return self._scale * t * np.array([x*x, y*y, x*y])
        else:
            return self._scale * t * np.array([x*x, y*y, z*z, x*y, y*z, x*z])

