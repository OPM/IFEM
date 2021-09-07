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

    def eval(self, x):
        return self._scale * x
