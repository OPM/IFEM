from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from math import *

SetFinalOutput("annulus.g2")

r1 = 5
r2 = 7

inner = CircleSegment([0,0,0],[-r1,0,0],pi,[0,0,1])
outer = CircleSegment([0,0,0],[-r2,0,0],pi,[0,0,1])

surf = LoftCurves([inner,outer])

FinalOutput(surf,True)
