import numpy as np
from scipy.optimize import fsolve

def solve(x):

    x3 = x[0]
    x2 = 0.5
    x1 = 0.5
    while True:
        temp_x2 = x2
        temp_x1 = x1
        x0 = x2
        x2 = ((x0 - x1**2)/x3)**(1/3)
        x1 = x2**0.5 - x0
        if abs((x2 - temp_x2) / (temp_x2 + 1e-08)) > 0.0001: continue
        if abs((x1 - temp_x1) / (temp_x1 + 1e-08)) > 0.0001: continue
        break
    y = [0.0] * 3
    y[0] = x2
    y[1] = x1
    y[2] = x0
    return y