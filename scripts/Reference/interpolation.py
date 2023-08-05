

def lerp(y0, yf, x):
    """ Linear interpolation between y0 and yf. 
        x is in [0,1]
        y0, yf could be scalar or np.array
    """
    assert x >= 0 and x <= 1
    return y0 + (yf - y0) * x


def CubicBezier(y0, yf, x):
    """ Cubic Bezier interpolation between y0 and yf. 
        x is in [0,1]
        y0, yf could be scalar or np.array
    """
    assert x >= 0 and x <= 1
    yDiff = yf - y0
    bezier = x**3 + 3.0 * (x*x*(1.0 - x))
    return y0 + bezier*yDiff

def CubicBezierFirstDerivative(y0, yf, x):
    """ Derivative of Cubic Bezier interpolation between y0 and yf. 
        x is in [0,1]
        y0, yf could be scalar or np.array
    """
    assert x >= 0 and x <= 1
    yDiff = yf - y0
    bezier = 6.0 * x * (1.0 - x)
    return bezier * yDiff