from this import d


def ftpers2mph(v):
    return v*0.681818

def mph2ftpers(v):
    return v/0.681818

def ftpers2mpers(v):
    return v*0.3048

def mpers2ftpers(v):
    return v/0.3048

def mph2mpers(v):
    return v*0.44704

def mpers2mph(v):
    return v/0.44704

def in2m(l):
    return l*0.0254

def m2in(l):
    return l/0.0254

def interpolate(x1, y1, x2, y2, x3):
    return (y2-y1)*(x3-x1)/(x2-x1)+y1