from this import d


def ftpers2mph(v):
    return v*0.681818

def mph2ftpers(v):
    return v*0.681818

def ftpers2mpers(v):
    return v*0.3048

def mpers2ftpers(v):
    return v/0.3048

def mph2mpers(v):
    return ftpers2mpers(mph2ftpers(v))

def mpers2mph(v):
    return ftpers2mph(mpers2ftpers(v))

def in2m(l):
    return l*0.0254

def m2in(l):
    return l/0.0254