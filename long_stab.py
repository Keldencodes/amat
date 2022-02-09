import matplotlib.pyplot as plt
import numpy as np
import conversions as c
import os
import pyxfoil as p

# define some constants

# distances in inches
xcg = 72 # center of gravity location from leading edge of fuselage
xc = 4
xf = 36
xw = 66

# spans in inches
bc = 24 # canard
bw = 95 # wing
bf = 24 # fuselage

# chords in inches
cc = 6
cf = 55
cw = 12

u_inf = 20 # mph

rho = 1.225 # kg/m3
mu = 1.802e-5 # kg/m*s

W = 55 # lbs

airfoil_c = 'naca0012'
airfoil_f = 'n24'
airfoil_w = 'naca2415'

def can_alpha_req(xcg, xc, xf, xw, bc, bw, bf, cc, cf, cw, airfoil_c, airfoil_f, airfoil_w, rho, mu, W, u_inf):
    # calculate required alpha of canard with varying fuselage/wing alpha
    # conversions
    cc = c.in2m(cc)
    cf = c.in2m(cf)
    cw = c.in2m(cw)
    
    bc = c.in2m(bc)
    bf = c.in2m(bf)
    bw = c.in2m(bw)
    
    xcg = c.in2m(xcg)
    xc = c.in2m(xc)
    xf = c.in2m(xf)
    xw = c.in2m(xw)
    
    u_inf = c.ftpers2mpers(c.mph2ftpers(u_inf))
    
    # surface areas
    Sc = cc * bc
    Sf = cf * bf
    Sw = cw * bw
    
    # reynolds of each 'wing'
    Re_c = rho * u_inf * cc / mu
    Re_f = rho * u_inf * cf / mu
    Re_w = rho * u_inf * cw / mu
    
    airfoil_c_path = os.path.join('data', airfoil_c)
    with open(airfoil_c_path, 'r') as infile:
        x_c, y_c = np.loadtxt(infile, unpack=True, skiprows=1)
        
    airfoil_f_path = os.path.join('data', airfoil_f)
    with open(airfoil_f_path, 'r') as infile:
        x_f, y_f = np.loadtxt(infile, unpack=True, skiprows=1)
    
    airfoil_w_path = os.path.join('data', airfoil_w)
    with open(airfoil_w_path, 'r') as infile:
        x_w, y_w = np.loadtxt(infile, unpack=True, skiprows=1)    
    
    
    
def total_lift(bc, bw, bf, u_inf, alf_main, alf_c):
    # calculate total lift of aircraft
    Sc = c.in2m(cc) * c.in2m(bc)
    Sf = c.in2m(cf) * c.in2m(bf)
    Sw = c.in2m(cw) * c.in2m(bw)
    
    
    
    
can_alpha_req(xcg, xc, xf, xw, bc, bw, bf, cc, cf, cw, airfoil_c, airfoil_f, airfoil_w, rho, mu, W, u_inf)