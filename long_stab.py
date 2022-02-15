import matplotlib.pyplot as plt
import numpy as np
import conversions as c
import os
import pyxfoil as p

# define some constants
# spans in inches
bc = 24 # canard
bw = 90 # wing
bf = 22 # fuselage with 2 columns of soccer balls

# chords in inches
cc = 12
cf = 10/.171 # using 1 row of soccer balls to define fuselage size
cw = 20

# distances in inches to xac and xcg
xc = cc/4 + 3
xf = cf/4
xw = cf - 3/4*cw # place trailing edge of wing coincident with TE of fuselage
xcg = 27.5 # center of gravity location from leading edge of fuselage
# should probably define the cg later using the xac=xcg large xbarac function
u_inf = 25 # mph

rho = 1.225 # kg/m3
mu = 1.802e-5 # kg/m*s

W = 55 # lbs

airfoil_c = 'naca0012'
airfoil_f = 'n24'
airfoil_w = 'naca2415'

def equilibrium(xcg, xc, xf, xw, bc, bw, bf, cc, cf, cw, airfoil_c, airfoil_f, airfoil_w, rho, mu, W, u_inf):
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
    
    u_inf = c.mph2mpers(u_inf)
    
    # surface areas in m2
    Sc = cc * bc
    Sf = cf * bf
    Sw = cw * bw
    
    # reynolds of each 'wing'
    Re_c = rho * u_inf * cc / mu
    Re_f = rho * u_inf * cf / mu
    Re_w = rho * u_inf * cw / mu
    
    airfoil_c_path = os.path.join('data', f'{airfoil_c}.dat')
    with open(airfoil_c_path, 'r') as infile:
        x_c, y_c = np.loadtxt(infile, unpack=True, skiprows=1)
        
    airfoil_f_path = os.path.join('data', f'{airfoil_f}.dat')
    with open(airfoil_f_path, 'r') as infile:
        x_f, y_f = np.loadtxt(infile, unpack=True, skiprows=1)
    
    airfoil_w_path = os.path.join('data', f'{airfoil_w}.dat')
    with open(airfoil_w_path, 'r') as infile:
        x_w, y_w = np.loadtxt(infile, unpack=True, skiprows=1)    
    
    #running xfoil
    start_c = -10
    stop_c = 15
    num_c = (stop_c-start_c)*10 + 1
    alphas_c = np.linspace(start=start_c, stop=stop_c, num=num_c)
    start_f = -10
    stop_f = 25
    num_f = (stop_f-start_f)*10 + 1
    alphas_f = np.linspace(start=start_f, stop=stop_f, num=num_f)
    start_w = -10
    stop_w = 25
    num_w = (stop_w-start_w)*10 + 1
    alphas_w = np.linspace(start=start_w, stop=stop_w, num=num_w)
    naca = False
    
    f_c = f'data/{airfoil_c}/{airfoil_c}_polar_Re{Re_c:.2e}a{start_c:.1f}-{stop_c:.1f}.dat'
    f_f = f'data/{airfoil_f}/{airfoil_f}_polar_Re{Re_f:.2e}a{start_f:.1f}-{stop_f:.1f}.dat'
    f_w = f'data/{airfoil_w}/{airfoil_w}_polar_Re{Re_w:.2e}a{start_w:.1f}-{stop_w:.1f}.dat'

    # if the xfoil data already exists, don't run it again
    if not os.path.exists(f_c):
        foil = f'data/{airfoil_c}.dat'
        p.GetPolar(foil, naca, alphas_c, Re_c, pane=True)
        print('did xfoil')

    if not os.path.exists(f_f):
        foil = f'data/{airfoil_f}.dat'
        p.GetPolar(foil, naca, alphas_f, Re_f, pane=True)
        print('did xfoil')

    if not os.path.exists(f_w):
        foil = f'data/{airfoil_w}.dat'
        p.GetPolar(foil, naca, alphas_w, Re_w, pane=True)
        print('did xfoil')
  
    canard_data = np.loadtxt(f_c, unpack=True, skiprows=12)
    fuselage_data = np.loadtxt(f_f, unpack=True, skiprows=12)
    wing_data = np.loadtxt(f_w, unpack=True, skiprows=12)

    # plot cl vs a for each
    plt.figure(figsize=(10,10))
    plt.plot(canard_data[0], canard_data[1], 'r', label=f'Canard: {airfoil_c}')
    plt.plot(fuselage_data[0], fuselage_data[1], 'g', label=f'Fuselage: {airfoil_f}')
    plt.plot(wing_data[0], wing_data[1], 'b', label=f'Wing: {airfoil_w}')
    plt.legend()
    plt.xlabel('alpha')
    plt.ylabel('Cl')
    plt.show()

    alpha_c = 9
    result = np.where(canard_data[0]==alpha_c)
    if result[0].size==0:
       alpha1_index = max(np.where(canard_data[0]<alpha_c)[0])
       alpha2_index = alpha1_index + 1
       alpha1 = canard_data[0][alpha1_index]
       alpha2 = canard_data[0][alpha2_index]
       cl1 = canard_data[1][alpha1_index]
       cl2 = canard_data[1][alpha2_index]
       clc = c.interpolate(x1=alpha1, y1=cl1, x2=alpha2, y2=cl2, x3=alpha)
    else:
       clc = canard_data[1][result[0][0]]

    alpha = 14.0
    result = np.where(fuselage_data[0]==alpha)
    if result[0].size==0:
       alpha1_index = max(np.where(fuselage_data[0]<alpha)[0])
       alpha2_index = alpha1_index + 1
       alpha1 = fuselage_data[0][alpha1_index]
       alpha2 = fuselage_data[0][alpha2_index]
       cl1 = fuselage_data[1][alpha1_index]
       cl2 = fuselage_data[1][alpha2_index]
       clf = c.interpolate(alpha1, cl1, alpha2, cl2, alpha)
    else:
       clf = fuselage_data[1][result[0][0]]
 
    result = np.where(wing_data[0]==alpha)
    if result[0].size==0:
       alpha1_index = max(np.where(wing_data[0]<alpha)[0])
       alpha2_index = alpha1_index + 1
       alpha1 = wing_data[0][alpha1_index]
       alpha2 = wing_data[0][alpha2_index]
       cl1 = wing_data[1][alpha1_index]
       cl2 = wing_data[1][alpha2_index]
       clw = c.interpolate(alpha1, cl1, alpha2, cl2, alpha)
    else:
       clw = wing_data[1][result[0][0]]
      

    Lc = .5 * clc * rho * (u_inf**2) * Sc
    Lw = .5 * clw * rho * (u_inf**2) * Sw
    Lf = .5 * clf * rho * (u_inf**2) * Sf

    l = (rho*(u_inf**2))/(2*(xcg-xc)) * (clw*Sw*(xw-xc) + clf*Sf*(xf-xc))
    print(f'Lift: {c.n2pound(l):.2f}')

    print(f'*********************************************** \n Results')
    print(f'Relative airspeed = {c.mpers2mph(u_inf):.1f}mph')
    print(f'Xcg = {c.m2in(xcg):.2f}')
    print(f'Fuselage size \n Chord = {(c.m2in(cf)):.1f}in \n Span = {(c.m2in(bf)):.1f}in \n Area = {(c.m2in(cf)*c.m2in(bf)):.2f}in2')
    print(f'Wing size \n Chord = {(c.m2in(cw)):.1f}in \n Span = {(c.m2in(bw)):.1f}in \n Area = {(c.m2in(cw)*c.m2in(bw)):.2f}in2')
    print(f'Canard size \n Chord = {(c.m2in(cc)):.1f}in \n Span = {(c.m2in(bc)):.1f}in \n Area = {(c.m2in(cc)*c.m2in(bc)):.2f}in2')
    print(f'Canard lift: Cl = {clc:.3f}, L = {Lc:.2f}N = {c.n2pound(Lc):.2f}lbs')
    print(f'Wing lift: CL = {clw:.3f}, L = {Lw:.2f}N = {c.n2pound(Lw):.2f}lbs')
    print(f'Fuselage lift: Cl = {clf:.3f}, L = {Lf:.2f}N = {c.n2pound(Lf):.2f}lbs')
    print(f'Total lift: {(Lc+Lw+Lf):.2f}N, {c.n2pound(Lc+Lw+Lf):.2f}lbs')
    
    
def total_lift(alpha_c, alpha_wf, bc, bw, bf, cc, cf, cw, airfoil_c, airfoil_f, airfoil_w, rho, mu, W, u_inf):
    # calculate total lift at given angles of attack
    # conversions
    cc = c.in2m(cc)
    cf = c.in2m(cf)
    cw = c.in2m(cw)
    
    bc = c.in2m(bc)
    bf = c.in2m(bf)
    bw = c.in2m(bw)
    
    u_inf = c.ftpers2mpers(c.mph2ftpers(u_inf))
    
    # surface areas
    Sc = cc * bc
    Sf = cf * bf
    Sw = cw * bw
    
    # reynolds of each 'wing'
    Re_c = rho * u_inf * cc / mu
    Re_f = rho * u_inf * cf / mu
    Re_w = rho * u_inf * cw / mu
    
    airfoil_c_path = os.path.join('data', f'{airfoil_c}.dat')
    with open(airfoil_c_path, 'r') as infile:
        x_c, y_c = np.loadtxt(infile, unpack=True, skiprows=1)
        
    airfoil_f_path = os.path.join('data', f'{airfoil_f}.dat')
    with open(airfoil_f_path, 'r') as infile:
        x_f, y_f = np.loadtxt(infile, unpack=True, skiprows=1)
    
    airfoil_w_path = os.path.join('data', f'{airfoil_w}.dat')
    with open(airfoil_w_path, 'r') as infile:
        x_w, y_w = np.loadtxt(infile, unpack=True, skiprows=1)    
    
    #running xfoil
    start = 0
    stop = 20
    num = 50
    alphas = np.linspace(start=start, stop=stop, num=num)
    naca = False
    
    f_c = f'data/{airfoil_c}/{airfoil_c}_polar_Re{Re_c:.2e}a{start:.1f}-{stop:.1f}.dat'
    f_f = f'data/{airfoil_f}/{airfoil_f}_polar_Re{Re_f:.2e}a{start:.1f}-{stop:.1f}.dat'
    f_w = f'data/{airfoil_w}/{airfoil_w}_polar_Re{Re_w:.2e}a{start:.1f}-{stop:.1f}.dat'

    if not os.path.exists(f_c):
        foil = f'data/{airfoil_c}.dat'
        p.GetPolar(foil, naca, alphas, Re_c, pane=True)
        print(f'did xfoil on {airfoil_c}')

    if not os.path.exists(f_f):
        foil = f'data/{airfoil_f}.dat'
        p.GetPolar(foil, naca, alphas, Re_f, pane=True)
        print(f'did xfoil on {airfoil_f}')

    if not os.path.exists(f_w):
        foil = f'data/{airfoil_w}.dat'
        p.GetPolar(foil, naca, alphas, Re_w, pane=True)
        print(f'did xfoil on {airfoil_w}')
  
    canard_data = np.loadtxt(f_c, unpack=True, skiprows=12)
    fuselage_data = np.loadtxt(f_f, unpack=True, skiprows=12)
    wing_data = np.loadtxt(f_w, unpack=True, skiprows=12)
    
    
    
    
    
equilibrium(xcg, xc, xf, xw, bc, bw, bf, cc, cf, cw, airfoil_c, airfoil_f, airfoil_w, rho, mu, W, u_inf)