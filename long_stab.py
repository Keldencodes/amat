import matplotlib.pyplot as plt
import numpy as np
import conversions as c
import os
import pyxfoil as p

# define some constants
# spans in inches
bc = 36 # canard
bw = 96 # wing
bf = 24 # fuselage with 2 columns of soccer balls

# chords in inches
cc = 8
cf = 10/.22 # using 1 row of soccer balls to define fuselage size
#cf2 = 10/.22
cw = 23

# distances in inches to xac and xcg
xc = cc/4
xf = cf/4
#xf2 = cf2/4
xw = cf - 3/4*cw # place trailing edge of wing coincident with TE of fuselage
xcg = np.linspace(30, 15, 41) # center of gravity location from leading edge of fuselage
# should probably define the cg later using the xac=xcg large xbarac function
u_inf = 50 # mph

rho = 1.225 # kg/m3
mu = 1.802e-5 # kg/m*s

w = 40 # lbs

airfoil_c = 'naca0012'
airfoil_f = '2422'
#airfoil_f2 = '2422'
airfoil_w = 'naca2415'

def equilibrium(xcg, xc, xf, xw, bc, bw, bf, cc, cf, cw, airfoil_c, airfoil_f, airfoil_w, rho, mu, w, u_inf):
    # calculate required alpha of canard with varying fuselage/wing alpha
    # conversions
    cc = c.in2m(cc)
    cf = c.in2m(cf)
    #cf2 = c.in2m(cf2)
    cw = c.in2m(cw)
    
    bc = c.in2m(bc)
    bf = c.in2m(bf)
    bw = c.in2m(bw)
    
    xcg = c.in2m(xcg)
    xc = c.in2m(xc)
    xf = c.in2m(xf)
    xw = c.in2m(xw)
    
    u_inf = c.mph2mpers(u_inf)
    w = c.pound2n(w)
    
    # surface areas in m2
    Sc = cc * bc
    Sf = cf * bf
    #Sf2 = cf2 * bf
    Sw = cw * bw
    
    # reynolds of each 'wing'
    Re_c = rho * u_inf * cc / mu
    Re_f = rho * u_inf * cf / mu
    #Re_f2 = rho * u_inf * cf2 / mu
    Re_w = rho * u_inf * cw / mu
    
    airfoil_c_path = os.path.join('data', f'{airfoil_c}.dat')
    with open(airfoil_c_path, 'r') as infile:
        x_c, y_c = np.loadtxt(infile, unpack=True, skiprows=1)
        
    # airfoil_f_path = os.path.join('data', f'{airfoil_f}.dat')
    # with open(airfoil_f_path, 'r') as infile:
    #     x_f, y_f = np.loadtxt(infile, unpack=True, skiprows=1)
    
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

    
    f_c = f'data/{airfoil_c}/{airfoil_c}_polar_Re{Re_c:.2e}a{start_c:.1f}-{stop_c:.1f}.dat'
    f_f = f'data/naca{airfoil_f}/naca{airfoil_f}_polar_Re{Re_f:.2e}a{start_f:.1f}-{stop_f:.1f}.dat'
    #f_f2 = f'data/naca{airfoil_f2}/naca{airfoil_f2}_polar_Re{Re_f2:.2e}a{start_f:.1f}-{stop_f:.1f}.dat'
    f_w = f'data/{airfoil_w}/{airfoil_w}_polar_Re{Re_w:.2e}a{start_w:.1f}-{stop_w:.1f}.dat'

    # if the xfoil data already exists, don't run it again
    if not os.path.exists(f_c):
        naca = False
        foil = f'data/{airfoil_c}.dat'
        p.GetPolar(foil, naca, alphas_c, Re_c, pane=True)
        print('did xfoil')

    if not os.path.exists(f_f):
        naca = True
        foil = airfoil_f
        p.GetPolar(foil, naca, alphas_f, Re_f, pane=True)
        print('did xfoil')

    # if not os.path.exists(f_f2):
    #     naca = True
    #     foil = airfoil_f2
    #     p.GetPolar(foil, naca, alphas_f, Re_f2, pane=True)
    #     print('did xfoil')
    
    if not os.path.exists(f_w):
        naca = False
        foil = f'data/{airfoil_w}.dat'
        p.GetPolar(foil, naca, alphas_w, Re_w, pane=True)
        print('did xfoil')
  
    canard_data = np.loadtxt(f_c, unpack=True, skiprows=12)
    fuselage_data = np.loadtxt(f_f, unpack=True, skiprows=12)
    #fuselage2_data = np.loadtxt(f_f2, unpack=True, skiprows=12)
    wing_data = np.loadtxt(f_w, unpack=True, skiprows=12)

    # plot cl vs a for each
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(canard_data[0], canard_data[1], 'r', label=f'Canard: {airfoil_c}')
    ax1.plot(fuselage_data[0], fuselage_data[1], 'g', label=f'Fuselage: {airfoil_f}')
    #ax1.plot(fuselage2_data[0], fuselage2_data[1], 'm', label=f'Fuselage2: {airfoil_f2}')
    ax1.plot(wing_data[0], wing_data[1], 'b', label=f'Wing: {airfoil_w}')
    ax1.legend()
    ax1.set_xlabel('alpha')
    ax1.set_ylabel('Cl')

    # plot cd vs a for each
    ax2.plot(canard_data[0], canard_data[2], 'r', label=f'Canard: {airfoil_c}')
    ax2.plot(fuselage_data[0], fuselage_data[2], 'g', label=f'Fuselage: {airfoil_f}')
    #ax2.plot(fuselage2_data[0], fuselage2_data[2], 'm', label=f'Fuselage2: {airfoil_f2}')
    ax2.plot(wing_data[0], wing_data[2], 'b', label=f'Wing: {airfoil_w}')
    ax2.legend()
    ax2.set_xlabel('alpha')
    ax2.set_ylabel('Cd')
    plt.show()

    alpha_c = 2 #canard_data[0][np.where(canard_data[1]==max(canard_data[1]))[0][0]] # stall angle
    clc, cdc = find_cl_cd(alpha_c, canard_data)
    cdc = cdc + (clc**2)/(np.pi*0.9*(bc/cc))

    alphaf = fuselage_data[0][np.where(fuselage_data[1]==max(fuselage_data[1]))[0][0]] # stall angle
    alphaw = wing_data[0][np.where(wing_data[1]==max(wing_data[1]))[0][0]] # stall angle
    alpha = 0.5 #min(alphaf, alphaw) # smaller of the wing and fuselage stall angles
    clf, cdf = find_cl_cd(alpha, fuselage_data)
    cdf = cdf + (clf**2)/(np.pi*0.9*(bf/cf))
    #clf2, cdf2 = find_cl_cd(alpha, fuselage2_data)

    clw, cdw = find_cl_cd(alpha, wing_data)
    cdw = cdw + (clw**2)/(np.pi*0.9*(bw/cw))
    
    Lc = .5 * clc * rho * (u_inf**2) * Sc
    Lw = .5 * clw * rho * (u_inf**2) * Sw
    Lf = .5 * clf * rho * (u_inf**2) * Sf
    #Lf2 = .5 * clf2 * rho * (u_inf**2) * Sf

    Dc = .5 * cdc * rho * (u_inf**2) * Sc
    Dw = .5 * cdw * rho * (u_inf**2) * Sw
    Df = .5 * cdf * rho * (u_inf**2) * Sf
    #Df2 = .5 * cdf2 * rho * (u_inf**2) * Sf

    if type(xcg)==np.ndarray:
        all_done = False
        for i in xcg:
            l = (rho*(u_inf**2))/(2*(i-xc)) * (clw*Sw*(xw-xc) + clf*Sf*(xf-xc)) # lift of wing and fuselage about x_canard
            if l-w >= 0: # want the wing and fuselage pitch about the x_canard to outweigh mass pitch about x_canard
                total_lift = Lc+Lw+Lf
                if total_lift >= w:
                    l_c = (Lw*(xw-i)+Lf*(xf-i))/(i-xc) # lift of wing required in moment eq about cg
                    if Lc <= l_c: # want the canard to be stalled before the others
                        # print(f'canard moment lift required: {l_c:.3f}N')
                        # print(f'max L_c: {Lc:.3f}N')
                        xcg_i = i
                        all_done = True
                        break
        if all_done == False:
            xcg_i = i
    else:
        l = (rho*(u_inf**2))/(2*(xcg-xc)) * (clw*Sw*(xw-xc) + clf*Sf*(xf-xc))
    
    print('**********************************************')
    print('Results')
    print(f'Aircraft mass: {c.n2pound(w)} lbs')
    print(f'Pitching Lift: {(c.n2pound(l-w)):.2f}')
    print(f'Total lift: {(Lc+Lw+Lf):.2f}N, {c.n2pound(Lc+Lw+Lf):.2f}lbs')
    print(f'Total drag = {(Dc+Dw+Df):.2f}N, {c.n2pound(Dc+Dw+Df):.2f}lbs')
    print(f'Relative airspeed = {c.mpers2mph(u_inf):.1f}mph')
    print(f'Xcg = {c.m2in(xcg_i):.2f}in')
    print(f'Fuselage size \n Chord = {(c.m2in(cf)):.1f}in \n Span = {(c.m2in(bf)):.1f}in \n Area = {(c.m2in(cf)*c.m2in(bf)):.2f}in2')
    #print(f'Fuselage2 size \n Chord = {(c.m2in(cf2)):.1f}in \n Span = {(c.m2in(bf)):.1f}in \n Area = {(c.m2in(cf2)*c.m2in(bf)):.2f}in2')
    print(f'Wing size \n Chord = {(c.m2in(cw)):.1f}in \n Span = {(c.m2in(bw)):.1f}in \n Area = {(c.m2in(cw)*c.m2in(bw)):.2f}in2')
    print(f'Canard size \n Chord = {(c.m2in(cc)):.1f}in \n Span = {(c.m2in(bc)):.1f}in \n Area = {(c.m2in(cc)*c.m2in(bc)):.2f}in2')
    print(f'Fuselage \n Lift: Cl = {clf:.3f}, L = {Lf:.2f}N = {c.n2pound(Lf):.2f}lbs \n Drag: Cd = {cdf:.3f}, D = {Df:.2f}N = {c.n2pound(Df):.2f}lbs \
        \n alpha: {alpha} deg \n crit a: {fuselage_data[0][np.where(fuselage_data[1]==max(fuselage_data[1]))[0][0]]} \n Re: {Re_f:.2f}')
    print(f'Wing \n Lift: CL = {clw:.3f}, L = {Lw:.2f}N = {c.n2pound(Lw):.2f}lbs \n Drag: Cd = {cdw:.3f}, D = {Dw:.2f}N = {c.n2pound(Dw):.2f}lbs \
        \n alpha: {alpha} deg \n crit a: {wing_data[0][np.where(wing_data[1]==max(wing_data[1]))[0][0]]} \n Re: {Re_w:.2f}')
    print(f'Canard \n Lift: Cl = {clc:.3f}, L = {Lc:.2f}N = {c.n2pound(Lc):.2f}lbs \n Drag: Cd = {cdc:.3f}, D = {Dc:.2f}N = {c.n2pound(Dc):.2f}lbs \
        \n alpha: {alpha_c} deg \n crit a: {canard_data[0][np.where(canard_data[1]==max(canard_data[1]))[0][0]]} \n Re: {Re_c:.2f}')
    
    # print(f'Fuselage2 lift: Cl = {clf2:.3f}, L = {Lf2:.2f}N = {c.n2pound(Lf2):.2f}lbs')
    # print(f'Fuselage2 drag: Cd = {cdf2:.3f}, L = {Df2:.2f}N = {c.n2pound(Df2):.2f}lbs')
    # print(f'Fuselage2 Re = {Re_f2:.2f}')
    
def find_cl_cd(alpha, data):
    result = np.where(data[0]==alpha)
    if result[0].size==0:
       alpha1_index = max(np.where(data[0]<alpha)[0])
       alpha2_index = alpha1_index + 1
       alpha1 = data[0][alpha1_index]
       alpha2 = data[0][alpha2_index]
       cl1 = data[1][alpha1_index]
       cl2 = data[1][alpha2_index]
       cl = c.interpolate(alpha1, cl1, alpha2, cl2, alpha)
       cd1 = data[2][alpha1_index]
       cd2 = data[2][alpha2_index]
       cd = c.interpolate(x1=alpha1, y1=cd1, x2=alpha2, y2=cd2, x3=alpha)
    else:
       cd = data[2][result[0][0]]
       cl = data[1][result[0][0]]
    return cl, cd    
    
equilibrium(xcg, xc, xf, xw, bc, bw, bf, cc, cf, cw, airfoil_c, airfoil_f, airfoil_w, rho, mu, w, u_inf)