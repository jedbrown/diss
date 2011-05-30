#!/usr/bin/env python

from __future__ import division

from collections import namedtuple

SIBase = namedtuple('SIBase', 'm kg s')

def printlines(args):
    print('\n'.join(args))
def reciprocal(*lst):
    return [1/x for x in lst]
def baseunits_si():
    return 1,1,1
def baseunits_natural():
    # Natural units are meter, bigmass, year
    timeunit = 31556926.
    # We write the SI units in terms of these natural units
    meter = 1.
    second = 1./timeunit
    rhog = 1                        # Bigmass is constructed so that rho*grav = 1
    massunit = (9.81 * meter * second**(-2) * 910 * meter**(-3)) / rhog
    return reciprocal(1, massunit, timeunit)
def baseunits_natural2():
    length = 100e3
    mass = 1e6
    time = 0.001
    return reciprocal(length, mass, time), 1e3
def baseunits_natural3():
    length = 1
    mass = 1e-8
    time = 1e-7
    return reciprocal(length, mass, time), 1e-12
def baseunits_natural4():
    length = 1
    # we want rho*grav*h = 10^7
    mass = 1e-9
    time = 1e-8
    return reciprocal(length, mass, time), 1e-15 # volume = 1e13

base, basescale = baseunits_natural4()
meter, kg, second = base
print('mks: %g %g %g ' % (meter,kg,second))
print ('rho g h = %g    rho u = %g' % (1e7*kg*meter**(-1)*second**(-2), 1e-1*kg*meter**(-2)*second**(-1)))

year = 31556926 * second
Newton = kg * meter / second**2
Pascal = Newton / meter**2
Joule = Newton * meter
Watt = Joule / second
Kelvin = 1

print('Newton:',Newton)
print('Pascal:',Pascal)

n = 3
A = 1.e-16 * Pascal**(-3) / year # EISMINT I number
rho = 910 * kg / meter**3
grav = 9.81 * meter * second**(-2)
B = A**(-1/n)
print('A:',A)
print('B:',B)

if True:
    c_i = 2009 * Joule / (kg*Kelvin)
    k_T = 2.1 * Watt / (meter*Kelvin)
    Latent = 3.34e5 * Joule / kg

if True:
    velocity = 10e3*meter/year
    momentum = rho*velocity
    pressure = rho*grav*1e3*meter
    energy = c_i * (1*Kelvin) * rho
    etasi = 1e13
    viscosity = etasi * Pascal*second
    strainrate = 1/year
    stress = viscosity * strainrate
    printlines(['Field values:',
                '  Density: %g' % (rho,),
                '  Velocity 10km/year: %g ~~ m s^-1' % (velocity,),
                '  *Momentum density at 10km/year: %g ~~ kg m^-2 s^-1' % (momentum,),
                '  *Hydrostatic pressure at 1km depth: %g ~~ kg m^-1 s^-2' % (pressure,),
                '  *Energy density change for 10K: %g ~~ kg m^-1 s^-2' % (energy,),
                '  Viscosity %g Pa s: %g' % (etasi,viscosity),
                '  Momentum viscosity %g Pa s: %g' % (etasi,viscosity/rho),
                '  Stress at 1/year: %g' % stress,
                ])
    for xlen,zlen,scale in ((1e3,1e3,1), (1e5,1e3,1e3), (1e5,1e3,basescale)):
        xlen *= meter; zlen *= meter
        volume = xlen**2*zlen
        vscale = scale * volume
        printlines([('Integrated over %.0g x %.0g x %.0g m^3 = (%g Length)^3 = %g Length^3, scale=%g'
                     % (xlen/meter,xlen/meter,zlen/meter,volume**(1/3),volume,scale)),
                    '  Mass: %g' % (rho*vscale,),
                    '  *Momentum at 10km/year: %g ~~ kg m s^-1' % (momentum*vscale),
                    '  *Pressure: %g ~~ kg m^2 s^-2' % (pressure*vscale),
                    '  *Energy in 10K change: %g ~~ kg m^2 s^-2' % (energy * vscale,),
                    '  Gravitational force: %g ~~ kg m s^-2' % (rho * grav * vscale),
                    '  Stress at 1/year, eta=%g Pa s: %g ~~ kg m^2 s^-2' % (etasi,stress*vscale),
                    '  Scaled momentum viscosity: %g' % (stress / (rho*strainrate) * vscale),
                    '  Scaled energy diffusivity: %g' % (vscale*1e-6*meter**2/second),
                    '  Scaled viscous heat production: %g' % (vscale*stress*strainrate),
                ])

def viscosity(Du,p,T):
    from numpy import exp
    R = 8.31441 # {\joule\per\mol\per\kelvin}
    Q = 6.0e4 # {\joule\per\mol} & Activation energy for creep
    V = -13.e-6 # {\metre\cubed\per\mol} & Activation volume for creep
    A_0 = 3.61e-13 # {\per\pascal\cubed\per\second} & Softness parameter
    A = A_0 * exp((-Q + p*V)/(R*T))
    n = 3
    gamma = 0.5*Du**2
    return A**(-1/n) * gamma**((1-n)/(2*n))
