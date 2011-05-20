#!/usr/bin/env python

from __future__ import division

# Natural units are meter, bigmass, year
timeunit = 31556926.

# We write the SI units in terms of these natural units
meter = 1.
second = 1./timeunit
rhog = 1                        # Bigmass is constructed so that rho*grav = 1
massunit = (9.81 * meter * second**(-2) * 910 * meter**(-3)) / rhog
kg = 1/massunit

print('mks:',meter,kg,second)

year = 31556926 * second
Newton = kg * meter / second**2
Pascal = Newton / meter**2

print('Newton:',Newton)
print('Pascal:',Pascal)

n = 3
A = 1.e-16 * Pascal**(-3) / year # EISMINT I number
rho = 910 * kg / meter**3
grav = 9.81 * meter * second**(-2)
B = A**(-1/n)
print('A:',A)
print('B:',B)


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
