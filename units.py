#!/usr/bin/env python

# Natural units are meter, bigmass, year
# We write the SI units in terms of these natural units
year = 1
meter = 1.
second = 1./31556926.
rhog = 1                        # Bigmass is constructed so that rho*grav = 1
kg = rhog / (9.81 * meter * second**(-2) * 910 * meter**(-3))

print('mks:',meter,kg,second)

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
