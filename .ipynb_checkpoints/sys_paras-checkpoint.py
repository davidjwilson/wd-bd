#!/usr/bin/env python3

import math
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def solveq(q,qf):
    s = q*((1+q)**2)
    return s-qf

# CONSTANTS
MSUN = 1.989E30
RSUN = 6.9599E8
MJUP = 1.898E27
G    = 6.673E-11
C    = 299792458.0

# BINARY PARAMETERS  --- EDIT THESE ---
Porb = 0.04737 # Orbital period (days)
K1   = 26.     # Radial velocity semi-amplitude of WD (km/s)
Ke   = 300.    # Radial velocity semi-amplitude of emission from BD (km/s)
M1   = 0.49    # WD mass (solar units)
R1   = 0.018   # WD radius (solar units)
R2   = 0.08    # BD radius (solar units)
R2a  = 0.13    # R2/a, scaled radius of BD

# Convert period to seconds
Ps = Porb * 86400.

# We will calculate parameters as a function of inclination
inc = np.linspace(1,90)

# Create empty arrays for mass ratio and secondary star mass
Q  = np.zeros(len(inc))
M2 = np.zeros(len(inc))

# Loop through inclination, calculate mass function then solve
# for mass ratio and thus M2
for i in range(len(inc)):
    sini = np.sin(np.radians(inc[i]))
    Mf = (M1*2.*math.pi*G*MSUN*(sini**3))/(((K1*1000.)**3)*Ps)
    q0 = [2.0]
    q  = fsolve(solveq,q0,args=(Mf))
    M2[i] = M1/q
    Q[i]  = 1./q

# Convert M2 to Jupiter masses
M2 *= MSUN/MJUP

# Since Ke is a lower limit on K2, we can set an upper limit on
# the mass ratio when K2=Ke
poss = Q<(K1/Ke)

# Calculate K2 based on uniform emission
f    = 0.5
diff = 1.
K2   = Ke
while math.fabs(diff)>0.0001:
  q    = K1/K2
  K2   = Ke / (1 - (f*(1+q)*R2a))
  diff = q - (K1/K2)

# Calculate M2 based on uniform emission
M2u = M1 * (K1/K2)

# Convert to Jupiter masses
M2uJ = M2u * (MSUN/MJUP)

print("Assuming uniform emission K2 =","%.*f" % (1,K2),"-> Companion mass =","%.*f" % (2,M2uJ),"MJup")

# Calculate all other parameters assuming uniform emission
uinc   = np.degrees(np.arcsin(((((1+q)**2)*Ps*((K2*1000.)**3))/(2.*math.pi*G*M1*MSUN))**(1./3.)))
print("Assuming uniform emission inclination =","%.*f" % (1,uinc))
sini   = np.sin(np.radians(uinc))
q      = K1/K2
vscale = (K1+K2)/sini
a      = 1000.*vscale/(2.*np.pi/Ps)/RSUN
z1     = 0.635*((M1/R1)+(M2u/a)) + (K1/sini)**2/(2.*C)
z2     = 0.635*((M2u/R2)+(M1/a)) + (K2/sini)**2/(2.*C)
z      = z1-z2

print("Implied white dwarf gravitational redshift =","%.*f" % (1,z),"km/s")
print("Inclination: >","%.*f" % (1,inc[poss].min()),"degrees")
print("Companion mass:","%.*f" % (2,M2[-1]),"->","%.*f" % (2,M2[poss].max()),"MJup")

plt.plot(inc,M2,'--')
plt.plot(inc[poss],M2[poss],'r')
plt.plot([inc[poss].min(),inc[poss].min()], [0,1000],'k--')
plt.plot([uinc,uinc], [0,1000],'g--')
plt.xlabel("Inclination (degrees)")
plt.ylabel("M2 (Jupiter masses)")
plt.xlim(0,90)
plt.ylim(0,200)
plt.show()
