#!/usr/bin/env 

"""
Find q and K2 with an mc.
"""

#first get the python modules we need
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const


def k2_calc(Kem, q, R2, a, f=0.5):
    """
    Converts the measured K velocity of an emmision line source to a K2 value
    """
    K2 = Kem / (1 - (f*(1+q)*(R2/a)))
    return K2

def mass_radius(mass, type='bd', error=False):
    """
    Mass radius relation from https://ui.adsabs.harvard.edu/abs/2017ApJ...834...17C/abstract
    """
    if type == 'bd':
        radius = mass ** -0.044
        edown = mass ** -(0.044 - 0.017)
        eup = mass ** -(0.044 + 0.019)
        error_array = np.sort([abs(radius -edown), abs(eup -radius)])
        
        # - 0.017, + 0.019
    if type == 'star':
        radius = mass ** 0.881
        edown = mass ** (0.881 - 0.024)
        eup = mass ** (0.881 + 0.025)
        error_array = np.sort([abs(radius -edown), abs(eup -radius)])
        # - 0.024, + 0.025
    if error:
        return radius, error_array
    else:
        return radius

def find_a(m1, m2, p):
    a = ((p.to(u.s))**2 * (const.G * (m1 +m2).to(u.kg))/(4*np.pi**2))**(1/3)
    return a


nsteps= 100

qs = []
k2s = []


kem, keme = 198.09305398172296, 0.4535365132628666
M1, M1e = 0.475735, 0.045966
k1, k1e = 32.86031937992142, 3.9970660761849968
period = 4.223


kems = np.random.normal(kem, keme, nsteps)
M1s = np.random.normal(M1, M1e, nsteps)
k1s = np.random.normal(k1, k1e, nsteps)

for ki in kems:
    for M1i in M1s:
        for k1i in k1s:
            q = k1i/ki
            a = find_a(M1i*u.Msun, (M1i*q)*u.Msun, period*u.h)
            r2 = mass_radius(((M1i*q)*u.Msun).to(u.Mjup).value, error=False)
            r2 = (r2*u.Rjup).to(u.Rsun).value
            diff = 1000
            n = 0
            while diff > 0.0001 and n < 100:
                n +=1
                k2 = k2_calc(ki, q, r2*u.Rsun, a.to(u.Rsun), f=0.5)
                nq = k1i/k2 
                diff = abs(nq -q)
                q = nq
                a = find_a(M1i*u.Msun, (M1i*q)*u.Msun, period*u.h)
                r2 = mass_radius(((M1i*q)*u.Msun).to(u.Mjup).value, error=False)
                r2 = (r2*u.Rjup).to(u.Rsun).value
            # print(n)
            qs.append(q)
            k2s.append(k2)



fig, ax = plt.subplots(ncols=2)                
fig.suptitle('Brown Dwarf')
ax[0].hist(qs)
ax[0].set_xlabel('q')
ax[0].set_ylabel('n')
ax[1].hist(k2s)
ax[1].set_xlabel('k2')
# ax[1].axvline(kem, ls='--', c='C2')
# plt.show()

print('Brown Dwarf')
qmean = np.mean(qs)
qmed = np.median(qs)
qstd = np.std(qs)
print(qmean, qmed, qstd)

k2mean = np.mean(k2s)
k2med = np.median(k2s)
k2std = np.std(k2s)
print(k2mean, k2med, k2std)




qs = []
k2s = []


for ki in kems:
    for M1i in M1s:
        for k1i in k1s:
            q = k1i/ki
            a = find_a(M1i*u.Msun, (M1i*q)*u.Msun, period*u.h)
            r2 = mass_radius(M1i*q, error=False, type='star')
            # r2 = (r2*u.Rjup).to(u.Rsun).value
            diff = 1000
            n = 0
            while diff > 0.0001 and n < 100:
                n +=1
                k2 = k2_calc(ki, q, r2*u.Rsun, a.to(u.Rsun), f=0.5)
                nq = k1i/k2 
                diff = abs(nq -q)
                q = nq
                a = find_a(M1i*u.Msun, (M1i*q)*u.Msun, period*u.h)
                r2 = mass_radius(M1i*q, error=False, type='star')
            # print(n)
            qs.append(q)
            k2s.append(k2)



fig, ax = plt.subplots(ncols=2)                
fig.suptitle('Star')
ax[0].hist(qs)
ax[0].set_xlabel('q')
ax[0].set_ylabel('n')
ax[1].hist(k2s)
ax[1].set_xlabel('k2')
# ax[1].axvline(kem, ls='--', c='C2')
# plt.show()

print('Star')
qmean = np.mean(qs)
qmed = np.median(qs)
qstd = np.std(qs)
print(qmean, qmed, qstd)

k2mean = np.mean(k2s)
k2med = np.median(k2s)
k2std = np.std(k2s)
print(k2mean, k2med, k2std)

plt.show()

print('Done')
