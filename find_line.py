import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from astropy.table import Table
from astropy.io import ascii
from astropy.convolution import convolve, Box1DKernel

x1 = []

def on_key(event):
    global x1
    if event.key == 'w':
        x1.append(event.xdata)
        print('%.3f' %event.xdata)
        plt.close()
    if event.key == 'e':
        x1.append(-1)
        print('%.3f' %event.xdata)
        plt.close()


"""
Plots each 1144 spectrum in turn to identify the line position. If no lines, marks it with -1

"""
lines = [8498.02,8542.09,8662.14]
path = '/home/david/work/1144_lt/spectra/nicola_2/WDJ114404.76+052951.77/' #note the / on the end

vispath = path + 'VIS_notell/'
sp = glob.glob(vispath+'*TAC.csv')

for s in sp:
    print(s)
    w, f, e = np.loadtxt(s, unpack=True, delimiter=',')
    mask = (w > 8450) & (w < 8700)
    w, f = w[mask], f[mask]
    f = convolve(f,Box1DKernel(6))
    fig =plt.figure()
    plt.plot(w, f)
    [plt.axvline(line, c='r', ls='--') for line in lines]
    cid = fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()
    
savedat = Table([sp, x1], names=['FILENAME', 'Xs'])
ascii.write(savedat, 'line_positions.ecsv', format='ecsv', overwrite=True)
    