import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt

lf = pd.read_fwf('lumfunc_n20.dat', names=['ngx', 'phi', 'Mlow', 'Mhigh'])
lf['M_r'] = 0.5*(lf['Mlow'] + lf['Mhigh'])

def schechter(x, logA, B, C):
    """
    fit func para log(phi)
    x = M: mag absoulta
    logA = phi^{*} ln(10) : ...
    B = M^{*} : magnitud de escala
    C = 1+alpha : pendiente
    returns log10(phi) 
    """
    return logA - 0.4*(x-B)*C - 10.0**(-0.4*(x-B))/np.log(10.0)

def fit_lumfunc():

    p0 = [np.log(10.0)*1.5e-2, -21.0, -1.2+1]
    popt, cov = curve_fit(schechter, lf['M_r'], np.log10(lf.phi*np.diff(lf['M_r'])[0]), p0=p0)
    print(f'{popt=}')

    return popt, cov

def plot_lf():
    popt, cov = fit_lumfunc()
    params = dict(phi_star = 10**(popt[0]), M_star = popt[1], alpha = popt[2]-1)
    print(f'{params=}')
    fig, ax = plt.subplots()
    c = ax.imshow(np.log10(np.abs(cov)))
    fig.colorbar(c)

    fig, ax2 = plt.subplots()
    ax2.plot(lf['M_r'], lf.phi*np.diff(lf['M_r'])[0], '.', c='k')
    ax2.plot(lf['M_r'], 10.0**schechter(lf['M_r'], *popt), c='r')
    ax2.semilogy()

if __name__ == '__main__':
    plot_lf()
    plt.show()