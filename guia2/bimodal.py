from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

_folder = 'figures/'

def normal(x, mu, sigma):
    A = 1.0/(np.sqrt(2.0*np.pi)*sigma)
    return A*np.exp(-0.5*((x-mu)/sigma)**2)

def bimodal(x, w1, mu1, sigma1, w2, mu2, sigma2):
    return w1*normal(x, mu1, sigma1)+w2*normal(x, mu2, sigma2)

def chisq(yobs, yfit):
    return np.sum((yobs-yfit)**2/yfit)

def fit_bimodal(dist, p0, weights=None, nbins=50, plots=False, verbose=False):

    y, xedges = np.histogram(dist, bins=nbins, density=True, weights=weights)
    x = 0.5*(xedges[:-1]+xedges[1:])
    
    popt, cov = curve_fit(
        bimodal, 
        x, y, 
        p0=p0,
        method='lm',
    )

    yfit = bimodal(x, *popt)
    #print(popt, cov)
    if verbose:
        print(f'{popt=}, \n{np.sqrt(np.diag(cov))=}')
        print(f'{chisq(y, yfit)=}')
        print(f'{popt[0]+popt[3]=}')

    if plots:
        fig, ax = plt.subplots()
        
        ax.stairs(edges=xedges, values=y, hatch='///', color='dimgray')
        ax.plot(x, yfit, label='Ajuste bimodal', c='k', lw=2)
        ax.plot(x, popt[0]*normal(x,popt[1],popt[2]), label='Comp 1', c='C0', ls='--')
        ax.plot(x, popt[3]*normal(x,popt[4],popt[5]), label='Comp 2', c='C3', ls='--')
        
        ax.set_xlabel('Dist')
        ax.set_ylabel('Densidad de objetos')

        ax.legend()
        #plt.show()
    else:
        ax = None

    return dict(x=x, xedges=xedges, y=y, yfit=yfit, popt=popt, cov=cov, ax=ax)

if __name__ == '__main__':

    G = pd.read_fwf('gals.dat')
    G.query('u_r < 4.0 and g_r < 1.3 and c9050 > 1.5 and c9050 < 4', inplace=True)

    g_r = fit_bimodal(G['g_r'], p0=[0.6, 0.3, 0.1, 0.4, 0.8, 0.1], plots=True)
    u_r = fit_bimodal(G['u_r'], p0=[0.5, 1.5, 0.3, 0.5, 2.5, 0.2], plots=True)

    plt.show()