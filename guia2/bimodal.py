from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import scienceplots

# plt.style.use('science')
# _params = {
#     'figure.figsize':(5.6,3.5),
#     'figure.dpi':100,
#     'font.family':'sans-serif',
#     'font.size':12,
#     'savefig.format':'png',
#     'savefig.transparent':True,
#     'xtick.direction':'in',
#     'xtick.top':True,
#     'ytick.direction':'in',
#     'ytick.right':True,
#     'errorbar.capsize':3,
#     'legend.loc':'upper right',
#     'legend.frameon':True,
#     'legend.fontsize':9,
# }
# plt.rcParams.update(_params)

_folder = 'figures/'

def normal(x, mu, sigma):
    A = 1.0/(np.sqrt(2.0*np.pi)*sigma)
    return A*np.exp(-0.5*((x-mu)/sigma)**2)

def bimodal(x, w1, mu1, sigma1, w2, mu2, sigma2):
    return w1*normal(x, mu1, sigma1)+w2*normal(x, mu2, sigma2)

def chisq(yobs, yfit):
    return np.sum((yobs-yfit)**2/yfit)

def fit_bimodal(dist, p0, nbins=50, plots=False):

    y, xedges = np.histogram(dist, bins=nbins, density=True)
    x = 0.5*(xedges[:-1]+xedges[1:])
    
    popt, cov = curve_fit(
        bimodal, 
        x, y, 
        p0=p0,
        method='lm',
    )

    yfit = bimodal(x, *popt)
    #print(popt, cov)
    print(f'{popt=}, \n{np.sqrt(np.diag(cov))=}')
    print(f'{chisq(y, yfit)=}')
    print(f'{popt[0]+popt[3]=}')

    if plots:
        plt.figure()
        
        plt.stairs(edges=xedges, values=y, hatch='///', color='dimgray')
        plt.plot(x, yfit, label='Ajuste bimodal', c='k', lw=2)
        plt.plot(x, popt[0]*normal(x,popt[1],popt[2]), label='Comp 1', c='C0', ls='--')
        plt.plot(x, popt[3]*normal(x,popt[4],popt[5]), label='Comp 2', c='C3', ls='--')
        
        plt.xlabel('Dist')
        plt.ylabel('Densidad de objetos')


        plt.legend()
        #plt.show()

    return dict(x=x, xedges=xedges, y=y, yfit=yfit, popt=popt, cov=cov)

if __name__ == '__main__':

    G = pd.read_fwf('gals.dat')
    G.query('u_r < 4.0 and g_r < 1.3 and c9050 > 1.5 and c9050 < 4', inplace=True)

    g_r = fit_bimodal(G['g_r'], p0=[0.6, 0.3, 0.1, 0.4, 0.8, 0.1], plots=True)
    u_r = fit_bimodal(G['u_r'], p0=[0.5, 1.5, 0.3, 0.5, 2.5, 0.2], plots=True)

    plt.show()