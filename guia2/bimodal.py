from scipy.optimize import curve_fit, minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
_params = {
    'figure.figsize':(5.6,3.5),
    'figure.dpi':100,
    'font.family':'sans-serif',
    'font.size':12,
    'savefig.format':'png',
    'savefig.transparent':True,
    'xtick.direction':'in',
    'xtick.top':True,
    'ytick.direction':'in',
    'ytick.right':True,
    'errorbar.capsize':3,
    'legend.loc':'upper right',
    'legend.frameon':True,
    'legend.fontsize':9,
}
plt.rcParams.update(_params)

_folder = 'figures/'

def normal(x, mu, sigma):
    A = 1.0/(np.sqrt(2.0*np.pi)*sigma)
    return A*np.exp(-0.5*((x-mu)/sigma)**2)

def bimodal(x, w1, mu1, sigma1, w2, mu2, sigma2):
    return w1*normal(x, mu1, sigma1)+w2*normal(x, mu2, sigma2)

def chisq(yobs, yfit):
    return np.sum((yobs-yfit)**2/yfit)

def fit_colors(color='u_r', nbins=50, plots=True):

    G = pd.read_fwf('gals.dat')
    G.query('u_r < 4.0 and g_r < 1.3 and c9050 > 1.5 and c9050 < 4', inplace=True)

    y, xedges = np.histogram(G[color], bins=nbins, density=True)
    x = 0.5*(xedges[:-1]+xedges[1:])

    if color == 'u_r':
        p0=[0.5, 1.5, 0.3, 0.5, 2.5, 0.2]
    else:
        p0=[0.6, 0.3, 0.1, 0.4, 0.8, 0.1]
    
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
        
        plt.stairs(edges=xedges, values=y, label='SDSS', hatch='///', color='dimgray')
        plt.plot(x, yfit, label='Ajuste bimodal', c='k', lw=2)
        plt.plot(x, popt[0]*normal(x,popt[1],popt[2]), label='Nube azul', c='C0', ls='--')
        plt.plot(x, popt[3]*normal(x,popt[4],popt[5]), label='Secuencia roja', c='C3', ls='--')
        
        plt.xlabel('$u-r$' if color=='u-r' else '$g-r$')
        plt.ylabel('Densidad de objetos')

        plt.legend(ncols=2)
        plt.show()

    return dict(x=x, xedges=xedges, y=y, yfit=yfit, popt=popt, cov=cov)

if __name__ == '__main__':

    g_r = fit_colors('g_r', plots=False)
    u_r = fit_colors('u_r', plots=False)

    np.savetxt(
        'colors_fit.dat', 
        np.vstack([
            u_r['popt'], 
            np.sqrt(np.diag(u_r['cov'])),
            g_r['popt'], 
            np.sqrt(np.diag(g_r['cov'])),
        ]).T,
        fmt='%8.6G', 
        header='u-r_popt u-r_perr g-r_popt g-r_perr',
        comments='#w1 mu1 sigma1 w2 mu2 sigma2 \n'
    )

    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(9,4))

    ax1.stairs(
        edges=g_r['xedges'],
        values=g_r['y'],
        label='SDSS',
        hatch='///',
        color='dimgray',
    )
    ax1.plot(
        g_r['x'],
        g_r['yfit'],
        label='Ajuste bimodal', 
        c='k', lw=1.8, alpha=0.8
    )
    ax1.plot(
        g_r['x'],
        g_r['popt'][0]*normal(g_r['x'], *g_r['popt'][1:3]),
        label='Nube azul', 
        c='b', ls='--'
    )
    ax1.plot(
        g_r['x'],
        g_r['popt'][3]*normal(g_r['x'], *g_r['popt'][4:]),
        label='Sec. roja', 
        c='r', ls='--'
    )

    ax1.set_xlabel('$g-r$')
    ax1.set_ylabel('Densidad de galaxias')
    #ax1.legend(ncols=2, loc='upper left')

    ax2.stairs(
        edges=u_r['xedges'],
        values=u_r['y'],
        label='SDSS',
        hatch='///',
        color='dimgray',
    )
    ax2.plot(
        u_r['x'],
        u_r['yfit'],
        label='Ajuste bimodal', 
        c='k', lw=1.8, alpha=0.8
    )
    ax2.plot(
        u_r['x'],
        u_r['popt'][0]*normal(u_r['x'], *u_r['popt'][1:3]),
        label='Nube azul', 
        c='b', ls='--'
    )
    ax2.plot(
        u_r['x'],
        u_r['popt'][3]*normal(u_r['x'], *u_r['popt'][4:]),
        label='Sec. roja', 
        c='r', ls='--'
    )

    ax2.set_xlabel('$u-r$')
    #ax2.set_ylabel('Densidad de galaxias')
    ax2.legend(ncols=2)
    
    fig.savefig(_folder+'colors_fit.png')
    #plt.show()