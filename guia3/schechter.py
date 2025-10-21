import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
_params = {
    'figure.figsize':(5,4),
    #'figure.figsize':(6.47,4), # golden ratio
    'figure.dpi':96,
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
    #'axes.axisbelow':True,
}
plt.rcParams.update(_params)

_folder = 'figures/'

def transform_ABC2schechter(logA, B, C):
    """
    logA = log(0.4 ln(10) phi^{*})
    B = M^{*} : magnitud de escala
    C = 1+alpha : pendiente
    ---
    phi_{*} = (10**logA)/(0.4 ln(10))
    M^{*} = B
    alpha = C-1
    """
    phi_star = (10**logA)/(0.4*np.log(10))
    M_star = B
    alpha = C-1.0
    return dict(phi_star=phi_star, M_star=M_star, alpha=alpha)

def transform_errs_ABC2schechter(logA, B, C, e_logA, e_B, e_C):
    """
    phi_{*} = (10**logA)/(0.4 ln(10))
    e_phi = ()/(0.4 ln(10))
    M^{*} = B
    alpha = C-1
    """

    e_phi_star = 2.5*(10**logA)*e_logA
    e_M_star = e_B
    e_alpha = e_C
    return dict(phi_star=e_phi_star, M_star=e_M_star, alpha=e_alpha)

def schechter(x, logA, B, C):
    """
    fit func para log(phi)
    x = M: mag absoulta
    logA = log(0.4 ln(10) phi^{*}) : ...
    B = M^{*} : magnitud de escala
    C = 1+alpha : pendiente
    returns log10(phi) 
    """
    return logA - 0.4*(x-B)*C - 10.0**(-0.4*(x-B))/np.log(10.0)

def fit_lumfunc(lf):
    p0 = [np.log(10.0)*1.5e-2, -21.0, -1.2+1]

    popt, cov = curve_fit(
        schechter, 
        lf['M_r'], 
        np.log10(lf['phi']*np.diff(lf['M_r'])[0]),
        p0=p0
    )
    return popt, cov

def plot_lf(lf, popt, cov):

    bbox_args = dict(boxstyle="round", fc="0.8")
    #params = dict(phi_star=10**(popt[0]), M_star=popt[1], alpha=popt[2]-1)
    params = transform_ABC2schechter(*popt)

    e = np.sqrt(np.diag(cov))
    errs = transform_errs_ABC2schechter(*popt, *e)

    params_text = (
        f'$\phi^{{*}} = ({10**4*params["phi_star"]:2.1f} \pm {10**4*errs["phi_star"]:2.1f}) \\times 10^{{-4}}$\n'
        f'$M^{{*}} = {params["M_star"]:2.2f} \pm {errs["M_star"]:2.2f}$\n'
        f'$\\alpha = {params["alpha"]:2.2f} \pm {errs["alpha"]:2.2f}$'
    )


    fig, ax = plt.subplots()
    ax.plot(lf['M_r'], lf['phi']*np.diff(lf['M_r'])[0], '.', c='k', label='Data')
    ax.plot(lf['M_r'], 10.0**schechter(lf['M_r'], *popt), c='r', label='Schechter')
    ax.semilogy()
    ax.set_xlabel(r'$M_r$')
    ax.set_ylabel(r'$\log(\phi)$')
    ax.legend(loc='upper left')

    ax.annotate(
        text=params_text,
        xy=(-18, 5e-5),
        bbox=bbox_args,
        horizontalalignment='center',
        linespacing=1.5
    )

if __name__ == '__main__':

    lf = pd.read_fwf('lumfunc_n20.dat', names=['ngx', 'phi', 'Mlow', 'Mhigh'])
    lf['M_r'] = 0.5*(lf['Mlow'] + lf['Mhigh'])
    
    popt, cov = fit_lumfunc(lf)
    params = transform_ABC2schechter(*popt)
    errs = np.sqrt(np.diag(cov))
    errs = transform_errs_ABC2schechter(*popt, *errs)

    print(f'{params=}')
    print(f'{errs=}')

    plot_lf(lf, popt, cov)
    plt.show()