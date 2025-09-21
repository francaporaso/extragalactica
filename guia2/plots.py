import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from bimodal import *
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
_params = {
    'figure.figsize':(5.6,3.5),
    'figure.dpi':100,
    'font.family':'sans-serif',
    'font.size':11,
    'savefig.format':'pdf',
    'savefig.transparent':True,
    'xtick.direction':'in',
    'xtick.top':True,
    'ytick.direction':'in',
    'ytick.right':True,
    'errorbar.capsize':3,
    'legend.frameon':True,
}
plt.rcParams.update(_params)

_folder = 'figures/'

def plot_problema2():

    df = pd.read_csv('SDSS_guia2_fmcaporaso.csv')

    # sky footprint
    fig, ax = plt.subplots(1,1)
    ax.scatter(df.ra, df.dec, s=3, alpha=0.8)
    ax.set_xlabel('Ascención recta [deg]')
    ax.set_ylabel('Declinación [deg]')
    fig.savefig('sky_footprint.pdf')
    #plt.show()

    # redshift distribution
    fig, ax = plt.subplots(1,1)
    ax.hist(df.redshift, bins=25, histtype='step')
    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel('Número de galaxias')
    fig.savefig('z_dist_25bins.pdf')
    #plt.show()

    # mag vs z
    fig, ax = plt.subplots(1,1)
    ax.scatter(df.petroMag_r, df.redshift, s=2, alpha=0.5, c='C0')
    ax.scatter([],[],c='C0',s=8,label='Full Sample')
    mask = (df.petroMag_r >= 14.5)&(df.petroMag_r <= 17.77)
    ax.scatter(df.petroMag_r[mask], df.redshift[mask], s=4, alpha=0.5, c='C1')
    ax.scatter([],[],c='C1',s=8,label='$14.5 \leq r \leq 17.77$')
    ax.set_xlabel('Magnitud petrosiana banda $r$')
    ax.set_ylabel('Redshift $z$')
    ax.legend(loc='upper right', frameon=True)
    fig.savefig('r_vs_redshift.pdf')
    #plt.show()

def plot_problema4():

    df = pd.read_fwf('gals.dat')

    # ========== Mr vs z
    # fig, ax = plt.subplots()
    # ax.scatter(df.z, df.M_pet_r, s=1, alpha=0.5)
    # ax.invert_yaxis()
    # ax.set_xlabel('Redshift $z$')
    # ax.set_ylabel('$M_r$ petrosiana')
    # #plt.show()
    # fig.savefig('mag_vs_z.pdf')
    
    # ========== dist colores
    # fig, (ax1,ax2) = plt.subplots(1,2, sharey=False, figsize=(9,4))
    # ax1.hist(df.u_r, 
    #          bins=np.linspace(0.0,4.0,50),
    #          histtype='step', hatch='///', density=True)
    # ax2.hist(df.g_r, 
    #          bins=np.linspace(0.0,1.3,50),
    #          histtype='step', hatch='///', density=True, color='C1')
    
    # ax1.set_xlabel('$u-r$')
    # ax2.set_xlabel('$g-r$')
    # ax1.set_ylabel('Densidad de objetos')
    # plt.show()

    # ========== Dist c9050 y fracdev_r
    # fig, (ax1,ax2) = plt.subplots(1,2, sharey=False, figsize=(9,4))
    # ax1.hist(df.c9050, 
    #          bins=np.linspace(1.4,4.0,50),
    #          #bins=100,
    #          histtype='step', hatch='///', density=True)
    # ax2.hist(df.fracDeV_r, 
    #          #bins=np.linspace(0.0,1.3,50),
    #          bins=50,
    #          histtype='step', hatch='///', density=True, color='C1')
    
    # ax1.set_xlabel('Indice de concentracion')
    # ax2.set_xlabel('fracDeV_r')
    # ax1.set_ylabel('Densidad de objetos')
    # plt.show()

    # ========== color magnitud
    #fig, ax = plt.subplots()
    #ax.scatter(df.M_pet_r, df.u_r, s=2, alpha=0.5)
    #plt.show()

    ## ajuste de doble gausiana -> probar con equal bin y equal number (quartiles)

def plot_bimodalcolors(save=False):
    G = pd.read_fwf('gals.dat')
    G.query('u_r < 4.0 and g_r < 1.3 and c9050 > 1.5 and c9050 < 4', inplace=True)

    g_r = fit_bimodal(G['g_r'], p0=[0.6, 0.3, 0.1, 0.4, 0.8, 0.1])
    u_r = fit_bimodal(G['u_r'], p0=[0.5, 1.5, 0.3, 0.5, 2.5, 0.2])

    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(9,4))

    ax1.stairs(edges=g_r['xedges'], values=g_r['y'],
               label='SDSS', hatch='///', color='dimgray')
    ax1.plot(g_r['x'], g_r['yfit'],
             label='Ajuste bimodal', c='k', lw=1.8, alpha=0.8)
    ax1.plot(g_r['x'], g_r['popt'][0]*normal(g_r['x'], *g_r['popt'][1:3]),
             label='Nube azul', c='b', ls='--')
    ax1.plot(g_r['x'], g_r['popt'][3]*normal(g_r['x'], *g_r['popt'][4:]),
             label='Sec. roja', c='r', ls='--')
    ax1.set_xlabel('$g-r$')
    ax1.set_ylabel('Densidad de galaxias')

    ax2.stairs(edges=u_r['xedges'], values=u_r['y'],
               label='SDSS', hatch='///', color='dimgray')
    ax2.plot(u_r['x'], u_r['yfit'],
             label='Ajuste bimodal', c='k', lw=1.8, alpha=0.8)
    ax2.plot(u_r['x'], u_r['popt'][0]*normal(u_r['x'], *u_r['popt'][1:3]),
             label='Nube azul', c='b', ls='--')
    ax2.plot(u_r['x'], u_r['popt'][3]*normal(u_r['x'], *u_r['popt'][4:]),
             label='Sec. roja', c='r', ls='--')
    ax2.set_xlabel('$u-r$')
    ax2.legend(ncols=2)
    fig.savefig(_folder+'colors_fit.png')

    if save:
        np.savetxt(
            _folder+'colors_fit.dat', 
            np.vstack([
                u_r['popt'], 
                np.sqrt(np.diag(u_r['cov'])),
                g_r['popt'], 
                np.sqrt(np.diag(g_r['cov'])),
            ]).T,
            fmt='%8.6G', 
            header='u-r_popt u-r_perr g-r_popt g-r_perr',
            comments='#w1 mu1 sigma1 w2 mu2 sigma2 \n')

if __name__=='__main__':
    pass
