import pandas as pd
import numpy as np
from bimodal import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scienceplots

plt.style.use('science')
_params = {
    'figure.figsize':(5,3),
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
}
plt.rcParams.update(_params)

_folder = 'figures/'

G = pd.read_fwf('gals.dat')
G.query('u_r < 4.0 and u_r > 0 and g_r < 1.3 and c9050 > 1.5 and c9050 < 4', inplace=True)

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

def plot_bimodalcolors(savename=None):

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

    if savename is not None:
        np.savetxt(
            _folder+savename+'colors_fit.dat',
            np.vstack([u_r['popt'], np.sqrt(np.diag(u_r['cov'])),
                       g_r['popt'], np.sqrt(np.diag(g_r['cov']))]).T,
            header='u-r_popt u-r_perr g-r_popt g-r_perr',
            comments='#w1 mu1 sigma1 w2 mu2 sigma2 \n',
            fmt='%8.6G', 
        )

def plot_conc_fracdev():
    fig, ax = plt.subplots(1,3, figsize=(12,3.8))
    ax[0].hist(G['c9050'], bins=50, density=True, histtype='step', hatch='///', color='dimgray')
    ax[1].hist(G['fracDeV_r'], bins=50, density=True, histtype='step', hatch='///', color='dimgray')
    ax[2].scatter(G['c9050'], G['fracDeV_r'], s=1, alpha=0.4)

    ax[0].set_xlabel('Índice de concentración $C$')
    ax[0].set_ylabel('Densidad de galaxias')
    ax[1].set_xlabel('$\\texttt{fracDeV_r}$')

    ax[2].set_xlabel('$C$')
    ax[2].set_ylabel('$\\texttt{fracDeV_r}$')


def plot_conc_u_r():
    fig, ax = plt.subplots(figsize=(5,5))
    ax.scatter(G['c9050'], G['u_r'], s=1, alpha=0.2)
    ax.axvline(2.5, ls=':', lw=1.5, c='dimgray') # 
    ax.axhline(2.0, ls=':', lw=1.5, c='dimgray') # color
    
    ax.set_xlabel('Índice de concentración $C$')
    ax.set_ylabel('$u-r$')
    
def plot_color_mag():
    color_min    = "b"
    color_center = "lightgray"
    color_max    = "r"
    cmap = colors.LinearSegmentedColormap.from_list(
        "cmap_name",
        [color_min, color_center, color_max]
    )

    mask = G['c9050'] > 2.5
    fig, ax = plt.subplots(figsize=(5,5))
    
    ## == color bar of concentration
    cmap = ax.scatter(G['M_pet_r'], G['u_r'], s=5, c=G['c9050'], cmap=cmap, norm=colors.CenteredNorm(vcenter=2.5), alpha=1, facecolor=None)
    fig.colorbar(cmap, label='$C$')

    ## == division by c=2.5
    # ax.scatter(G['M_pet_r'][~mask], G['u_r'][~mask], s=12, marker='^', edgecolor='C0', facecolor='none', alpha=0.5)
    # ax.scatter(G['M_pet_r'][mask], G['u_r'][mask], s=12, marker='o', edgecolor='C3', facecolor='none', alpha=0.5)
    #ax.scatter(G['M_pet_r'], G['u_r'], s=12, marker='o' if G['c9050']<2.5 else 's', edgecolor='C3' if G['c9050']<2.5 else 'C0', facecolor='none', alpha=0.5)

    ## == hexbin
    #ax.hexbin(G['M_pet_r'], G['u_r'], gridsize=50, bins='log', cmap='binary')

    ax.axhline(2.0, ls='--', c='k', alpha=0.8)

    #ax.invert_xaxis()
    ax.set_xlabel('$M_r$ petrosiana')
    ax.set_ylabel('$u-r$')
    return fig, ax

def fit_colormag():
    q = np.quantile(G['M_pet_r'], [0.25,0.50,0.75])
    print(q)
    fig, ax = plot_color_mag()
    
    for i in range(3):
        ax.axvline(q[i], ls='--', c='k', alpha=0.8)

    mask = [
        G['M_pet_r']<q[0], 
        (G['M_pet_r']>q[0])&(G['M_pet_r']<q[1]), 
        (G['M_pet_r']>q[1])&(G['M_pet_r']<q[2]), 
        G['M_pet_r']>q[2]
    ]

    quantiles = [f'$M_r<{q[0]:2.2f}$',f'${q[0]:2.2f}<M_r<{q[1]:2.2f}$',f'${q[1]:2.2f}<M_r<{q[2]:2.2f}$',f'$M_r>{q[2]:2.2f}$']

    fits = [fit_bimodal(G['u_r'][mask[i]], p0=[0.5, 1.5, 0.3, 0.5, 2.5, 0.2], nbins=np.linspace(0.65,3.25,50)) for i in range(4)]

    fig2, axes = plt.subplots(2,2, sharex=True, sharey=True, figsize=(5,5))
    fig2.subplots_adjust(wspace=0.1, hspace=0.2)
    axesflat = axes.flatten()
    for i in range(4):
        axesflat[i].axvline(2.0,ls='--',c='k',alpha=0.8)
                
        axesflat[i].stairs(edges=fits[i]['xedges'], values=fits[i]['y'],
                hatch='///', color='dimgray')
        axesflat[i].plot(fits[i]['x'], fits[i]['yfit'],
                c='k', lw=1.8, alpha=0.8)
        axesflat[i].plot(fits[i]['x'], fits[i]['popt'][0]*normal(fits[i]['x'], *fits[i]['popt'][1:3]),
                c='b', ls='--')
        axesflat[i].plot(fits[i]['x'], fits[i]['popt'][3]*normal(fits[i]['x'], *fits[i]['popt'][4:]),
                c='r', ls='--')
        axesflat[i].set_title(quantiles[i], fontsize=9)

    axesflat[1].plot([],[], c='r', label='Sec. roja')
    axesflat[1].plot([],[], c='b', label='Nube azul')
    axesflat[1].plot([],[], c='k', label='Ajuste bimodal')
    #axesflat[1].stairs([],[], c='dimgray', hatch='///', label='SDSS')
    axesflat[1].legend()

    fig2.text(0.5, 0.04, '$u-r$', ha='center')
    fig2.text(0.04, 0.5, 'Densidad de galaxias', va='center', rotation='vertical')


def plot_sizemag():

    s_early='#81c97f'
    s_late='#b980c2'
    s_blue='#6da5d3'
    s_red='#ed5e60'
    
    early = '#4daf4a'
    late = '#984ea3'
    blue = "#136cb4"
    red = '#e41a1c'

    fig, axes = plt.subplots(1,3, figsize=(12,4))

    axes[0].scatter(G['M_pet_r'], np.log10(G['r50']), s=10, c='C0', alpha=0.3)
    axes[1].scatter(G['M_pet_r'][G['u_r']<2.0], np.log10(G['r50'][G['u_r']<2.0]), s=10, marker='^', c=s_blue, alpha=0.5)
    axes[1].scatter(G['M_pet_r'][G['u_r']>2.0], np.log10(G['r50'][G['u_r']>2.0]), s=10, marker='s', c=s_red, alpha=0.3)
    axes[2].scatter(G['M_pet_r'][G['c9050']<2.5], np.log10(G['r50'][G['c9050']<2.5]), s=10, marker='^', c=s_late, alpha=0.5)
    axes[2].scatter(G['M_pet_r'][G['c9050']>2.5], np.log10(G['r50'][G['c9050']>2.5]), s=10, marker='s', c=s_early, alpha=0.3)

    axes[1].plot(x:=np.arange(G['M_pet_r'].min(),G['M_pet_r'].max(),1), -0.01*x, c=blue)
    axes[1].plot(x:=np.arange(G['M_pet_r'].min(),G['M_pet_r'].max(),1), -0.01*x+0.3, c=red)

    axes[2].plot(x:=np.arange(G['M_pet_r'].min(),G['M_pet_r'].max(),1), -0.01*x, c=late)
    axes[2].plot(x:=np.arange(G['M_pet_r'].min(),G['M_pet_r'].max(),1), -0.01*x+0.3, c=early)

if __name__=='__main__':

    #plot_color_mag()
    #fit_colormag()
    plot_sizemag()
    plt.show()
    print('no está listo...')