import pandas as pd
import numpy as np
from bimodal import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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

#SDSS = pd.read_csv('SDSS_guia2_fmcaporaso.csv')

G = pd.read_fwf('gals.dat')
G.query('u_r < 4.0 and u_r > 0 and g_r < 1.3 and c9050 > 1.5 and c9050 < 4', inplace=True)

def plot_SDSSraw():

    # sky footprint
    fig, ax = plt.subplots(figsize=(6.47,4),subplot_kw={'projection': 'mollweide'})
    ax.grid(True, alpha=0.5)
    ax.scatter(np.deg2rad(SDSS.ra)-np.pi, np.deg2rad(SDSS.dec), s=1, alpha=1, c='C1')
    ax.set_xlabel('Ascención recta [deg]')
    ax.set_ylabel('Declinación [deg]')
    ax.set_xticklabels(np.arange(30,360,30))
    #fig.savefig('sky_footprint.pdf')
    #plt.show()

    # redshift distribution
    fig, ax = plt.subplots()
    ax.hist(SDSS.redshift, bins=30, histtype='step', color='dimgray', hatch='///')
    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel('Número de galaxias')
    #fig.savefig('z_dist_25bins.pdf')
    #plt.show()

    # mag vs z
    fig, ax = plt.subplots()
    ax.scatter(SDSS.redshift, SDSS.petroMag_r, s=1, alpha=0.5, c='C1')
    #ax.scatter([],[],c='C0',s=8,label='SDSS')
    ax.axhline(14.5, ls='--', c='k')
    ax.text(0.051,18.4,'$\\boldsymbol{r = 17.77}$',ha='right')
    ax.text(0.051,14.4,'$\\boldsymbol{r = 14.5}$',ha='right')
    ax.axhline(17.77, ls='--', c='k')
    ax.invert_yaxis()
    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel('Magnitud aparente $r$')
    #ax.legend(loc='upper right', frameon=True)
    #fig.savefig('r_vs_redshift.pdf')
    #plt.show()

def plot_volumecomplete():
    fig, ax = plt.subplots()
    ax.scatter(G.z, G.M_pet_r, s=1, alpha=0.9, c='C1')
    ax.invert_yaxis()
    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel('Magnitud absoluta $M_r$')

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
    #fig.savefig(_folder+'colors_fit.png')

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
    ax[2].scatter(G['c9050'], G['fracDeV_r'], s=1, alpha=0.4, c='C1')

    ax[0].set_xlabel('Índice de concentración $C$')
    ax[0].set_ylabel('Densidad de galaxias')
    ax[1].set_xlabel('$\\texttt{fracDeV_r}$')

    ax[2].set_xlabel('$C$')
    ax[2].set_ylabel('$\\texttt{fracDeV_r}$')


def plot_conc_u_r():
    fig, ax = plt.subplots(figsize=(5,5))
    ax.scatter(G['c9050'], G['u_r'], s=1, alpha=0.2, c='C1')
    ax.axvline(2.5, ls='--', lw=1., c='k') #
    ax.axhline(2.0, ls='--', lw=1., c='k') # color

    ax.set_xlabel('Índice de concentración $C$')
    ax.set_ylabel('$u-r$')

def plot_color_mag():
    blue = "#136cb4"
    red = "#e0282b"

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
    # cmap = ax.scatter(G['M_pet_r'], G['u_r'], s=5, c=G['c9050'], cmap=cmap, norm=colors.CenteredNorm(vcenter=2.5), alpha=1, facecolor=None)
    # fig.colorbar(cmap, label='$C$')

    ## == division by c=2.5
    ax.scatter(G['M_pet_r'][~mask], G['u_r'][~mask], s=1, marker='o', edgecolor=blue, facecolor='none', alpha=0.5)
    ax.scatter(G['M_pet_r'][mask], G['u_r'][mask], s=1, marker='o', edgecolor=red, facecolor='none', alpha=0.5)
    ### ax.scatter(G['M_pet_r'], G['u_r'], s=12, marker='o' if G['c9050']<2.5 else 's', edgecolor='C3' if G['c9050']<2.5 else 'C0', facecolor='none', alpha=0.5)

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
    mean_mag = [0.5*(G['M_pet_r'].min()+q[0]), 0.5*(q[0]+q[1]), 0.5*(q[1]+q[2]), 0.5*(q[2]+G['M_pet_r'].max())]

    fits = [fit_bimodal(G['u_r'][mask[i]], p0=[0.5, 1.5, 0.3, 0.5, 2.5, 0.2], nbins=np.linspace(0.65,3.25,50)) for i in range(4)]
    for i in range(4):
        ax.errorbar(mean_mag[i],fits[i]['popt'][1],fits[i]['popt'][2], fmt='o', markerfacecolor='b', elinewidth=1.5, ecolor='k', markeredgecolor='k')
        ax.plot(mean_mag[i],fits[i]['popt'][1],c='r')
        ax.errorbar(mean_mag[i],fits[i]['popt'][4],fits[i]['popt'][5], fmt='o', markerfacecolor='r', elinewidth=1.5, ecolor='k', markeredgecolor='k')
        ax.plot(mean_mag[i],fits[i]['popt'][4],c='b')

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


def plot_luminositysize():

    linear = lambda x, a,b: a*x+b

    s_early = '#4daf4a'
    s_late = '#984ea3'
    s_blue = "#136cb4"
    s_red = '#e41a1c'

    colors = [
    "#0b375c",
    "#690e10",
    "#5c1d66",
    "#1a6e17",]

    ### C grande -> luz concentrada en el centro -> early
    ### C chico -> luz desparramada -> late
    ### u-r grande -> mayor mag in r -> roja
    ### u-r chico -> mayor mag in u -> azul

    samples = [G[G['u_r']<2.0],
               G[G['u_r']>=2.0],
               G[G['c9050']<2.5],
               G[G['c9050']>=2.5]]

    fig, axes = plt.subplots(1,3, figsize=(12,4), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.1, hspace=0.2)

    axes[0].scatter(G['M_pet_r'], np.log10(G['r50']), s=1, c='dimgray', alpha=0.3)
    for i in range(4):
        (a,b), _ = curve_fit(linear, samples[i]['M_pet_r'], np.log10(samples[i]['r50']), p0=[1.0,-2.0])
        if i<2:
            axes[1].plot(samples[i]['M_pet_r'], samples[i]['M_pet_r']*a+b, c=colors[i], lw=2.0)
            axes[1].text(-22,0.95-0.08*i,f'${a=:2.2f}, {b=:2.2f}$',fontsize=10, c=colors[i])
        else:
            axes[2].plot(samples[i]['M_pet_r'], samples[i]['M_pet_r']*a+b, c=colors[i], lw=2.0)
            axes[2].text(-22,0.95-0.08*(i-2),f'${a=:2.2f}, {b=:2.2f}$',fontsize=10, c=colors[i])

    axes[1].scatter(samples[1]['M_pet_r'], np.log10(samples[1]['r50']), s=0.8, marker='o', c=s_red, alpha=0.3)
    axes[1].scatter(samples[0]['M_pet_r'], np.log10(samples[0]['r50']), s=0.8, marker='o', c=s_blue, alpha=0.3)
    axes[2].scatter(samples[3]['M_pet_r'], np.log10(samples[3]['r50']), s=0.8, marker='o', c=s_early, alpha=0.3)
    axes[2].scatter(samples[2]['M_pet_r'], np.log10(samples[2]['r50']), s=0.8, marker='o', c=s_late, alpha=0.3)

    for i in range(3):
        axes[i].set_xlabel('$M_r$ petrosiana')
    axes[0].set_ylabel(r'$\log_{10}(r_{50}/\mathrm{kpc})$')

    axes[1].plot([],[],'o', markersize=5, c=s_blue,label='$u-r<2.0$')
    axes[1].plot([],[],'o', markersize=5, c=s_red,label='$u-r\geq2.0$')

    axes[2].plot([],[],'o', markersize=5, c=s_late,label='$C<2.5$')
    axes[2].plot([],[],'o', markersize=5, c=s_early,label='$C\geq2.5$')

    axes[1].legend()
    axes[2].legend()

    fig2, ax = plt.subplots(figsize=(5,5))

    earlyred = samples[1][G['c9050']>=2.5]
    lateblue = samples[0][G['c9050']<2.5]

    ax.scatter(earlyred['M_pet_r'], np.log10(earlyred['r50']), s=0.8, marker='o', c='r', alpha=0.5)
    (a,b), _ = curve_fit(linear, earlyred['M_pet_r'], np.log10(earlyred['r50']), p0=[1.0,-2.0])
    ax.plot(earlyred['M_pet_r'], earlyred['M_pet_r']*a+b, c=colors[1], lw=2.0)
    ax.text(-22,0.95,f'${a=:2.2f}, {b=:2.2f}$',fontsize=12, c=colors[1])

    ax.scatter(lateblue['M_pet_r'], np.log10(lateblue['r50']), s=0.8, marker='o', c='b', alpha=0.5)
    (a,b), _ = curve_fit(linear, lateblue['M_pet_r'], np.log10(lateblue['r50']), p0=[1.0,-2.0])
    ax.plot(lateblue['M_pet_r'], lateblue['M_pet_r']*a+b, c=colors[0], lw=2.0)
    ax.text(-22,0.87,f'${a=:2.2f}, {b=:2.2f}$',fontsize=12, c=colors[0])

    ax.plot([],[],'o', markersize=5, c='r', label='Rojo + early')
    ax.plot([],[],'o', markersize=5, c='b', label='Azul + late')

    ax.set_ylabel(r'$\log\left({r_{50}/\mathrm{kpc}}\right)$')
    fig2.text(0.5, 0.01, '$M_r$ petrosiana', ha='center')

    ax.legend()


def plot_kormendy():
    # x = logr50, y = mu50

    linear = lambda x, a,b: a*x+b

    s_early = '#4daf4a'
    s_late = '#984ea3'
    s_blue = "#136cb4"
    s_red = '#e41a1c'

    colors = [
    "#0b375c",
    "#690e10",
    "#5c1d66",
    "#1a6e17",]

    ### C grande -> luz concentrada en el centro -> early
    ### C chico -> luz desparramada -> late
    ### u-r grande -> mayor mag in r -> roja
    ### u-r chico -> mayor mag in u -> azul

    samples = [G[(G['u_r']<2.0)&(G['mu50']<24.5)],
               G[(G['u_r']>=2.0)&(G['mu50']<24.5)],
               G[(G['c9050']<2.5)&(G['mu50']<24.5)],
               G[(G['c9050']>=2.5)&(G['mu50']<24.5)]]

    fig, axes = plt.subplots(1,3, figsize=(12,4), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.1, hspace=0.2)
    for i in range(3): axes[i].invert_yaxis()

    axes[0].scatter(np.log10(G['r50'][(G['mu50']<24.5)]), G['mu50'][(G['mu50']<24.5)], s=1, c='dimgray', alpha=0.3)
    # for i in range(4):
    #     (a,b), _ = curve_fit(linear, np.log10(samples[i]['r50']), samples[i]['mu50'], p0=[4.0,20.0])
    #     if i<2:
    #         axes[1].plot(np.log10(samples[i]['r50']), np.log10(samples[i]['r50'])*a+b, c=colors[i], lw=2.0)
    #         # axes[1].text(-22,0.95-0.08*i,f'${a=:2.2f}, {b=:2.2f}$',fontsize=10, c=colors[i])
    #     else:
    #         axes[2].plot(np.log10(samples[i]['mu50']), np.log10(samples[i]['mu50'])*a+b, c=colors[i], lw=2.0)
    #         # axes[2].text(-22,0.95-0.08*(i-2),f'${a=:2.2f}, {b=:2.2f}$',fontsize=10, c=colors[i])

    axes[1].scatter(np.log10(samples[1]['r50']), samples[1]['mu50'], s=0.8, marker='o', c=s_red, alpha=0.3)
    axes[1].scatter(np.log10(samples[0]['r50']), samples[0]['mu50'], s=0.8, marker='o', c=s_blue, alpha=0.3)
    axes[2].scatter(np.log10(samples[3]['r50']), samples[3]['mu50'], s=0.8, marker='o', c=s_early, alpha=0.3)
    axes[2].scatter(np.log10(samples[2]['r50']), samples[2]['mu50'], s=0.8, marker='o', c=s_late, alpha=0.3)

    for i in range(3):
        axes[i].set_xlabel(r'$\log_{10}(r_{50}/\mathrm{kpc})$')
    axes[0].set_ylabel(r'$\mu_{50} \, \mathrm{[mag/arcsec^2]}$')

    axes[1].plot([],[],'o', markersize=5, c=s_blue,label='$u-r<2.0$')
    axes[1].plot([],[],'o', markersize=5, c=s_red,label='$u-r\geq2.0$')

    axes[2].plot([],[],'o', markersize=5, c=s_late,label='$C<2.5$')
    axes[2].plot([],[],'o', markersize=5, c=s_early,label='$C\geq2.5$')

    axes[1].legend()
    axes[2].legend()

    fig2, ax = plt.subplots(figsize=(5,5))
    ax.invert_yaxis()

    earlyred = samples[1][G['c9050']>=2.5]
    lateblue = samples[0][G['c9050']<2.5]

    ax.scatter(np.log10(earlyred['r50']), earlyred['mu50'], s=0.8, marker='o', c='r', alpha=0.5)
    # (a,b), _ = curve_fit(linear, np.log10(earlyred['r50']), earlyred['mu50'], p0=[1.0,-2.0])
    # ax.plot(np.log10(earlyred['r50']), np.log10(earlyred['r50'])*a+b, c=colors[1], lw=2.0)
    # ax.text(-22,0.95,f'${a=:2.2f}, {b=:2.2f}$',fontsize=12, c=colors[1])

    ax.scatter(np.log10(lateblue['r50']), lateblue['mu50'], s=0.8, marker='o', c='b', alpha=0.5)
    # (a,b), _ = curve_fit(linear, np.log10(lateblue['r50']), lateblue['mu50'], p0=[4.0,20.0])
    # ax.plot(np.log10(lateblue['r50']), np.log10(lateblue['r50'])*a+b, c=colors[0], lw=2.0)
    # ax.text(-22,0.87,f'${a=:2.2f}, {b=:2.2f}$',fontsize=12, c=colors[0]k)

    ax.plot([],[],'o', markersize=5, c='r', label='Rojo + early')
    ax.plot([],[],'o', markersize=5, c='b', label='Azul + late')

    ax.set_ylabel(r'$\mu_{50} \, \mathrm{[mag/arcsec^2]}$')
    fig2.text(0.5, 0.01, r'$\log\left({r_{50}/\mathrm{kpc}}\right)$', ha='center')

    ax.legend()



if __name__=='__main__':

    #plot_SDSSraw()
    #plot_volumecomplete()
    plot_bimodalcolors('fig5_')
    #plot_conc_fracdev()
    #plot_conc_u_r()
    #plot_color_mag()
    #fit_colormag()
    #plot_luminositysize()
    #plot_kormendy()
    plt.show()

