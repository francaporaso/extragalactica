import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots

plt.style.use('science')
_params = {
    'figure.figsize':(5,4),
    #'figure.figsize':(6.47,4), # golden ratio
    'figure.dpi':150,
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

def plot_kcorr(g):
    nrow, ncol = 5, 2
    band = ['u','g','r','i','z']
    fig, axes = plt.subplots(nrow, ncol, sharex=True, figsize=(7,5))
    for i in range(nrow):
       axes[i,0].plot(g['z_gal'], g[f'rk_p({i})'], ',', c='C0')
       axes[i,1].plot(g['z_gal'], g[f'rks_p({i})'], ',', c='C1')

    fig.supxlabel('Redshift $z$')
    # for i in range(2):
    #     axes[-1,i].set_xlabel('Redshift $z$')
    for j in range(5):
        axes[j,0].set_ylabel(f'$K_{{{band[j]}}}(z)$')
        axes[j,1].set_ylabel(f'$K_{{0.1{band[j]}}}(z)$')

    fig.align_ylabels(axes[:,0])
    fig.align_ylabels(axes[:,1])
    fig.subplots_adjust(wspace=0.3)

def plot_magz(g):
    fig, ax = plt.subplots()
    ax.plot(g['z_gal'], g['pet_r'], '.')

if __name__ == '__main__':
    columns = ['z_gal', 'pet_r', 'ext_r', 'r50']+[f'rk_p({k})' for k in range(5)]+[f'rks_p({k})' for k in range(5)]
    g = pd.read_fwf('mgs.dat', names=columns)
    g.query('z_gal<0.3', inplace=True)

    # plot_kcorr(g)
    # plt.savefig(_folder+'kcorr.png')

    plot_magz(g)
    plt.show()
