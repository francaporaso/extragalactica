from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots

parser = ArgumentParser()
parser.add_argument('--save', action='store_true')
args = parser.parse_args()

plt.style.use('science')
_params = {
    'figure.figsize':(5,5),
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

df = pd.read_fwf('dressler_data.dat', names=['log_den_low', 'log_den_high', 'log_den', 'tot', 'frac_ell', 'frac_s0', 'frac_spir'])
edges = np.concatenate([df['log_den_low'], [df['log_den_high'].to_numpy()[-1]]])

fig, ax = plt.subplots(2, 1, sharex=True, height_ratios=[0.4,1])

ax[0].stairs(df['tot'], edges, color='k')
ax[0].set_ylabel('Número total')
ax[0].semilogy()

ax[1].plot(df['log_den'], df['frac_ell'], '.-', c="#FF1F2A", label='E')
ax[1].plot(df['log_den'], df['frac_s0'], '.--', c="#00CD44", label='S0')
ax[1].plot(df['log_den'], df['frac_spir'], '.-.', c="#007ADE", label='S+Irr')
ax[1].set_xlabel(r'$\log \rho_{\mathrm{proj}}$')
ax[1].set_ylabel('Fracción de galaxias')
ax[1].legend()
fig.subplots_adjust(hspace=0.05)
fig.align_ylabels()

if args.save:
    fig.savefig(_folder+'dressler_plot.png')
else:
    plt.show()
