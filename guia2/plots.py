import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

params = {
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
plt.rcParams.update(params)

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

if __name__=='__main__':
    plot_problema2()
