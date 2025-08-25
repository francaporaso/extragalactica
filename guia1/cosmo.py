import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

def cosmological_distances():
    cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3)

    z = np.arange(0.0,6.0,0.001)
    df = pd.DataFrame({
        'z':z,
        'chi':cosmo.comoving_distance(z),
        'd_L':cosmo.luminosity_distance(z),
        'd_A':cosmo.angular_diameter_distance(z)
    })

    df.to_csv('cosmo_astropy.dat', index=False, sep=',')

    return

if __name__ == '__main__':
    cosmological_distances()
