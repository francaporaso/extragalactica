import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

g = pd.read_fwf('distances.dat', names=['proj_den', 'ty'])

nbins = 10

fig, ax = plt.subplots()
ax.hist(np.log10(g.proj_den), bins=10)

E = g.query('ty == 1')
S0 = g.query('ty == 2')
S = g.query('ty == 3')

gx, b = np.histogram(np.log10(g['proj_den']), bins=nbins)

ell, _ = np.histogram(np.log10(E['proj_den']), bins=b)
lent, _ = np.histogram(np.log10(S0['proj_den']), bins=b)
spir, _ = np.histogram(np.log10(S['proj_den']), bins=b)

x = 0.5*(b[1:]+b[:-1])

print(ell/gx)

fig, ax = plt.subplots()
ax.plot(x, ell/gx, '.-', label='E')
ax.plot(x, lent/gx, '.-', label='S0')
ax.plot(x, spir/gx, '.-', label='S')

ax.legend()
plt.show()