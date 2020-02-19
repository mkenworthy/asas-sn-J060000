import numpy as np
from pathlib import Path
from astropy.table import Table, vstack, unique
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.time import Time
from datetime import datetime
import j0600


fig, ax = plt.subplots()

(t, dm) = j0600.j0600()

# reject ones with large errors
t = t[t['Uncertainty']<0.06]


# read in the proposed straight line fit numbers

# run the emcee to fit the data

# use the data from the gp fit as starting positions for the other wavelengths




# pick out the gp data

tg = t[np.where(t['Band']=='gp')]
tg = tg[np.where(tg['Magnitude']>13.5)]


ax.errorbar(tg['MJD'],tg['Magcorr'], yerr=tg['Uncertainty'],fmt='b.')

#ax.vlines(xlines,0,20)
#ax.plot(xlines,ylines,'r')

plt.draw()
plt.show()

