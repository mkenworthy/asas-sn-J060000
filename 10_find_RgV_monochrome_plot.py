import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from astropy.stats import sigma_clip
import numpy.ma as ma

from astropy.io import ascii

fin = 'j0600_gVRcorr.ecsv'
print('reading in {}'.format(fin))
t = ascii.read(fin)

# reject ones with large errors
t = t[t['Uncertainty']<0.05]

# get a list of the unique bandpasses
t_by_band = t.group_by('Band')
print('all observed photometric bands:')
print(t_by_band.groups.keys)

# get list of unique observers from AAVSO
t_by_observer = t.group_by('Observer Code')

print('all observers:')
print(t_by_observer.groups.keys)

# this band is the one we scale all the other bands to
align_band = "gp"

# select all photometry in this band only
t_band = t[t['Band']==align_band]

print("Looking at band {}".format(align_band))

# get all unique names for observers in that waveband
t_band_by_obs = t_band.group_by('Observer Code').groups.keys

fig, (ax,ax2,ax3,ax4) = plt.subplots(4,1,figsize=(20,11), sharex=True)

print("All unique observers from band {}".format(align_band))
print(t_band_by_obs)

for currobs in t_band_by_obs:
    mask = (t['Band']==align_band) & (t['Observer Code']==currobs[0])
    tm = t[mask]
    ax.errorbar(tm['MJD'], tm['Magcorr'], yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=0.5, mec='none')

ax.legend()
ax.set_ylim(15.2,14.0)
ax.set_xlim(58500., 59300.)
ax.set_title('Photometry from band {}'.format(align_band))


xband = "R"

# select all photometry in this band only
t_xband = t[t['Band']==xband]

print("Looking at band {}".format(xband))

# get all unique names for observers in that waveband
t_xband_by_obs = t_xband.group_by('Observer Code').groups.keys

for currobs in t_xband_by_obs:
    mask = (t['Band']==xband) & (t['Observer Code']==currobs[0])
    tm = t[mask]
    ax2.errorbar(tm['MJD'], tm['Magcorr'], yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=0.5, mec='none')

ax2.legend()
ax2.set_ylim(13.5,12.8)
#ax2.set_xlim(58500., 59300.)
ax2.set_title('Photometry from band {}'.format(xband))

def mag_to_tx(m,m0,scale=1.0):
    expo = ((m0-m)*scale)/2.5
    x = np.power(10.,expo)
    return x

# first select all align_band
t_band = t[t['Band']==align_band]
tx_align_band = mag_to_tx(t_band['Magcorr'], 14.2)
tx_xband = mag_to_tx(t_xband['Magcorr'], 12.85,1.6)

ax3.errorbar(t_band['MJD'], tx_align_band, yerr=t_band['Uncertainty'],
             fmt='.',
             color='blue',
             label='gp normalised',
             alpha=0.5,
             mec='none')

ax3.errorbar(t_xband['MJD'], tx_xband, yerr=t_xband['Uncertainty'],
             fmt='.',
             color='brown',
             label='R band normalised',
             alpha=0.5,
             mec='none')

ax3.set_ylim(0.0, 1.1)
#ax3.set_xlim(58500., 59300.)
ax3.set_title('Photometry scaled to unity transmission')
ax3.legend()

ax4.errorbar(t_band['MJD'], tx_align_band, yerr=t_band['Uncertainty'],
             fmt='.',
             color='black',
             alpha=0.5,
             mec='none')

ax4.errorbar(t_xband['MJD'], tx_xband, yerr=t_xband['Uncertainty'],
             fmt='.',
             color='black',
             alpha=0.5,
             mec='none')

ax4.set_ylim(0.0, 1.1)
#ax4.set_xlim(58500., 59300.)
ax4.set_title('Scaled photometry all combined')


import sys, os, datetime

pathname = os.path.dirname(sys.argv[0])

tagline = os.path.abspath(pathname) + '/' + \
        sys.argv[0] + ' on ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

ax.text(0.99, 0.01, tagline, ha='right', va='bottom', \
        transform=ax.transAxes, color='black', fontsize=10)
plt.savefig('plot_10_find_monochrome_1.pdf')

fig2, (ax5) = plt.subplots(1,1,figsize=(20,6))



ax5.errorbar(t_band['MJD'], tx_align_band, yerr=t_band['Uncertainty'],
             fmt='.',
             color='black',
             alpha=0.3,
             mec='none')

ax5.errorbar(t_xband['MJD'], tx_xband, yerr=t_xband['Uncertainty'],
             fmt='.',
             color='black',
             alpha=0.3,
             mec='none')

ax5.set_ylim(0.2, 1.1)
ax5.set_xlim(58700., 59300.)
ax5.set_xlabel('Epoch [MJD]')
ax5.set_title('J0600 scaled and normalised photometry')


ax5.text(0.99, 0.01, tagline, ha='right', va='bottom', \
        transform=ax5.transAxes, color='black', fontsize=10)

plt.savefig('plot_10_find_monochrome_2.pdf')


plt.show()
