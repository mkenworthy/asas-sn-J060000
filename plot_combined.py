import numpy as np
from astropy import constants as c
from astropy import units as u
import matplotlib.pyplot as plt
from datetime import datetime
from astropy.time import Time
from pathlib import Path
from astropy.table import Table, vstack, unique
from astropy.io import ascii
import matplotlib.pyplot as plt


params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)



dm = Table(
      dtype=('U8', 'f8', 'f8', 'f8', 'f8',  'f8', 'f8', 'f8', 'f8'),    meta={'name': 'delta magnitude table for Observers and J0600'},
       names=('Obs',  'B',  'V',  'R',  'I',  'SI', 'SR', 'gp', 'ip'))

#                      B     V     R     I     SI    SR    gp    ip
dm.add_row(['DFS',    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
dm.add_row(['GSA',    2.30, 0.00, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00])
dm.add_row(['HMB',    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
dm.add_row(['MAK',    0.00,-0.20, 0.00,-0.70, 0.00, 0.00,-0.04,-1.00])
dm.add_row(['MLF',    0.00, 0.00, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00])
dm.add_row(['NLX',    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
dm.add_row(['TTG',    0.00, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
dm.add_row(['VARA'  , 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,-0.40])
dm.add_row(['ASASSN', 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
dm.add_row(['MECK',   0.00, 0.00, 0.00, 0.00, 0.00, 0.00,-0.65,-0.42])

# how soon is now?
now = Time(datetime.utcnow(),format='datetime')

print('current MJD is {}'.format(now.mjd))

# aavso
from astropy.io import ascii
import matplotlib
aavso_file = 'data/aavso/aavsodata_5e3872681d138.txt'
t = ascii.read(aavso_file)

# read in JD into a Time object for easy conversion

t['MJD'] = t['JD'] - 2400000.5

times = Time(t['JD'], format='jd')

for f in Path('./').rglob('obs_*.ecsv'):
    print(f)

    tin = ascii.read(f)

    t = vstack([t,tin])

    
t.sort('MJD')
# g should be renamed to gp
t['Band'][np.where(t['Band']=='SI')] = 'ip'
t['Band'][np.where(t['Band']=='SR')] = 'rp'

# reject ones with large errors
t = t[t['Uncertainty']<0.05]

# SALT spectroscopc observations [HD to MJD]


# The observation spanned 27 Jan 22:54 - 23:26 UTC
#
salt_epoch = np.array([2458838.307512, 2458840.30739, 2458841.30588, 2458842.548183, 2458846.531863, 2458850.525301]) - 2400000.5

#plt.xlim(58830, now.mjd)
#plt.ylim(15.15,14.00)

#plt.vlines(salt_epoch, 15.15,14.00, linestyle='dotted', color='red')


# get a list of the unique bandpasses
t_by_band = t.group_by('Band')
print('all observed photometric bands:')
print(t_by_band.groups.keys)

# get list of unique observers from AAVSO
t_by_observer = t.group_by('Observer Code')

print('all observers:')
print(t_by_observer.groups.keys)

mybands = ('B','gp','V','R','I','ip') # ordered bands to plot

mag0 = {'U':14.0,'gp':14.2, 'B':14.88, 'V':13.615, 'R':12.887, 'I':12.25, 'ip':12.25}

band_wlen = {'U':400, 'gp':475, 'B':445, 'V':551, 'R':658, 'I':806, 'ip':769.8}
band_color = {'U':'blue', 'gp':'purple', 'B':'blue', 'V':'green', 'R':'red', 'I':'brown','ip':'brown'}

n_bands = len(mybands)

t_min = 58820
t_max = now.mjd

fig, axes = plt.subplots(n_bands, 1, figsize=(10, 16), sharex=True) # , sharey=True)
fig.subplots_adjust(hspace=0.05, wspace=0.05)



for (ax, band) in zip(axes, mybands): # loop through all the bands we want to plot
        print('current band is {}'.format(band))
        # photometric band label
        ax.text(0.1, 0.9, band[0], ha='center', va='center', fontsize=24, transform=ax.transAxes)
        
        ax.hlines(mag0[band], t_min, t_max, linestyle='dashed', color=band_color[band])

        # AAVSO table
        # we are plotting 'band' so let's get all the observers who have observed in that band
        t_currentband = t[t['Band']==band]
        #print('total number of measurements in {} band is {}'.format(band, len(ta_currentband)))
        
        if len(t_currentband)>0:
            # get unique names of all observers for that currentband
            t_currentband_by_obs = t_currentband.group_by('Observer Code').groups.keys
            
            for currobs in t_currentband_by_obs:
                mask = (t['Band']==band) & (t['Observer Code']==currobs[0])
                tm = t[mask]
                co = currobs[0]
                delta_mag = dm[dm['Obs']==co][band]

                ax.errorbar(tm['MJD'],tm['Magnitude']+delta_mag, yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=0.5, mec='none')

        ax.set_ylim(mag0[band]+0.9, mag0[band]-0.1)
        ax.legend(loc='lower left')

#        ax.plot(c['MJD'],c['g']-mag0['g']+mag0[band], color='grey')
#        ax.plot(c['MJD'],(c['g']-mag0['g'])*(band_wlen['g']/band_wlen[band])+mag0[band], color='red')


# set axis labels
axes[-1].set_xlabel('Time [MJD]')

ax.set_xlim(t_min, now.mjd)
  
fig.suptitle('J0600 Photometry', fontsize='x-large')

out = 'ASASSNV_J0600_MULTI_{:.1f}.pdf'.format(now.mjd)
plt.savefig(out,bbox_inches='tight')
