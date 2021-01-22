import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from astropy.stats import sigma_clip
import numpy.ma as ma

from astropy.io import ascii

fin = 'j0600_gVcorr.ecsv'
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

align_band = "R"


## CAUTION: resetting POBS deltamag so that fitting works. we'll fix it later.
mmm = (t['Observer Code']=="POBS")*(t['Band']==align_band)

t['deltamag'][mmm] = 0.0


# select all photometry in this band only
t_band = t[t['Band']==align_band]

print("Looking at band {}".format(align_band))

# get all unique names for observers in that waveband
t_band_by_obs = t_band.group_by('Observer Code').groups.keys

fig, (ax, ax2,ax3,ax4) = plt.subplots(4,1,figsize=(10,11))

print("All unique observers from band {}".format(align_band))
print(t_band_by_obs)

for currobs in t_band_by_obs:
    mask = (t['Band']==align_band) & (t['Observer Code']==currobs[0])
    tm = t[mask]
    ax.errorbar(tm['MJD'], tm['Magcorr'], yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=0.5, mec='none')

ax.legend()
ax.set_ylim(13.5,13.0)
ax.set_title('Initial uncorrected photometry from band {}'.format(align_band))

def find_offset(t, obs1, obs2, t_match, fiddle=0.0):
    # assume obs1 is the reference photometry data set, with more points
    # obs2 has only a few points
    # interpolate over the obs1 data set and move the obs2 data
    #### at the obs2 epochs, work out what the obs1 magnitude would be


    # MAKE SURE THAT obs1 has dm = 0!!!! otherwise 
    
    tin = t[(t['MJD']>t_match[0]) * (t['MJD']<t_match[1])* (t['Band']==align_band)]

    
    obs1_t = tin[(tin['Observer Code']==obs1)]

    print(obs1_t['deltamag'])

    # now, we might have points within the tin time range that are outside the range of values that obs1 has, which
    # means that the interpolation won't work outside these values

    obs1_t_max = np.max(obs1_t['MJD'])
    obs1_t_min = np.min(obs1_t['MJD'])

    # we clip out the obs2 observations that are within the obs1 observations so that the interpolation is guaranteed to be valid
    tin_obs2 = t[(t['MJD']>obs1_t_min) * (t['MJD']<obs1_t_max)* (t['Band']==align_band)]
    
    obs2_t = tin_obs2[(tin_obs2['Observer Code']==obs2)]

    print('number of points in {} lightcurve is {}'.format(obs1, len(obs1_t)))
    print('number of points in {} lightcurve is {}'.format(obs2, len(obs2_t)))

#    ax2.axhline(0.0, color='black', alpha=0.1)

    obs2_flux_at_obs1 = np.interp(obs2_t['MJD'], obs1_t['MJD'], obs1_t['Magnitude'], left=None, right=None, period=None)

    ax3.scatter(obs2_t['MJD'], obs2_flux_at_obs1)
    
    ax2.errorbar(obs2_flux_at_obs1,obs2_t['Magnitude']-obs2_flux_at_obs1, yerr=obs2_t['Uncertainty'], fmt='.')

    filtered_data = sigma_clip(obs2_t['Magnitude']-obs2_flux_at_obs1, sigma=2, maxiters=3)
    dm = filtered_data.mean()
    dm_err = filtered_data.std()

    str1 = 'dm = {:.3f} \pm {:.3f}$'.format(dm, dm_err)
    print(str1)
    
    ax2.axhline(dm, color='red')
    ax2.axhline(dm-dm_err, linestyle='dashed', color='red')
    ax2.axhline(dm+dm_err, linestyle='dashed', color='red')
    ax2.set_title('Measuring the magnitude offset between {} and {}'.format(obs1,obs2))

    # select all photometry in the band to be adjusted only
    mask = (t['Band']==align_band) * (t['Observer Code']==obs2)

    # MODIFY the table with the new deltamag and magcorr columns
    t['deltamag'][mask] = - dm + fiddle
    t['Magcorr'][mask]  = t['Magnitude'][mask] - dm + fiddle


    # plot the adjusted fits to check that they align
    tm = t[(t['Observer Code']==obs1)*(t['Band']==align_band)]
    ax3.errorbar(tm['MJD'], tm['Magcorr'], yerr=tm['Uncertainty'], fmt='.', label=obs1, alpha=0.5, mec='none')

    tm = t[(t['Observer Code']==obs2)*(t['Band']==align_band)]    
    ax3.errorbar(tm['MJD'], tm['Magcorr'], yerr=tm['Uncertainty'], fmt='.', label=obs2, alpha=0.5, mec='none')

    (ymi,yma) = ax3.get_ylim()
    sampl = Rectangle((t_match[0], ymi), t_match[1]-t_match[0], (yma-ymi),facecolor="black", alpha=0.1, zorder=-10)
    ax3.add_patch(sampl)
    ax3.legend()
    
    return (dm, dm_err) 

# no      ASASSN
# y         DFS
# y          HMB
# n          MAK
# n         NLX
# y        POBS
# n         SDM
# ref         TTG

(dm, dm_err) = find_offset(t, "POBS", "DFS", (58880,58920))
#(dm, dm_err) = find_offset(t, obs1,"DFS" , (58860, 58920))
#(dm, dm_err) = find_offset(t, obs1, "POBS", (58862,58884))
ax.get_shared_x_axes().join(ax, ax3)

def remove_obs_band(t,obs,band):
    mask = (t['Observer Code']==obs)*(t['Band']==band)
    tnew = t[~mask]
    return tnew

# remove observers where there are too few/incorrect points
t = remove_obs_band(t,"MLF","R")
t = remove_obs_band(t,"HMB","R")
t = remove_obs_band(t,"ASTEP","R")
t = remove_obs_band(t,"GSA","R")

# select all photometry in this band only
t_band = t[t['Band']==align_band]
# get all unique names for observers in that waveband
t_band_by_obs = t_band.group_by('Observer Code').groups.keys

# check that the Magcorr has the correct values in it!
for currobs in t_band_by_obs:
    mask = (t['Band']==align_band) & (t['Observer Code']==currobs[0])
    tm = t[mask]
    co = currobs[0]
    ax4.errorbar(tm['MJD'], tm['Magcorr'], yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=0.5, mec='none')

ax.get_shared_x_axes().join(ax, ax4)
    
ax4.legend()
ax4.set_ylim(15.2,14.0)
ax4.set_title('Corrected photometry from band {}'.format(align_band))

plt.draw()
plt.show()

#outf='j0600_gVRcorr.ecsv'
#print('Writing out {}.'.format(outf))
#t.write(outf, format='ascii.ecsv', overwrite=True)

