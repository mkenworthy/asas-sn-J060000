import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from astropy.stats import sigma_clip
import numpy.ma as ma
from astropy.io import ascii


fin = 'J0600_all.ecsv'
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

mybands = ('B','gp','V','R','I','ip') # ordered bands to plot

mag0 = {'U':14.0,'gp':14.203, 'B':14.88, 'V':13.615, 'R':12.887, 'I':12.25, 'ip':12.25}

band_wlen = {'U':400, 'gp':475, 'B':445, 'V':551, 'R':658, 'I':806, 'ip':769.8}
band_color = {'U':'blue', 'gp':'purple', 'B':'blue', 'V':'green', 'R':'red', 'I':'brown','ip':'brown'}

n_bands = len(mybands)

align_band = "gp"

# select all photometry in this band only
t_band = t[t['Band']==align_band]

print("Looking at band {}".format(align_band))

# get all unique names for observers in that waveband
t_band_by_obs = t_band.group_by('Observer Code').groups.keys

fig, (ax, ax2,ax3) = plt.subplots(3,1,figsize=(10,11))

print("All unique observers from band {}".format(align_band))
print(t_band_by_obs)

for currobs in t_band_by_obs:
    mask = (t['Band']==align_band) & (t['Observer Code']==currobs[0])
    tm = t[mask]
    co = currobs[0]

    ax.errorbar(tm['MJD'], tm['Magnitude'], yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=0.5, mec='none')
#    ax.errorbar(tm['MJD'], tm['Magcorr'], yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=0.5, mec='none')

ax.legend()
ax.set_ylim(15.2,14.0)
ax.set_title('Initial uncorrected photometry from band {}'.format(align_band))
# identify region to do overlap

t_match = (58840, 58880)
t_match = (58840, 59200)

t_region = t[(t['MJD']>t_match[0]) * (t['MJD']<t_match[1])* (t['Band']==align_band)]
print(t_region)
#ax.errorbar(t_region['MJD'], t_region['Magnitude'], yerr=t_region['Uncertainty'], fmt='.', alpha=0.5, mec='none')

# draw a grey rectangle to identify overlap

obs1 = "ASASSN"
obs2 = "MAK"

obs1_t = t_region[(t_region['Observer Code']==obs1)]
obs2_t = t_region[(t_region['Observer Code']==obs2)]

(ymi,yma) = ax.get_ylim()


sampl = Rectangle((t_match[0], ymi), t_match[1]-t_match[0], (yma-ymi),facecolor="black", alpha=0.1, zorder=-10)

# Add collection to axes
ax.add_patch(sampl)

# count which observer has the more data points - they will have an interpolator to the observer with fewer points

print('number of points in {} lightcurve is {}'.format(obs1, len(obs1_t)))

print('number of points in {} lightcurve is {}'.format(obs2, len(obs2_t)))

ax2.axhline(0.0, color='black', alpha=0.1)

#obs1_flux_at_obs2 = np.interp(obs1_t['MJD'], obs2_t['MJD'], obs2_t['Magcorr'], left=None, right=None, period=None)
#ax2.errorbar(obs1_flux_at_obs2,obs1_t['Magcorr']-obs1_flux_at_obs2,yerr=obs1_t['Uncertainty'],fmt='.')
#filtered_data = sigma_clip(obs1_t['Magcorr']-obs1_flux_at_obs2, sigma=3, maxiters=3)

obs1_flux_at_obs2 = np.interp(obs1_t['MJD'], obs2_t['MJD'], obs2_t['Magnitude'], left=None, right=None, period=None)
ax2.errorbar(obs1_flux_at_obs2,obs1_t['Magnitude']-obs1_flux_at_obs2,yerr=obs1_t['Uncertainty'],fmt='.')

filtered_data = sigma_clip(obs1_t['Magnitude']-obs1_flux_at_obs2, sigma=3, maxiters=3)
mag0_asas_gp = filtered_data.mean()
mag0_asas_gp_err = filtered_data.std()

str1 = 'ASAS $g\' = {:.3f} \pm {:.3f}$'.format(mag0_asas_gp, mag0_asas_gp_err)
print(str1)
    
ax2.axhline(mag0_asas_gp, color='red')
ax2.axhline(mag0_asas_gp-mag0_asas_gp_err, linestyle='dashed', color='red')
ax2.axhline(mag0_asas_gp+mag0_asas_gp_err, linestyle='dashed', color='red')
ax2.set_title('Measuring the magnitude offset between {} and {}'.format(obs1,obs2))


# select the obs2 data in that filter, put the value into the deltamag column, then put the corrected value into the Magcorr column

# select all photometry in the band to be adjusted only
mask = (t['Band']==align_band) * (t['Observer Code']==obs2)

t['deltamag'][mask] = mag0_asas_gp

t['Magcorr'][mask] = t['Magnitude'][mask] + mag0_asas_gp

#outf='tmp1.ecsv'
#print('Writing out {}.'.format(outf))
#t.write(outf, format='ascii.ecsv', overwrite=True)

# check that the Magcorr has the correct values in it!
for currobs in t_band_by_obs:
    mask = (t['Band']==align_band) & (t['Observer Code']==currobs[0])
    tm = t[mask]
    co = currobs[0]
    ax3.errorbar(tm['MJD'], tm['Magcorr'], yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=0.5, mec='none')

ax.get_shared_x_axes().join(ax, ax3)
    
ax3.legend()
ax3.set_ylim(15.2,14.0)
ax3.set_title('Corrected photometry from band {}'.format(align_band))

plt.draw()
plt.show()

outf='j0600_gcorr.ecsv'
print('Writing out {}.'.format(outf))
t.write(outf, format='ascii.ecsv', overwrite=True)

