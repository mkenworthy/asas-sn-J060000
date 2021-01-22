import numpy as np
from j0600 import *
import matplotlib.pyplot as plt

(t, dm) = j0600()

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

# get the ASAS-SN gp band data - we're treating this as the baseline calibrated data
t_asas = t[(t['Band']=='gp') & (t['Observer Code']=='ASASSN')]
fig3, (ax3) = plt.subplots(1,1,figsize=(10,6))

ax3.errorbar(t_asas['MJD'], t_asas['Magcorr'], yerr=t_asas['Uncertainty'], fmt='.', label="ASAS $g\'$", alpha=0.5, mec='none')
ax3.legend()
ax3.set_ylim(15.2,14.0)

#plt.show()

align_band = "V"

# select all photometry in this band only
t_band = t[t['Band']==align_band]

print("Looking at band {}".format(align_band))

# get all unique names for observers in that waveband
t_band_by_obs = t_band.group_by('Observer Code').groups.keys

fig, (ax, ax2) = plt.subplots(2,1,figsize=(10,6))

print("All unique observers from band {}".format(align_band))

# offset plot for all the observers
for (n, currobs) in enumerate(t_band_by_obs):
    mask = (t['Band']==align_band) & (t['Observer Code']==currobs[0])
    tm = t[mask]
    co = currobs[0]

    ax.errorbar(tm['MJD'], tm['Magcorr']+0.5*n, yerr=tm['Uncertainty'], fmt='.', label=currobs[0], alpha=1.0, mec='none')

ax.set_title('all observed photometry in {} band offset by observer names'.format(align_band))
ax.legend()
ax.set_ylim(17.6,13.5)
ax.set_xlim(58300, 59300)
plt.show()
quit()

# choose the most accurate photometry and align everyone to that person
#
# TTG has the best photometry with the greatest coverage
#
align_obs = "TTG"


#
#
#
#
#
#
#
#
#

# identify region to do overlap

t_match = (58840, 58880)

t_region = t[(t['MJD']>t_match[0]) * (t['MJD']<t_match[1])* (t['Band']==align_band)]
print(t_region)
ax.errorbar(t_region['MJD'], t_region['Magcorr'], yerr=t_region['Uncertainty'], fmt='.', alpha=0.5, mec='none')

# draw a red rectangle to identify overlap

obs1 = "ASASSN"
obs2 = "MAK"

obs1_t = t_region[(t_region['Observer Code']==obs1)]
obs2_t = t_region[(t_region['Observer Code']==obs2)]

(ymi,yma) = ax.get_ylim()
from matplotlib.patches import Rectangle

    
from astropy.stats import sigma_clip
import numpy.ma as ma

sampl = Rectangle((t_match[0], ymi), t_match[1]-t_match[0], (yma-ymi),facecolor="black", alpha=0.1, zorder=-5)

# Add collection to axes
ax.add_patch(sampl)

# count which observer has the more data points - they will have an interpolator to the observer with fewer points

print(len(obs1_t))

print(len(obs2_t))

ax2.axhline(0.0, color='black', alpha=0.1)

if (len(obs2_t)>len(obs1_t)):
    obs1_flux_at_obs2 = np.interp(obs1_t['MJD'], obs2_t['MJD'], obs2_t['Magcorr'], left=None, right=None, period=None)
    ax2.errorbar(obs1_flux_at_obs2,obs1_t['Magcorr']-obs1_flux_at_obs2,yerr=obs1_t['Uncertainty'],fmt='.')

    filtered_data = sigma_clip(obs1_t['Magcorr']-obs1_flux_at_obs2, sigma=3, maxiters=3)

    mag0_asas_gp = filtered_data.mean()
    mag0_asas_gp_err = filtered_data.std()

    str1 = 'ASAS $g\' = {:.3f} \pm {:.3f}$'.format(mag0_asas_gp, mag0_asas_gp_err)
    print(str1)
    
    ax2.axhline(mag0_asas_gp, color='red')
    ax2.axhline(mag0_asas_gp-mag0_asas_gp_err, linestyle='dashed', color='red')
    ax2.axhline(mag0_asas_gp+mag0_asas_gp_err, linestyle='dashed', color='red')

plt.show()

