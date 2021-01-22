import numpy as np
from j0600 import *
import matplotlib.pyplot as plt
from kenworthy import tag_plot

(t, dm) = j0600()

# reject ones with large errors
t = t[t['Uncertainty']<0.05]

align_band = "gp"

# select all photometry in this band only

t_band = t[t['Band']==align_band]

print("Looking at band {}".format(align_band))

# select only the ASASSN photometry

t_asas_gp = t_band[t_band['Observer Code']=='ASASSN']

fig, (ax) = plt.subplots(1,1,figsize=(10,6))

ax.errorbar(t_asas_gp['MJD'], t_asas_gp['Magcorr'], yerr=t_asas_gp['Uncertainty'], fmt='.', label="ASAS-SN $g\'$", alpha=0.5, mec='none')
ax.set_ylim(15.25,14.00)
ax.legend()

# identify region to do statistics on
t_match = (58400, 58600)
t_match = (58000, 58750)

# plot that on the graph

(ymi,yma) = ax.get_ylim()
from matplotlib.patches import Rectangle

sampl = Rectangle((t_match[0], ymi), t_match[1]-t_match[0], (yma-ymi),facecolor="black", alpha=0.1)

# Add collection to axes
ax.add_patch(sampl)
ax.set_xlabel('Epoch [MJD]')
ax.set_ylabel('ASAS SN $g\'$')
t_region = t_asas_gp[(t_asas_gp['MJD']>t_match[0]) * (t_asas_gp['MJD']<t_match[1])]

from astropy.stats import sigma_clip
import numpy.ma as ma

filtered_data = sigma_clip(t_region['Magcorr'], sigma=2, maxiters=5)

mag0_asas_gp = filtered_data.mean()
mag0_asas_gp_err = filtered_data.std()

str1 = 'ASAS $g\' = {:.3f} \pm {:.3f}$'.format(mag0_asas_gp, mag0_asas_gp_err)
print(str1)

ax.axhline(mag0_asas_gp, color='red')
ax.axhline(mag0_asas_gp-mag0_asas_gp_err, linestyle='dashed', color='red')
ax.axhline(mag0_asas_gp+mag0_asas_gp_err, linestyle='dashed', color='red')

import sys, os, datetime

pathname = os.path.dirname(sys.argv[0])

tagline = os.path.abspath(pathname) + '/' + \
        sys.argv[0] + ' on ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

ax.text(0.99, 0.01, tagline, ha='right', va='bottom', \
        transform=ax.transAxes, color='black', fontsize=10)

ax.text(0.1, 0.5, str1, ha='left', va='center', \
        transform=ax.transAxes, color='black', fontsize=16)

plt.draw()
plt.savefig('find_asassn_gp_baseline.pdf')
plt.show()

quit()

