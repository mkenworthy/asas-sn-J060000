import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.io import ascii

# TODO set the

# UC4 296-009008
# (90.140913, -30.99747)
# ASASSN ID AP 18326421 V=12.79

def find_min(xobj,yobj, xlist, ylist):
    'find the index of the star closest to position xobj, yobj'
    dr = (xlist-xobj)**2 + (ylist-yobj)**2
    # find min dr 
    ind = np.argmin(dr)
    print(np.sqrt(dr[ind]))
    return ind

from astropy.table import Table
from pathlib import Path
from astropy.io import ascii

# read in aavso table as a template
d = ascii.read('00demo_aavso.txt')

outf = 'obs_LCOGT_MAK.ecsv'

# get a list of all the AAVSO files in the directory

my_file = Path(outf)

if my_file.is_file():
    # file exists
    print('{} exists! Reading in...'.format(my_file))
    t = ascii.read(outf)
else: # if it does not exist, make it for the first time

    print('{} does not exist. Processing all the files.'.format(my_file))
    t = ascii.read('00demo_aavso.txt')
    t['lco_filename']=['verylongfilenameusedtomakestaaasduffgoon']
    t['MJD'] = 12345.00
    t['LAT'] = 0.00
    t['LON'] = 0.00
    t.replace_column('Band',['3lf'])
    t.replace_column('Star Name',['ASASSN-V J060000.76-310027.83'])
    t.replace_column('Comp Star 1',['UC4 296-009008'])
    t.remove_column('Charts')
    t.remove_row(0)
#    print(t.columns)


    # process all the files in the list and add them to the table

for imname in Path('./').rglob('*e91.fits.fz'):

    # is imname already in the table under the lco_filename
    if (t['lco_filename'].size > 0):
        if (np.any(t['lco_filename'] == imname.name)):
            print('{} already processed. skipping...'.format(imname))
            continue

    
    hdu1 = fits.open(imname)
    with fits.open(imname) as hdul:
        print('processing: ',imname)

        imh = hdu1['SCI'].header
        
        epoch = imh['MJD-OBS']
        print('{:.5f} MJD'.format(epoch))
        print(imname)
        t.add_row()
        t[-1]['Band']         = imh['FILTER']
        t[-1]['MJD']          = imh['MJD-OBS']
        t[-1]['lco_filename'] = imname
        t[-1]['Observer Code'] = 'MAK'
        t[-1]['Airmass']      = 0.0
        t[-1]['Observer Affiliation'] = imh['SITEID']
        t[-1]['LAT'] = imh['LATITUDE']
        t[-1]['LON'] = imh['LONGITUD']
        t[-1]['Airmass'] = imh['AIRMASS']
        t[-1]['Star Name'] = 'ASASSN-V J060000.76-310027.83'
        t[-1]['Comp Star 1'] = 'UC4 296-009008'

        phot = hdu1['CAT'].data
        # calculate the flux of the star relative to the STD star
        
        # Parse the WCS keywords in the primary HDU
        w = wcs.WCS(hdu1['SCI'].header)

        # Print out the "name" of the WCS, as defined in the FITS header
        #print(w.wcs.name)
               # Print out all of the settings that were parsed from the header
        #w.wcs.print_contents()

        # star coords in pix
        star_crd = w.wcs_world2pix([[90.0033, -31.0076]], 0) # star
        # RA
        #90.0033 (06:00:00.792)
        #Dec
        #-31.0076 (-31:00:27.36)
        #Epoch
        #2000

        # reference star UC4 296-009008
        refe_crd = w.wcs_world2pix([[90.1409130, -30.9974700]], 0) # 296-009008
        #print(refe_crd)
        star_arg = find_min(star_crd[0][0], star_crd[0][1], phot['x'], phot['y'])
        #print(star_arg)
        refe_arg = find_min(refe_crd[0][0], refe_crd[0][1], phot['x'], phot['y'])
        #print(refe_arg)

        targflux = phot['flux'][star_arg]
        targfluxerr = phot['fluxerr'][star_arg]
        refflux = phot['flux'][refe_arg]
        reffluxerr = phot['fluxerr'][refe_arg]

        print('star flux is {:.1f} +- {:.1f}'.format(targflux, targfluxerr))
        print('refe flux is {:.1f} +- {:.1f}'.format(refflux, reffluxerr))

        dm = -2.5 * np.log10(targflux/refflux)
        dmerr = np.sqrt( (targfluxerr/targflux)**2 + (reffluxerr/refflux)**2) * dm
        print()
        print('delta mag is {:.3f}+-{:.3f} mag'.format(dm, dmerr))

        mag = dm + 13.054 # 13.054 is g' magnitude eyeballed from ASAS-SN

        t[-1]['Magnitude'] = mag
        t[-1]['Uncertainty'] = dmerr

t.sort('MJD')


# write out the table
t.write(outf, format='ascii.ecsv', overwrite=True)
print(t)


