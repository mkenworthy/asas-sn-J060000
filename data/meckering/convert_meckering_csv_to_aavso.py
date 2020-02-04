import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, vstack, unique
from pathlib import Path

# read in aavso table as a template
d = ascii.read('00demo_aavso.txt')

outf = 'obs_MECK.ecsv'

# get a list of all the AAVSO files in the directory

t = ascii.read('00demo_aavso.txt')
t['lco_filename']=['verylongfilenameusedtomakestaaasduffgoon']
t['MJD'] = 12345.00
t['LAT'] = 0.00
t['LON'] = 0.00
t['Airmass']= 0.000
t.replace_column('Band',['3lf'])
t.replace_column('Star Name',['ASASSN-V J060000.76-310027.83'])
t.replace_column('Comp Star 1',['UC4 296-009008'])
t.remove_column('Charts')
t.remove_row(0)

# process all the files in the list and add them to the table
# grep all of gary saccos text files
for f in Path('./').rglob('simplified.csv'):
    
    print(f.name)

    # convert input table to columns in template table
    tin = ascii.read(f, format='csv', data_start=2)
    tin.remove_rows(np.where(tin['Mean'].mask))
#    tin.remove_column('Limit')
    tin.rename_column('Filter','Band')
    tin.rename_column('Mean','Magnitude')
    tin['mag_err']=0.01
    
    tin['MJD'] = tin['JD (Mean)'] - 2400000.5
    tin['Observer Affiliation'] = 'MECK'
    tin['Observer Code'] = 'MECK'

    tin['Star Name'] = 'ASASSN-V J060000.76-310027.83'


    tin['Band'][np.where(tin['Band']=='rprime')] = 'rp'
    tin['Band'][np.where(tin['Band']=='iprime')] = 'ip'
    print(tin['Airmass'].dtype)

    t = vstack([t,tin])

t = unique(t, 'MJD')

t.sort('MJD')

# write out the table
t.write(outf, format='ascii.ecsv', overwrite=True)
print(t)
print('wrote out table to {}'.format(outf))


